use anyhow::{Context, Result};
use rust_htslib::bam::Record;
use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    },
};
#[derive(Clone, Default)]
pub struct IntervalSet {
    by_chrom: BTreeMap<String, Vec<(i64, i64)>>,
}

impl IntervalSet {
    pub fn load(path: &Path) -> Result<Self> {
        let f = BufReader::new(
            File::open(path).with_context(|| format!("opening {}", path.display()))?,
        );
        let mut s = Self::default();
        for line in f.lines() {
            let l = line?;
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let p: Vec<_> = l.split('\t').collect();
            if p.len() < 3 {
                continue;
            }
            let start = p[1].parse()?;
            let end = p[2].parse()?;
            s.by_chrom
                .entry(p[0].into())
                .or_default()
                .push((start, end));
        }
        for v in s.by_chrom.values_mut() {
            v.sort_unstable()
        }
        Ok(s)
    }
    pub fn contains(&self, c: &str, p: i64) -> bool {
        let Some(v) = self.by_chrom.get(c) else {
            return false;
        };
        let i = v.partition_point(|(s, _)| *s <= p);
        i > 0 && p < v[i - 1].1
    }

    pub fn overlaps(&self, c: &str, start0: i64, end0: i64) -> bool {
        self.by_chrom.get(c).is_some_and(|intervals| {
            intervals
                .iter()
                .any(|(start, end)| *start < end0 && *end > start0)
        })
    }
}

#[derive(Clone, Default)]
pub struct FilterContext {
    pub genomic: Option<IntervalSet>,
    pub donor: Option<IntervalSet>,
    pub lowpower: Option<IntervalSet>,
    pub problematic: Option<IntervalSet>,
    pub snv: HashSet<(String, i64, String, String)>,
    pub indel: HashSet<(String, i64, i64, String, String)>,
    pub(crate) snv_blacklist_hits: Arc<AtomicU64>,
    pub(crate) indel_blacklist_hits: Arc<AtomicU64>,
}

impl FilterContext {
    pub fn fragment_filtered(&self, chrom: &str, start0: i64, end0: i64) -> bool {
        [
            &self.genomic,
            &self.donor,
            &self.lowpower,
            &self.problematic,
        ]
        .into_iter()
        .flatten()
        .any(|intervals| intervals.overlaps(chrom, start0, end0))
    }

    pub fn position_filtered(&self, c: &str, p: i64) -> bool {
        self.genomic.as_ref().is_some_and(|x| x.contains(c, p))
            || self.donor.as_ref().is_some_and(|x| x.contains(c, p))
            || self.lowpower.as_ref().is_some_and(|x| x.contains(c, p))
            || self.problematic.as_ref().is_some_and(|x| x.contains(c, p))
    }

    pub fn snv_blacklisted(&self, chrom: &str, pos0: i64, obs_ref: &str, obs_alt: &str) -> bool {
        let hit = self.snv.contains(&(
            chrom.to_string(),
            pos0,
            obs_ref.to_ascii_uppercase(),
            obs_alt.to_ascii_uppercase(),
        ));
        if hit {
            self.snv_blacklist_hits.fetch_add(1, Ordering::Relaxed);
        }
        hit
    }

    pub fn indel_blacklisted(
        &self,
        chrom: &str,
        start0: i64,
        end0: i64,
        obs_ref: &str,
        obs_alt: &str,
    ) -> bool {
        let hit = self.indel.contains(&(
            chrom.to_string(),
            start0,
            end0,
            obs_ref.to_ascii_uppercase(),
            obs_alt.to_ascii_uppercase(),
        ));
        if hit {
            self.indel_blacklist_hits.fetch_add(1, Ordering::Relaxed);
        }
        hit
    }

    pub fn blacklist_hit_counts(&self) -> (u64, u64) {
        (
            self.snv_blacklist_hits.load(Ordering::Relaxed),
            self.indel_blacklist_hits.load(Ordering::Relaxed),
        )
    }
}

fn header_index(headers: &csv::StringRecord, aliases: &[&str]) -> Option<usize> {
    headers.iter().position(|name| {
        // Variant blacklist headers use names such as "#CHR". Strip a single
        // leading comment marker from header fields, but do not configure the
        // CSV reader to skip '#' lines because the header itself starts with '#'.
        let normalized = name
            .trim()
            .trim_start_matches('#')
            .trim()
            .to_ascii_lowercase();
        aliases.iter().any(|alias| normalized == *alias)
    })
}

fn required_header_index(
    headers: &csv::StringRecord,
    aliases: &[&str],
    label: &str,
    path: &Path,
) -> Result<usize> {
    header_index(headers, aliases).ok_or_else(|| {
        anyhow::anyhow!(
            "{} is missing required {} column; accepted names={:?}",
            path.display(),
            label,
            aliases
        )
    })
}

fn normalize_chrom(chrom: &str) -> String {
    let chrom = chrom.trim();
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{chrom}")
    }
}

fn resolve_variant_file(base: &Path, patient_id: &str, suffix: &str) -> Result<PathBuf> {
    let patient_dir = base.join(patient_id);
    let exact = patient_dir.join(format!("{patient_id}_{suffix}_calls.tsv"));
    if exact.is_file() {
        return Ok(exact);
    }

    let patient_lower = patient_id.to_ascii_lowercase();
    let mut patient_matches = std::fs::read_dir(base)
        .with_context(|| format!("reading SNP blacklist base directory {}", base.display()))?
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.path().is_dir())
        .filter(|entry| entry.file_name().to_string_lossy().to_ascii_lowercase() == patient_lower)
        .collect::<Vec<_>>();

    if patient_matches.len() != 1 {
        anyhow::bail!(
            "expected exactly one patient directory for patient_id={patient_id:?} under {}; found {}",
            base.display(),
            patient_matches.len()
        );
    }

    let dir = patient_matches.pop().unwrap().path();
    let expected_lower = format!("{}_{}_calls.tsv", patient_id, suffix).to_ascii_lowercase();
    let mut matches = std::fs::read_dir(&dir)
        .with_context(|| format!("reading patient SNP blacklist directory {}", dir.display()))?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| path.is_file())
        .filter(|path| {
            path.file_name()
                .map(|name| name.to_string_lossy().to_ascii_lowercase() == expected_lower)
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();

    if matches.len() != 1 {
        anyhow::bail!(
            "expected exactly one {suffix} blacklist file named {patient_id}_{suffix}_calls.tsv under {}; found {}",
            dir.display(),
            matches.len()
        );
    }
    Ok(matches.pop().unwrap())
}

fn parse_pos0(record: &csv::StringRecord, headers: &csv::StringRecord, path: &Path) -> Result<i64> {
    // BED-style blacklist files use zero-based START coordinates.
    if let Some(index) = header_index(headers, &["start", "start0", "pos0", "position0"]) {
        return record[index]
            .trim()
            .parse::<i64>()
            .with_context(|| format!("parsing 0-based position in {}", path.display()));
    }

    // Generic POS/POSITION columns are interpreted as one-based, matching VCF-
    // style conventions. Explicit position_1based/pos1 names are also accepted.
    let index = required_header_index(
        headers,
        &["pos", "position", "position_1based", "pos1"],
        "position",
        path,
    )?;
    let pos1 = record[index]
        .trim()
        .parse::<i64>()
        .with_context(|| format!("parsing 1-based position in {}", path.display()))?;
    if pos1 <= 0 {
        anyhow::bail!(
            "1-based position must be positive in {}: {pos1}",
            path.display()
        );
    }
    Ok(pos1 - 1)
}

fn load_snv_blacklist(path: &Path) -> Result<HashSet<(String, i64, String, String)>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening SNV blacklist {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let chrom_i = required_header_index(
        &headers,
        &["chrom", "chromosome", "chr"],
        "chromosome",
        path,
    )?;
    let ref_i = required_header_index(
        &headers,
        &["ref", "reference", "reference_allele"],
        "reference allele",
        path,
    )?;
    let alt_i = required_header_index(
        &headers,
        &["alt", "alternate", "alternative", "alt_allele"],
        "alternate allele",
        path,
    )?;
    let mut out = HashSet::new();
    for row in reader.records() {
        let row = row?;
        let chrom = normalize_chrom(&row[chrom_i]);
        let pos0 = parse_pos0(&row, &headers, path)?;
        let obs_ref = row[ref_i].trim().to_ascii_uppercase();
        for alt in row[alt_i].split(',') {
            let obs_alt = alt.trim().to_ascii_uppercase();
            if obs_ref.len() == 1 && obs_alt.len() == 1 && obs_ref != obs_alt {
                out.insert((chrom.clone(), pos0, obs_ref.clone(), obs_alt));
            }
        }
    }
    Ok(out)
}

fn load_indel_blacklist(path: &Path) -> Result<HashSet<(String, i64, i64, String, String)>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening indel blacklist {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let chrom_i = required_header_index(
        &headers,
        &["chrom", "chromosome", "chr"],
        "chromosome",
        path,
    )?;
    let ref_i = required_header_index(
        &headers,
        &["ref", "reference", "reference_allele"],
        "reference allele",
        path,
    )?;
    let alt_i = required_header_index(
        &headers,
        &["alt", "alternate", "alternative", "alt_allele"],
        "alternate allele",
        path,
    )?;
    let end0_i = header_index(&headers, &["end0"]);
    let end_i = header_index(&headers, &["end", "stop"]);
    let mut out = HashSet::new();
    for row in reader.records() {
        let row = row?;
        let chrom = normalize_chrom(&row[chrom_i]);
        let start0 = parse_pos0(&row, &headers, path)?;
        let obs_ref = row[ref_i].trim().to_ascii_uppercase();
        for alt in row[alt_i].split(',') {
            let obs_alt = alt.trim().to_ascii_uppercase();
            if obs_ref == obs_alt || obs_ref.is_empty() || obs_alt.is_empty() {
                continue;
            }
            let end0 = if let Some(index) = end0_i {
                row[index]
                    .trim()
                    .parse::<i64>()
                    .with_context(|| format!("parsing end0 in {}", path.display()))?
            } else if let Some(index) = end_i {
                row[index]
                    .trim()
                    .parse::<i64>()
                    .with_context(|| format!("parsing end in {}", path.display()))?
            } else {
                start0 + obs_ref.len() as i64
            };
            out.insert((chrom.clone(), start0, end0, obs_ref.clone(), obs_alt));
        }
    }
    Ok(out)
}

pub fn load_patient_variant_blacklists(
    base: &Path,
    patient_id: &str,
) -> Result<(
    HashSet<(String, i64, String, String)>,
    HashSet<(String, i64, i64, String, String)>,
)> {
    if !base.is_dir() {
        anyhow::bail!("SNP blacklist base directory not found: {}", base.display());
    }
    let snv_path = resolve_variant_file(base, patient_id, "snvs")?;
    let indel_path = resolve_variant_file(base, patient_id, "indels")?;
    let snv = load_snv_blacklist(&snv_path)?;
    let indel = load_indel_blacklist(&indel_path)?;
    if snv.is_empty() && indel.is_empty() {
        anyhow::bail!(
            "SNP blacklist directory was provided but no usable variants were loaded for patient_id={patient_id:?}"
        );
    }
    Ok((snv, indel))
}

pub fn read_passes(r: &Record, min_mapq: u8) -> bool {
    !(r.is_unmapped()
        || r.is_secondary()
        || r.is_supplementary()
        || r.is_duplicate()
        || r.is_quality_check_failed()
        || !r.is_paired())
        && r.mapq() >= min_mapq
        && !r.qual().is_empty()
}
