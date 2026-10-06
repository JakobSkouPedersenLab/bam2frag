use crate::{
    arrays,
    atac::TissueAtacMasks,
    cli::ExtractArgs,
    config, features,
    filters::{load_patient_variant_blacklists, read_passes, FilterContext, IntervalSet},
    model::{ChromStats, ReservoirKey, ValidatedFragmentCandidate},
    output,
    pair_eligibility::{construct_pair_candidate, FragmentObservationStats, PairRejectionCounters},
    reference::ReferenceCache,
    reservoir::{KeyedReservoir, ReservoirStats},
    scan_plan,
};
use anyhow::{bail, Context, Result};
use hashbrown::HashMap;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record};
use std::{collections::BTreeSet, fs, path::Path, time::Instant};
use tracing::info;

struct RegionWork {
    region_id: i64,
    candidates: Vec<ValidatedFragmentCandidate>,
}

fn group_candidates_by_region(
    mut candidates: Vec<ValidatedFragmentCandidate>,
    cfg: &config::ResolvedConfig,
) -> Vec<RegionWork> {
    candidates.sort_by_key(|candidate| {
        (
            candidate.key.region_id,
            candidate.anchor_rpos0,
            candidate.pair.fragment_hash,
        )
    });

    let mut out = Vec::new();
    for candidate in candidates {
        let region_id = candidate.key.region_id;
        let append = out
            .last()
            .is_some_and(|region: &RegionWork| region.region_id == region_id);
        if append {
            out.last_mut().unwrap().candidates.push(candidate);
        } else {
            out.push(RegionWork {
                region_id,
                candidates: vec![candidate],
            });
        }
    }
    out
}

fn materialize_and_pack_region(
    work: RegionWork,
    cfg: &config::ResolvedConfig,
    reference: &mut ReferenceCache,
    atac_masks: Option<&TissueAtacMasks>,
) -> Result<arrays::PackedRegion> {
    let expected_fragments = work.candidates.len();
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    let mut scratch = features::MaterializeScratch::with_width(width);
    let mut builder = arrays::RegionPayloadBuilder::new(work.region_id, expected_fragments, cfg);
    for candidate in work.candidates {
        let scalars = features::materialize_validated_into_scratch(
            candidate,
            cfg,
            reference,
            &mut scratch,
            atac_masks,
        )?;
        builder.push(&scratch, scalars, cfg)?;
    }
    builder.finish(cfg)
}
struct ChromResult {
    candidates: Vec<ValidatedFragmentCandidate>,
    reservoir_stats: Vec<ReservoirStats>,
    stats: ChromStats,
}

/// Apply a fallible transformation in parallel while creating one mutable
/// worker-local state object per Rayon worker. Because the input uses an
/// indexed parallel iterator and the result is collected into a `Vec`, output
/// order matches input order.
///
/// This helper is public only so the concurrency contract can be exercised by
/// an integration test without requiring BAM/FASTA fixtures.
#[doc(hidden)]
pub fn parallel_map_init_ordered<T, S, U, Init, Map>(
    items: Vec<T>,
    init: Init,
    map: Map,
) -> Result<Vec<U>>
where
    T: Send,
    S: Send,
    U: Send,
    Init: Fn() -> Result<S> + Sync + Send,
    Map: Fn(&mut S, T) -> Result<U> + Sync + Send,
{
    items
        .into_par_iter()
        .map_init(init, |state_result, item| {
            let state = state_result
                .as_mut()
                .map_err(|error| anyhow::anyhow!("{error:#}"))?;
            map(state, item)
        })
        .collect()
}

fn need_file(p: &Path, label: &str) -> Result<()> {
    if !p.is_file() {
        bail!("{label} not found: {}", p.display())
    }
    Ok(())
}

fn load_filters(
    args: &ExtractArgs,
    patient_id: &str,
    cfg: &config::ResolvedConfig,
) -> Result<FilterContext> {
    if cfg.requires_snp_blacklist() && args.snp_blacklist.is_none() {
        bail!(
            "configuration requires patient SNP/indel blacklist, but --snp-blacklist was not provided"
        );
    }
    let (snv, indel) = match args.snp_blacklist.as_deref() {
        Some(base) => load_patient_variant_blacklists(base, patient_id)?,
        None => Default::default(),
    };
    info!(
        "patient_id={} snv_blacklist={} indel_blacklist={} lowpower={} problematic={}",
        patient_id,
        snv.len(),
        indel.len(),
        args.lowpower_blacklist.is_some(),
        args.problematic_blacklist.is_some(),
    );
    Ok(FilterContext {
        genomic: args
            .genomic_blacklist
            .as_deref()
            .map(IntervalSet::load)
            .transpose()?,
        donor: args
            .donor_blacklist
            .as_deref()
            .map(IntervalSet::load)
            .transpose()?,
        lowpower: args
            .lowpower_blacklist
            .as_deref()
            .map(IntervalSet::load)
            .transpose()?,
        problematic: args
            .problematic_blacklist
            .as_deref()
            .map(IntervalSet::load)
            .transpose()?,
        snv,
        indel,
        ..Default::default()
    })
}

fn peak_rss_gb() -> Option<f64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmHWM:"))?;
    let kb = line.split_whitespace().nth(1)?.parse::<f64>().ok()?;
    Some(kb / 1024.0 / 1024.0)
}

fn log_stage_done(stage: &str, started: Instant) {
    info!(
        "stage={} elapsed_sec={:.1} peak_rss_gb={:.2}",
        stage,
        started.elapsed().as_secs_f64(),
        peak_rss_gb().unwrap_or(f64::NAN),
    );
}

pub fn run_extract(args: ExtractArgs) -> Result<()> {
    need_file(&args.bam, "BAM")?;
    need_file(&args.reference, "reference FASTA")?;
    let cfg = config::load_resolved(&args.configs)?;
    let patient = args
        .patient_id
        .clone()
        .unwrap_or_else(|| args.sample_id.clone());
    let filters = load_filters(&args, &patient, &cfg)?;
    if args.output_dir.exists() {
        if !args.overwrite {
            bail!(
                "output exists; pass --overwrite: {}",
                args.output_dir.display()
            )
        }
        fs::remove_dir_all(&args.output_dir)?
    }
    let parent = args.output_dir.parent().unwrap_or(Path::new("."));
    let tmp = parent.join(format!(
        "{}.incomplete.{}",
        args.output_dir.file_name().unwrap().to_string_lossy(),
        std::process::id()
    ));
    if tmp.exists() {
        fs::remove_dir_all(&tmp)?
    }
    fs::create_dir_all(&tmp)?;
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.max(1))
        .thread_name(|index| format!("bam2frag-extract-{index}"))
        .build()
        .context("building extraction Rayon thread pool")?;
    info!(
        "parallel_stage=pool_ready requested_threads={} actual_threads={}",
        args.threads.max(1),
        thread_pool.current_num_threads()
    );
    let header = bam::IndexedReader::from_path(&args.bam)?
        .header()
        .to_owned();
    let mut chroms: Vec<(String, u32, u64)> = Vec::new();
    let mut missing_chromosomes = Vec::new();
    for chrom in &cfg.typed.coordinate_system.chromosomes {
        if let Some(tid) = header.tid(chrom.as_bytes()) {
            let bam_len = header.target_len(tid).unwrap_or(0);
            chroms.push((chrom.clone(), tid, bam_len));
        } else {
            missing_chromosomes.push(chrom.clone());
        }
    }
    if !missing_chromosomes.is_empty() {
        bail!(
            "configured chromosomes missing from BAM header: {:?}",
            missing_chromosomes
        );
    }
    let t0 = Instant::now();
    const REGIONS_IN_FLIGHT_MAX: usize = 4;

    // `chunks` panics on zero, and a batch of zero chromosomes would be a
    // silent no-op in any case.
    if args.chromosome_batch_size == 0 {
        bail!("--chromosome-batch-size must be at least 1");
    }
    let chromosome_batch_size = args.chromosome_batch_size;
    if chromosome_batch_size > thread_pool.current_num_threads() {
        // Not an error: a larger batch is sometimes wanted to keep the packer
        // fed. It does cost memory for no extra scan parallelism, though, so
        // say so rather than let it look deliberate in hindsight.
        info!(
            "note: --chromosome-batch-size={} exceeds --threads={}; the surplus chromosomes \
             wait in memory without being scanned any sooner",
            chromosome_batch_size,
            thread_pool.current_num_threads(),
        );
    }

    info!(
        "parallel_stage=chromosome_scan_start chromosomes={} threads={} chromosome_batch_size={} sampling_mode={} stored_fragments_per_region={} peak_rss_gb={:.2}",
        chroms.len(),
        thread_pool.current_num_threads(),
        chromosome_batch_size,
        "natural",
        cfg.typed.reservoir_sampling.stored_fragments_per_region,
        peak_rss_gb().unwrap_or(f64::NAN),
    );

    // Bounded producer/consumer path: one batch of chromosome reservoirs
    // holding BAM records is alive at a time, and each completed batch is
    // materialized, bit-packed and released before the next one is scanned.

    let mut scan_elapsed_sec = 0.0f64;
    let mut materialize_elapsed_sec = 0.0f64;
    let mut rst = Vec::new();
    let mut stats = Vec::new();
    let mut active_keys = BTreeSet::<i64>::new();
    let mut total_retained = 0usize;
    let mut total_regions = 0usize;
    let mut materialize_processed = 0usize;
    let reference_path = args.reference.clone();
    let regions_in_flight = thread_pool
        .current_num_threads()
        .clamp(1, REGIONS_IN_FLIGHT_MAX);
    let mut bank_writer = arrays::FragmentBankWriter::create(&tmp, 0)?;

    info!(
        "annotation_stage=materialize_pack_start rayon_threads={} regions_in_flight={} bounded_chromosome_streaming=true",
        thread_pool.current_num_threads(),
        regions_in_flight,
    );

    for (batch_index, chrom_batch) in chroms.chunks(chromosome_batch_size).enumerate() {
        let batch_started = Instant::now();
        let batch_scan_started = Instant::now();
        let batch_base = batch_index * chromosome_batch_size;
        let results: Vec<Result<ChromResult>> = thread_pool.install(|| {
            chrom_batch
                .par_iter()
                .enumerate()
                .map(|(offset, (chrom, tid, len))| {
                    let chrom_index = batch_base + offset;
                    process_chrom(
                        chrom,
                        *tid,
                        *len,
                        &args,
                        &cfg,
                        &filters,
                        args.seed + 104729 * chrom_index as u64,
                    )
                })
                .collect()
        });
        scan_elapsed_sec += batch_scan_started.elapsed().as_secs_f64();

        for result in results {
            let result = result?;
            let chrom_name = result.stats.chrom.clone();
            let retained_before = result.candidates.len();

            let chrom_keys = arrays::region_keys_from_candidates(&result.candidates, &cfg);
            active_keys.extend(chrom_keys.iter().copied());
            total_regions += chrom_keys.len();
            total_retained += retained_before;
            rst.extend(result.reservoir_stats);
            stats.push(result.stats);

            let chrom_materialize_started = Instant::now();
            let atac_masks = cfg
                .typed
                .fragment_atac_covariates
                .enabled
                .then(|| {
                    let root = args.chromatin_dir.as_deref().context(
                        "--chromatin-dir is required when fragment ATAC covariates are enabled",
                    )?;
                    TissueAtacMasks::load(root, &chrom_name, &cfg)
                })
                .transpose()?;
            let mut region_work = group_candidates_by_region(result.candidates, &cfg);
            while !region_work.is_empty() {
                let take = regions_in_flight.min(region_work.len());
                let batch: Vec<RegionWork> = region_work.drain(..take).collect();
                let packed: Vec<Result<arrays::PackedRegion>> = thread_pool.install(|| {
                    batch
                        .into_par_iter()
                        .map_init(
                            || ReferenceCache::open(&reference_path),
                            |reference_result, work| {
                                let reference = reference_result
                                    .as_mut()
                                    .map_err(|error| anyhow::anyhow!("{error:#}"))?;
                                materialize_and_pack_region(
                                    work,
                                    &cfg,
                                    reference,
                                    atac_masks.as_ref(),
                                )
                            },
                        )
                        .collect()
                });
                for region in packed {
                    let region = region?;
                    materialize_processed += region.fragment_count();
                    bank_writer.write_region(region)?;
                }
            }

            materialize_elapsed_sec += chrom_materialize_started.elapsed().as_secs_f64();
            info!(
                "parallel_stage=chromosome_release chrom={} retained={} materialized={} cumulative_retained={} peak_rss_gb={:.2}",
                chrom_name,
                retained_before,
                materialize_processed,
                total_retained,
                peak_rss_gb().unwrap_or(f64::NAN),
            );
            // `result.candidates` has been moved into region_work and all
            // BAM Records for this chromosome are now dropped.
        }

        info!(
            "parallel_stage=chromosome_batch_done batch={} chromosomes={} elapsed_sec={:.1} peak_rss_gb={:.2}",
            batch_index + 1,
            chrom_batch.len(),
            batch_started.elapsed().as_secs_f64(),
            peak_rss_gb().unwrap_or(f64::NAN),
        );
    }
    info!(
        "stage=chromosome_scan elapsed_sec={:.1} peak_rss_gb={:.2}",
        scan_elapsed_sec,
        peak_rss_gb().unwrap_or(f64::NAN),
    );
    info!(
        "stage=materialize_pack elapsed_sec={:.1} fragments={} rate={:.1}/s peak_rss_gb={:.2}",
        materialize_elapsed_sec,
        materialize_processed,
        materialize_processed as f64 / materialize_elapsed_sec.max(1e-9),
        peak_rss_gb().unwrap_or(f64::NAN),
    );

    let (payload_bytes, stored_bytes) = bank_writer.byte_totals();
    info!(
        "stage=fragment_bank_compression codec={} level={} payload_bytes={} stored_bytes={} ratio={:.2}",
        arrays::BANK_COMPRESSION_CODEC,
        arrays::BANK_COMPRESSION_LEVEL,
        payload_bytes,
        stored_bytes,
        payload_bytes as f64 / (stored_bytes as f64).max(1.0),
    );
    let index_rows = bank_writer.finish(&tmp)?;
    let keys = active_keys.into_iter().collect::<Vec<_>>();
    if index_rows.len() != keys.len() {
        bail!(
            "fragment bank region count {} != active region key count {}",
            index_rows.len(),
            keys.len()
        );
    }

    let summarize_started = Instant::now();
    let am = arrays::write_region_available_counts(&tmp, &index_rows, &rst)?;
    log_stage_done("region_available_counts", summarize_started);

    info!(
        "annotation_stage=materialize_pack_done total={} regions={} rate={:.1}/s elapsed_sec={:.1} peak_rss_gb={:.2}",
        materialize_processed,
        total_regions,
        materialize_processed as f64 / materialize_elapsed_sec.max(1e-9),
        materialize_elapsed_sec,
        peak_rss_gb().unwrap_or(f64::NAN),
    );

    let write_started = Instant::now();
    output::write_stats(
        &tmp.join("extraction_stats.csv"),
        &tmp.join("extraction_stats.json"),
        &stats,
        &args.sample_id,
        &patient,
    )?;
    output::write_meta(
        &tmp.join("sample_meta.json"),
        &args.sample_id,
        &patient,
        &cfg,
        &am,
    )?;
    fs::rename(&tmp, &args.output_dir)?;
    log_stage_done("final_output_write", write_started);

    let (snv_blacklist_hits, indel_blacklist_hits) = filters.blacklist_hit_counts();
    info!(
        "completed sample={} sampling_mode={} regions={} active_regions={} region_output_mode=active_regions_only kept={} snv_blacklist_hits={} indel_blacklist_hits={} elapsed_min={:.2} peak_rss_gb={:.2}",
        args.sample_id,
        "natural",
        am.n_regions,
        am.n_regions,
        am.kept,
        snv_blacklist_hits,
        indel_blacklist_hits,
        t0.elapsed().as_secs_f64() / 60.0,
        peak_rss_gb().unwrap_or(f64::NAN),
    );
    Ok(())
}

fn process_chrom(
    chrom: &str,
    tid: u32,
    chrom_len: u64,
    args: &ExtractArgs,
    cfg: &config::ResolvedConfig,
    filters: &FilterContext,
    seed: u64,
) -> Result<ChromResult> {
    let intervals = scan_plan::scan_intervals(chrom, chrom_len, cfg)?;
    let reservoir_cfg = &cfg.typed.reservoir_sampling;
    let mut res = KeyedReservoir::new(
        reservoir_cfg.stored_fragments_per_region,
        seed,
        reservoir_cfg.max_region_keys,
    );
    let mut stats = ChromStats::new(chrom);
    let mut reference = ReferenceCache::open(&args.reference)?;
    let mut pair_rejections = PairRejectionCounters::default();
    let mut observations = FragmentObservationStats::default();
    let mut bam = bam::IndexedReader::from_path(&args.bam)?;
    let mut processed = 0u64;
    let started = Instant::now();
    for (start, end) in intervals {
        let pad = cfg.typed.reservoir_sampling.fetch_padding_bp as i64;
        let fs = (start as i64 - pad).max(0);
        let fe = end as i64 + pad;
        bam.fetch((tid, fs, fe))?;
        let mut mates: HashMap<Vec<u8>, Record> = HashMap::new();
        let mut order = std::collections::VecDeque::new();
        for record in bam.records() {
            let r = record?;
            if !read_passes(&r, cfg.min_mapq) {
                stats.skipped += 1;
                continue;
            }
            let q = r.qname().to_vec();
            let Some(m) = mates.remove(&q) else {
                while mates.len() >= cfg.max_mate_buffer {
                    if let Some(old) = order.pop_front() {
                        if mates.remove(&old).is_some() {
                            stats.evicted += 1
                        }
                    } else {
                        break;
                    }
                }
                order.push_back(q.clone());
                mates.insert(q, r);
                continue;
            };
            let (a, b) = if r.is_first_in_template() && m.is_last_in_template() {
                (r, m)
            } else if r.is_last_in_template() && m.is_first_in_template() {
                (m, r)
            } else {
                stats.skipped += 1;
                continue;
            };
            stats.reads += 1;
            processed += 1;
            if args.debug_max_fragments > 0 && processed > args.debug_max_fragments {
                break;
            }
            let result =
                construct_pair_candidate(&a, &b, chrom, chrom_len as i64, cfg, &mut reference)?;
            // End-motif base quality is counted, never rejected: the fragment
            // stays eligible for ordinary span reconstruction.
            features::record_pair_decision(
                &mut pair_rejections,
                &mut observations,
                &result,
                &a,
                &b,
            );
            match result {
                Ok(pair) => {
                    let midpoint =
                        pair.left_5p_genomic + (pair.right_5p_genomic - pair.left_5p_genomic) / 2;
                    if midpoint < start as i64 || midpoint >= end as i64 {
                        continue;
                    }
                    if filters.fragment_filtered(
                        chrom,
                        pair.left_5p_genomic,
                        pair.right_5p_genomic + 1,
                    ) {
                        stats.skipped += 1;
                        continue;
                    }
                    let chr_index = cfg
                        .typed
                        .coordinate_system
                        .chromosomes
                        .iter()
                        .position(|name| name == chrom)
                        .context("chromosome missing from configured coordinate system")?
                        as i64;
                    let region_width = cfg.typed.coordinate_system.region_window_bp as i64;
                    let region_id = chr_index * cfg.typed.coordinate_system.region_stride as i64
                        + midpoint.div_euclid(region_width);
                    res.consider(ValidatedFragmentCandidate {
                        key: ReservoirKey { region_id },
                        pair,
                        anchor_rpos0: midpoint,
                    })?;
                }
                Err(_) => stats.skipped += 1,
            }
            if args.progress_every_fragments > 0 && processed % args.progress_every_fragments == 0 {
                info!(
                    "sample={} chrom={} fragments_examined={} rate={:.0}/s elapsed_min={:.1}",
                    args.sample_id,
                    chrom,
                    processed,
                    processed as f64 / started.elapsed().as_secs_f64().max(1e-9),
                    started.elapsed().as_secs_f64() / 60.0
                )
            }
        }
    }
    let (candidates, reservoir_stats) = res.into_parts();
    stats.pair_rejections = pair_rejections.clone();
    stats.observations = observations;
    info!(
        "sample={} chrom={} provisional_rows={} keys={} pair_rejections={:?} elapsed_min={:.2}",
        args.sample_id,
        chrom,
        candidates.len(),
        reservoir_stats.len(),
        pair_rejections,
        started.elapsed().as_secs_f64() / 60.0
    );
    Ok(ChromResult {
        candidates,
        reservoir_stats,
        stats,
    })
}
