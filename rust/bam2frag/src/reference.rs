use anyhow::{Context, Result};
use rust_htslib::faidx;
use std::{collections::HashMap, path::Path};
pub struct ReferenceCache {
    reader: faidx::Reader,
    window: i64,
    cache: HashMap<String, (i64, i64, Vec<u8>)>,
}

impl ReferenceCache {
    pub fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            reader: faidx::Reader::from_path(path)
                .with_context(|| format!("opening FASTA {}", path.display()))?,
            window: 1_000_000,
            cache: HashMap::new(),
        })
    }
    pub fn fetch_span_into(&mut self, chrom: &str, start0: i64, end0: i64, out: &mut Vec<u8>) {
        let len = (end0 - start0).max(0) as usize;
        out.clear();
        out.resize(len, b'N');
        if len == 0 {
            return;
        }
        let query_start = start0.max(0);
        let query_end = end0.max(query_start);
        if query_end <= query_start {
            return;
        }
        let need_reload = self
            .cache
            .get(chrom)
            .map(|(lo, hi, _)| !(*lo <= query_start && query_end <= *hi))
            .unwrap_or(true);
        if need_reload {
            let lo = (query_start - self.window / 2).max(0);
            let hi = lo.saturating_add(self.window.max(query_end - query_start));
            let seq = self
                .reader
                .fetch_seq_string(chrom, lo as usize, (hi - 1) as usize)
                .unwrap_or_default()
                .into_bytes();
            self.cache
                .insert(chrom.into(), (lo, lo + seq.len() as i64, seq));
        }
        if let Some((lo, hi, seq)) = self.cache.get(chrom) {
            let copy_start = query_start.max(*lo);
            let copy_end = query_end.min(*hi);
            if copy_end > copy_start {
                let src = (copy_start - *lo) as usize;
                let dst = (copy_start - start0) as usize;
                let n = (copy_end - copy_start) as usize;
                out[dst..dst + n].copy_from_slice(&seq[src..src + n]);
            }
        }
    }

    pub fn fetch_span(&mut self, chrom: &str, start0: i64, end0: i64) -> Vec<u8> {
        let mut out = Vec::new();
        self.fetch_span_into(chrom, start0, end0, &mut out);
        out
    }

    /// Fetch an exact zero-based, half-open reference interval.
    ///
    /// Unlike the legacy cached helpers, this method never pads an invalid or
    /// truncated request with `N`; callers can distinguish chromosome-edge
    /// rejection from a genuine FASTA/index error before calling it.
    pub fn fetch_span_exact(&mut self, chrom: &str, start0: i64, end0: i64) -> Result<Vec<u8>> {
        if start0 < 0 || end0 <= start0 {
            anyhow::bail!("invalid reference interval {chrom}:{start0}-{end0}");
        }
        let expected_len = (end0 - start0) as usize;
        let sequence = self
            .reader
            .fetch_seq_string(chrom, start0 as usize, (end0 - 1) as usize)
            .with_context(|| format!("fetching exact reference interval {chrom}:{start0}-{end0}"))?
            .into_bytes();
        if sequence.len() != expected_len {
            anyhow::bail!(
                "reference interval {chrom}:{start0}-{end0} returned {} bp, expected {expected_len}",
                sequence.len()
            );
        }
        Ok(sequence)
    }
    pub fn base(&mut self, chrom: &str, pos0: i64) -> u8 {
        self.fetch_span(chrom, pos0, pos0 + 1)
            .first()
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase()
    }
}
