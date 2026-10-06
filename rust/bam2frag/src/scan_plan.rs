use crate::config::ResolvedConfig;
use anyhow::{bail, Result};

/// The genomic spans to scan for one chromosome.
///
/// Regions tile each chromosome contiguously at `region_window_bp` and every
/// region is scanned in full, so there is nothing to plan: the span is the
/// whole chromosome. A fragment is assigned to its region by midpoint as it is
/// read, in `pipeline.rs`.
///
/// This replaced a generic block/window subsampler that could scan a sample of
/// windows per region. Every configuration ran it in whole-region mode, where
/// it built one window per region and then merged those contiguous windows
/// straight back into the single span returned here.
pub fn scan_intervals(
    chrom: &str,
    chrom_len: u64,
    cfg: &ResolvedConfig,
) -> Result<Vec<(u64, u64)>> {
    let cs = &cfg.typed.coordinate_system;
    let n_regions = chrom_len.div_ceil(cs.region_window_bp);

    // `region_stride` is the width of each chromosome's slice of the region-id
    // space, not a genomic distance: ids are `chr_idx * region_stride +
    // region_bin`. If one chromosome ever produces more bins than that, its
    // high bins collide with the next chromosome's low ones and two
    // chromosomes' fragments silently share a reservoir. At the shipped 1 Mb
    // region width the margin is ~4000x, but a smaller width closes it fast:
    // below about 249 bp, chr1 alone overflows a stride of 1,000,000.
    if n_regions > cs.region_stride {
        bail!(
            "{chrom} yields {n_regions} regions at region_window_bp={}, which exceeds \
             region_stride={}; region ids would collide with the next chromosome",
            cs.region_window_bp,
            cs.region_stride,
        );
    }

    Ok(vec![(0, chrom_len)])
}
