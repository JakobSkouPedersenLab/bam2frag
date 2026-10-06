use crate::config::ResolvedConfig;
use anyhow::{bail, Context, Result};
use ndarray::Array1;
use ndarray_npy::NpzReader;
use std::{fs::File, path::Path};

pub const ATAC_TISSUE_FEATURE_NAMES: [&str; 3] = [
    "atac_blood_overlap_fraction",
    "atac_colon_overlap_fraction",
    "atac_bladder_overlap_fraction",
];

fn load_sparse_coverage_mask(path: &Path) -> Result<Vec<f32>> {
    let mut npz = NpzReader::new(File::open(path)?)?;
    let data: Array1<f32> = npz.by_name("data.npy")?;
    let indices: Array1<i32> = npz.by_name("indices.npy")?;
    let indptr: Array1<i32> = npz.by_name("indptr.npy")?;
    let shape: Array1<i64> = npz.by_name("shape.npy")?;
    if shape.len() != 2 {
        bail!("expected 2D sparse shape in {}", path.display());
    }
    let rows = usize::try_from(shape[0])?;
    let cols = usize::try_from(shape[1])?;
    let mut out = vec![0.0_f32; rows.checked_mul(cols).context("sparse shape overflow")?];
    let is_csr = indptr.len() == rows + 1;
    let is_csc = indptr.len() == cols + 1;
    if !is_csr && !is_csc {
        bail!("cannot infer CSR/CSC for {}", path.display());
    }
    if is_csr {
        for row in 0..rows {
            for k in usize::try_from(indptr[row])?..usize::try_from(indptr[row + 1])? {
                let col = usize::try_from(indices[k])?;
                let value = data[k];
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    bail!(
                        "ATAC coverage value {value} is outside [0,1] in {}",
                        path.display()
                    );
                }
                if col < cols {
                    out[row * cols + col] = value;
                }
            }
        }
    } else {
        for col in 0..cols {
            for k in usize::try_from(indptr[col])?..usize::try_from(indptr[col + 1])? {
                let row = usize::try_from(indices[k])?;
                let value = data[k];
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    bail!(
                        "ATAC coverage value {value} is outside [0,1] in {}",
                        path.display()
                    );
                }
                if row < rows {
                    out[row * cols + col] = value;
                }
            }
        }
    }
    Ok(out)
}

/// Chromosome-specific aggregate tissue ATAC masks.
pub struct TissueAtacMasks {
    /// Fraction of each genomic bin covered by the configured ATAC atlas group.
    masks: [Vec<f32>; 3],
    bin_size_bp: u64,
}

impl TissueAtacMasks {
    /// Load and union configured ATAC atlas members for one chromosome.
    pub fn load(root: &Path, chrom: &str, cfg: &ResolvedConfig) -> Result<Self> {
        let definition = &cfg.typed.fragment_atac_covariates;
        if !definition.enabled {
            bail!("cannot load tissue ATAC masks when fragment_atac_covariates is disabled");
        }
        let mut masks: [Vec<f32>; 3] = std::array::from_fn(|_| Vec::new());
        for (index, group) in definition.groups.iter().enumerate() {
            let mut aggregate = Vec::<f32>::new();
            for member in &group.members {
                let path = root
                    .join("chromatin_masks/ATAC/sparse_overlaps_by_chromosome")
                    .join(member)
                    .join(format!("{chrom}.npz"));
                let member_mask = load_sparse_coverage_mask(&path).with_context(|| {
                    format!(
                        "loading ATAC tissue group {:?} member {:?} for {chrom}",
                        group.name, member
                    )
                })?;
                if aggregate.len() < member_mask.len() {
                    aggregate.resize(member_mask.len(), 0.0);
                }
                for (target, source) in aggregate.iter_mut().zip(member_mask) {
                    // The atlas members describe alternative measurements for the
                    // same tissue group. Their union is the greatest per-bin
                    // coverage, not a sum (which could exceed one).
                    *target = (*target).max(source);
                }
            }
            if aggregate.is_empty() {
                bail!(
                    "ATAC tissue group {:?} produced an empty mask for {chrom}",
                    group.name
                );
            }
            masks[index] = aggregate;
        }
        Ok(Self {
            masks,
            bin_size_bp: definition.mask_bin_size_bp,
        })
    }

    /// Return base-weighted ATAC overlap fractions in blood, colon, bladder order.
    ///
    /// Each atlas value is the covered fraction of a fixed-width genomic bin.
    /// Boundary bins are weighted by the number of fragment bases that fall in
    /// the bin. Because the source only retains bin-level coverage, coverage
    /// inside a partially covered boundary bin is necessarily treated as
    /// uniform within that bin.
    pub fn overlap_fractions(&self, start0: i64, end0: i64) -> [f32; 3] {
        let start0 = start0.max(0) as u64;
        let end0 = end0.max(start0 as i64) as u64;
        let fragment_len = end0.saturating_sub(start0);
        if fragment_len == 0 {
            return [0.0; 3];
        }
        let first_bin = (start0 / self.bin_size_bp) as usize;
        let end_bin = end0.div_ceil(self.bin_size_bp) as usize;
        if end_bin <= first_bin {
            return [0.0; 3];
        }
        std::array::from_fn(|feature_index| {
            let mask = &self.masks[feature_index];
            let stop = end_bin.min(mask.len());
            if stop <= first_bin {
                return 0.0;
            }
            let mut covered_bases = 0.0_f64;
            for (bin_index, &coverage) in mask[first_bin..stop].iter().enumerate() {
                let bin_index = first_bin + bin_index;
                let bin_start = (bin_index as u64) * self.bin_size_bp;
                let bin_end = bin_start + self.bin_size_bp;
                let overlap_start = start0.max(bin_start);
                let overlap_end = end0.min(bin_end);
                let overlap_bp = overlap_end.saturating_sub(overlap_start);
                covered_bases += f64::from(coverage) * overlap_bp as f64;
            }
            (covered_bases / fragment_len as f64) as f32
        })
    }
}

#[cfg(test)]
mod tests {
    use super::TissueAtacMasks;

    #[test]
    fn binned_overlap_fractions_preserve_the_fixed_feature_order() {
        let masks = TissueAtacMasks {
            masks: [
                vec![1.0, 0.0, 1.0, 0.0],
                vec![0.0, 1.0, 1.0, 0.0],
                vec![0.0, 0.0, 0.0, 1.0],
            ],
            bin_size_bp: 10,
        };

        assert_eq!(masks.overlap_fractions(5, 25), [0.5, 0.75, 0.0]);
        assert_eq!(masks.overlap_fractions(30, 40), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn overlap_fractions_preserve_fractional_bin_coverage() {
        let masks = TissueAtacMasks {
            masks: [vec![0.2, 0.8], vec![0.0, 0.0], vec![1.0, 1.0]],
            bin_size_bp: 10,
        };

        assert_eq!(masks.overlap_fractions(0, 20), [0.5, 0.0, 1.0]);
        assert_eq!(masks.overlap_fractions(5, 15), [0.5, 0.0, 1.0]);
        assert_eq!(masks.overlap_fractions(10, 10), [0.0; 3]);
    }
}
