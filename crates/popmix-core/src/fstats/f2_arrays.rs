//! F2 blocks as their arrays: read through accessors, and built from arrays
//! (an f2 directory, popmix-io's `f2_storage`) with their shapes checked.

use ndarray::Array2;

use super::f2_blocks::{F2Blocks, FstBlocks};
use crate::error::{PopmixError, Result};
use crate::genotypes::PopulationId;

/// The arrays of an [`F2Blocks`], for [`F2Blocks::from_arrays`]. Pairs are
/// in [`F2Blocks::pair_index`] order: (0,0), (0,1), ..., (1,1), ...
#[derive(Debug, Clone)]
pub struct F2Arrays {
    /// Bias-corrected f2 sums: [blocks, pairs].
    pub block_sums: Array2<f64>,
    /// Squared frequency differences, not bias-corrected: [blocks, pairs].
    pub block_raw_sums: Array2<f64>,
    /// SNPs with data, per block and pair: [blocks, pairs].
    pub block_weights: Array2<f64>,
    /// Heterozygosity sums, p(1 - p): [blocks, populations].
    pub block_het_sums: Array2<f64>,
    /// SNPs with data, per block and population: [blocks, populations].
    pub block_het_weights: Array2<f64>,
    pub populations: Vec<PopulationId>,
    /// SNPs per block.
    pub block_lengths: Vec<usize>,
    /// FST on its own SNP set and blocks, if known.
    pub fst: Option<FstBlocks>,
}

fn check_shape(what: &str, arr: &Array2<f64>, want: (usize, usize)) -> Result<()> {
    if arr.dim() == want {
        return Ok(());
    }
    Err(PopmixError::DimensionMismatch {
        context: what.to_string(),
        expected: format!("{want:?}"),
        actual: format!("{:?}", arr.dim()),
    })
}

impl F2Blocks {
    /// F2 blocks from their arrays, which must agree in shape.
    pub fn from_arrays(a: F2Arrays) -> Result<Self> {
        let n = a.populations.len();
        let nb = a.block_lengths.len();
        let pairs = (nb, n * (n + 1) / 2);
        check_shape("block_sums", &a.block_sums, pairs)?;
        check_shape("block_raw_sums", &a.block_raw_sums, pairs)?;
        check_shape("block_weights", &a.block_weights, pairs)?;
        check_shape("block_het_sums", &a.block_het_sums, (nb, n))?;
        check_shape("block_het_weights", &a.block_het_weights, (nb, n))?;
        if let Some(fst) = &a.fst {
            check_shape("FST blocks", &fst.num, (fst.block_lengths.len(), pairs.1))?;
        }
        Ok(F2Blocks {
            block_sums: a.block_sums,
            block_raw_sums: a.block_raw_sums,
            block_weights: a.block_weights,
            block_het_sums: a.block_het_sums,
            block_het_weights: a.block_het_weights,
            populations: a.populations,
            n_blocks: nb,
            block_lengths: a.block_lengths,
            fst: a.fst,
            provenance: None,
        })
    }

    /// Bias-corrected f2 sums: [blocks, pairs].
    pub fn block_sums(&self) -> &Array2<f64> {
        &self.block_sums
    }

    /// Squared frequency differences, not bias-corrected: [blocks, pairs].
    pub fn block_raw_sums(&self) -> &Array2<f64> {
        &self.block_raw_sums
    }

    /// SNPs with data, per block and pair: [blocks, pairs].
    pub fn block_weights(&self) -> &Array2<f64> {
        &self.block_weights
    }

    /// Heterozygosity sums, p(1 - p): [blocks, populations].
    pub fn block_het_sums(&self) -> &Array2<f64> {
        &self.block_het_sums
    }

    /// SNPs with data, per block and population: [blocks, populations].
    pub fn block_het_weights(&self) -> &Array2<f64> {
        &self.block_het_weights
    }

    /// The populations, in the order of the arrays.
    pub fn populations(&self) -> &[PopulationId] {
        &self.populations
    }

    pub fn n_blocks(&self) -> usize {
        self.n_blocks
    }

    /// SNPs per block.
    pub fn block_lengths(&self) -> &[usize] {
        &self.block_lengths
    }

    /// FST on its own SNP set and blocks, or None if not known (f2
    /// directories written before storage version 4).
    pub fn fst_blocks(&self) -> Option<&FstBlocks> {
        self.fst.as_ref()
    }
}

impl FstBlocks {
    /// FST blocks from their block lengths and their numerator and
    /// denominator block means, [blocks, pairs] each.
    pub fn new(block_lengths: Vec<usize>, num: Array2<f64>, den: Array2<f64>) -> Result<Self> {
        if num.nrows() != block_lengths.len() {
            return Err(PopmixError::DimensionMismatch {
                context: "FST numerator blocks".to_string(),
                expected: block_lengths.len().to_string(),
                actual: num.nrows().to_string(),
            });
        }
        check_shape("FST denominator", &den, num.dim())?;
        Ok(FstBlocks {
            block_lengths,
            num,
            den,
        })
    }

    /// SNPs per block.
    pub fn block_lengths(&self) -> &[usize] {
        &self.block_lengths
    }

    /// Block means of the numerator: [blocks, pairs].
    pub fn num(&self) -> &Array2<f64> {
        &self.num
    }

    /// Block means of the denominator: [blocks, pairs].
    pub fn den(&self) -> &Array2<f64> {
        &self.den
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arrays(n_pops: usize, n_blocks: usize) -> F2Arrays {
        let pairs = n_pops * (n_pops + 1) / 2;
        let a = |cols| Array2::from_shape_fn((n_blocks, cols), |(b, j)| (b * 10 + j) as f64);
        F2Arrays {
            block_sums: a(pairs),
            block_raw_sums: a(pairs),
            block_weights: a(pairs),
            block_het_sums: a(n_pops),
            block_het_weights: a(n_pops),
            populations: (0..n_pops)
                .map(|i| PopulationId::from(format!("P{i}").as_str()))
                .collect(),
            block_lengths: vec![100; n_blocks],
            fst: None,
        }
    }

    #[test]
    fn f2_blocks_hold_the_arrays_they_are_built_from() {
        let f2 = F2Blocks::from_arrays(arrays(3, 4)).unwrap();
        assert_eq!(f2.n_blocks(), 4);
        assert_eq!(f2.populations().len(), 3);
        assert_eq!(f2.block_sums()[[2, 5]], 25.0);
        assert_eq!(f2.block_het_sums().dim(), (4, 3));
        assert!(f2.fst_blocks().is_none());
    }

    #[test]
    fn arrays_of_the_wrong_shape_are_an_error() {
        let mut a = arrays(3, 4);
        a.block_weights = Array2::zeros((4, 5));
        let err = F2Blocks::from_arrays(a).unwrap_err();
        assert!(err.to_string().contains("block_weights"), "{err}");
        let err =
            FstBlocks::new(vec![1, 2], Array2::zeros((3, 6)), Array2::zeros((3, 6))).unwrap_err();
        assert!(
            matches!(err, PopmixError::DimensionMismatch { .. }),
            "{err}"
        );
    }
}
