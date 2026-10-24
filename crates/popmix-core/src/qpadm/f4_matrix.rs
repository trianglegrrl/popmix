//! F4 matrix computation for qpWave/qpAdm.
//!
//! The F4 matrix is computed between left populations L = {L1, ..., Ln}
//! and right populations R = {R1, ..., Rm}:
//!
//! F4[i,j] = f4(Li+1, L1; Rj+1, R1) for i in 0..n-1, j in 0..m-1
//!
//! This matrix has dimensions (n-1) x (m-1) because L1 and R1 serve as
//! reference populations.
//!
//! The covariance matrix is the weighted block jackknife of admixtools2
//! (`est_to_loo_nafix` + `jack_pairarr_stats`), with blocks weighted by
//! their length in SNPs; see [`block_jackknife_cov`].

use ndarray::{Array1, Array2};

use crate::fstats::{block_jackknife_cov, leave_one_out_nafix, F2Blocks};
use crate::genotypes::PopulationId;
use crate::Result;

/// F4 matrix with associated covariance for statistical tests.
#[derive(Debug, Clone)]
pub struct F4Matrix {
    /// The F4 values matrix, shape (n_left-1, n_right-1)
    pub values: Array2<f64>,
    /// Standard errors for each F4 value
    pub std_errors: Array2<f64>,
    /// Covariance matrix of the vectorized F4 values
    /// Shape: ((n_left-1)*(n_right-1), (n_left-1)*(n_right-1))
    pub covariance: Array2<f64>,
    /// Leave-one-block-out f4 matrices, one per block (admixtools2
    /// `est_to_loo_nafix`), used for qpAdm weight standard errors.
    pub leave_one_out: Vec<Array2<f64>>,
    /// Left populations (first is reference)
    pub left_pops: Vec<PopulationId>,
    /// Right populations (first is reference)
    pub right_pops: Vec<PopulationId>,
}

impl F4Matrix {
    /// Compute F4 matrix from F2 blocks.
    ///
    /// # Arguments
    /// * `f2_blocks` - Precomputed F2 block statistics
    /// * `left` - Left population IDs (first serves as reference)
    /// * `right` - Right population IDs (first serves as reference)
    ///
    /// # Returns
    /// F4Matrix with values, standard errors, and covariance matrix
    pub fn from_f2_blocks(f2_blocks: &F2Blocks, left: &[&str], right: &[&str]) -> Result<Self> {
        if left.len() < 2 {
            return Err(crate::PopmixError::Validation(
                "qpWave/qpAdm requires at least 2 left populations".to_string(),
            ));
        }
        if right.len() < 2 {
            return Err(crate::PopmixError::Validation(
                "qpWave/qpAdm requires at least 2 right populations".to_string(),
            ));
        }

        let n_left = left.len();
        let n_right = right.len();
        let n_rows = n_left - 1;
        let n_cols = n_right - 1;

        // Compute F4 matrix values and standard errors
        let mut values = Array2::zeros((n_rows, n_cols));
        let mut std_errors = Array2::zeros((n_rows, n_cols));

        for i in 0..n_rows {
            for j in 0..n_cols {
                // f4(Li+1, L0; Rj+1, R0)
                if let Some(result) = f2_blocks.f4(left[i + 1], left[0], right[j + 1], right[0]) {
                    values[[i, j]] = result.estimate;
                    std_errors[[i, j]] = result.std_err;
                } else {
                    return Err(crate::PopmixError::Validation(format!(
                        "Missing data for f4({}, {}; {}, {})",
                        left[i + 1],
                        left[0],
                        right[j + 1],
                        right[0]
                    )));
                }
            }
        }

        // Compute full covariance matrix using jackknife
        let (covariance, leave_one_out) = Self::compute_covariance(f2_blocks, left, right)?;

        Ok(Self {
            values,
            std_errors,
            covariance,
            leave_one_out,
            left_pops: left.iter().map(|s| PopulationId::from(*s)).collect(),
            right_pops: right.iter().map(|s| PopulationId::from(*s)).collect(),
        })
    }

    /// Compute the covariance matrix of all F4 statistics using weighted jackknife.
    ///
    /// Uses the admixtools2 formula which accounts for variable block sizes:
    /// h_k = total_snps / block_snps_k
    /// xtau_k = (est * h_k - loo_k * (h_k - 1) - block_k) / sqrt(h_k - 1)
    /// var = sum(xtau * xtau') / numblocks
    ///
    /// The covariance is needed for proper chi-squared tests in qpWave.
    fn compute_covariance(
        f2_blocks: &F2Blocks,
        left: &[&str],
        right: &[&str],
    ) -> Result<(Array2<f64>, Vec<Array2<f64>>)> {
        let n_rows = left.len() - 1;
        let n_cols = right.len() - 1;
        let n_stats = n_rows * n_cols;
        let n_blocks = f2_blocks.n_blocks;

        if n_blocks < 2 {
            // Not enough blocks for jackknife
            return Ok((Array2::zeros((n_stats, n_stats)), Vec::new()));
        }

        // Compute block-level F4 values for each statistic
        // f4(A,B;C,D) = (f2(A,D) + f2(B,C) - f2(A,C) - f2(B,D)) / 2
        let mut block_f4s: Vec<Array1<f64>> = Vec::with_capacity(n_stats);

        for i in 0..n_rows {
            for j in 0..n_cols {
                let pop_a = left[i + 1];
                let pop_b = left[0];
                let pop_c = right[j + 1];
                let pop_d = right[0];

                // Get block-level f2 values
                let f2_ad = f2_blocks.get_block_f2_values(pop_a, pop_d);
                let f2_bc = f2_blocks.get_block_f2_values(pop_b, pop_c);
                let f2_ac = f2_blocks.get_block_f2_values(pop_a, pop_c);
                let f2_bd = f2_blocks.get_block_f2_values(pop_b, pop_d);

                match (f2_ad, f2_bc, f2_ac, f2_bd) {
                    (Some(ad), Some(bc), Some(ac), Some(bd)) => {
                        let f4_blocks: Array1<f64> = (&ad + &bc - &ac - &bd) / 2.0;
                        block_f4s.push(f4_blocks);
                    }
                    _ => {
                        return Err(crate::PopmixError::Validation(format!(
                            "Missing block data for f4({}, {}; {}, {})",
                            pop_a, pop_b, pop_c, pop_d
                        )));
                    }
                }
            }
        }

        let stats: Vec<Vec<f64>> = block_f4s.iter().map(|b| b.to_vec()).collect();
        let lengths: Vec<f64> = f2_blocks.block_lengths.iter().map(|&l| l as f64).collect();
        let covariance = match block_jackknife_cov(&stats, &lengths) {
            Some((_, cov)) => cov,
            None => Array2::zeros((n_stats, n_stats)), // Not enough valid blocks
        };

        // Leave-one-out f4 matrices, stats in row-major order
        let loo: Vec<Vec<f64>> = stats
            .iter()
            .map(|x| leave_one_out_nafix(x, &lengths))
            .collect();
        let leave_one_out = (0..n_blocks)
            .map(|k| Array2::from_shape_fn((n_rows, n_cols), |(i, j)| loo[i * n_cols + j][k]))
            .collect();

        Ok((covariance, leave_one_out))
    }

    /// Get the number of rows (n_left - 1).
    pub fn n_rows(&self) -> usize {
        self.values.nrows()
    }

    /// Get the number of columns (n_right - 1).
    pub fn n_cols(&self) -> usize {
        self.values.ncols()
    }

    /// Flatten the F4 matrix to a vector (row-major order).
    pub fn to_vector(&self) -> Array1<f64> {
        Array1::from_iter(self.values.iter().copied())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use ndarray::Array2;

    fn create_test_f2_blocks() -> F2Blocks {
        // Create synthetic allele frequencies for 4 populations
        let pops = vec!["Pop1", "Pop2", "Pop3", "Pop4"];
        let n_snps = 1000;
        let n_blocks = 10;

        // Create frequencies with known relationships
        // Pop1 and Pop2 are closer, Pop3 and Pop4 are outgroups
        let mut rng_state: u64 = 12345;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 4));
        let mut counts = Array2::from_elem((n_snps, 4), 50u32);

        for i in 0..n_snps {
            // Base frequency
            let base = rand() * 0.8 + 0.1;

            // Pop1 and Pop2 share ancestry (small drift)
            freqs[[i, 0]] = (base + rand() * 0.05).clamp(0.0, 1.0);
            freqs[[i, 1]] = (base + rand() * 0.05).clamp(0.0, 1.0);

            // Pop3 is more divergent
            freqs[[i, 2]] = (base + rand() * 0.2 - 0.1).clamp(0.0, 1.0);

            // Pop4 is most divergent (outgroup)
            freqs[[i, 3]] = rand();
        }

        // Create chromosomes and positions for block assignment
        let mut chromosomes = Vec::with_capacity(n_snps);
        let mut positions = Vec::with_capacity(n_snps);

        for i in 0..n_snps {
            chromosomes.push(crate::genotypes::Chromosome::Autosome(1));
            // Spread SNPs across 0.5 Morgan (10 blocks of 0.05 Morgan)
            positions.push(i as f64 * 0.0005);
        }

        let af = AlleleFrequencies::new(
            freqs,
            counts,
            pops.iter()
                .map(|s| crate::genotypes::PopulationId::from(*s))
                .collect(),
            chromosomes,
            positions,
        );

        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_f4_matrix_dimensions() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2"];
        let right = vec!["Pop3", "Pop4"];

        let f4_matrix = F4Matrix::from_f2_blocks(&f2_blocks, &left, &right).unwrap();

        // (2-1) x (2-1) = 1 x 1
        assert_eq!(f4_matrix.n_rows(), 1);
        assert_eq!(f4_matrix.n_cols(), 1);
        assert_eq!(f4_matrix.covariance.shape(), &[1, 1]);
    }

    #[test]
    fn test_f4_matrix_larger() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2", "Pop3"];
        let right = vec!["Pop3", "Pop4"];

        let f4_matrix = F4Matrix::from_f2_blocks(&f2_blocks, &left, &right).unwrap();

        // (3-1) x (2-1) = 2 x 1
        assert_eq!(f4_matrix.n_rows(), 2);
        assert_eq!(f4_matrix.n_cols(), 1);
        assert_eq!(f4_matrix.covariance.shape(), &[2, 2]);
    }

    #[test]
    fn test_f4_matrix_requires_min_pops() {
        let f2_blocks = create_test_f2_blocks();

        // Only 1 left population should fail
        let result = F4Matrix::from_f2_blocks(&f2_blocks, &["Pop1"], &["Pop3", "Pop4"]);
        assert!(result.is_err());

        // Only 1 right population should fail
        let result = F4Matrix::from_f2_blocks(&f2_blocks, &["Pop1", "Pop2"], &["Pop3"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_f4_values_match_direct_computation() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2"];
        let right = vec!["Pop3", "Pop4"];

        let f4_matrix = F4Matrix::from_f2_blocks(&f2_blocks, &left, &right).unwrap();

        // Compare with direct f4 computation
        let direct_f4 = f2_blocks.f4("Pop2", "Pop1", "Pop4", "Pop3").unwrap();

        assert!(
            (f4_matrix.values[[0, 0]] - direct_f4.estimate).abs() < 1e-10,
            "F4 matrix value {} doesn't match direct computation {}",
            f4_matrix.values[[0, 0]],
            direct_f4.estimate
        );
    }

    #[test]
    fn test_covariance_is_symmetric() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2", "Pop3"];
        let right = vec!["Pop3", "Pop4"];

        let f4_matrix = F4Matrix::from_f2_blocks(&f2_blocks, &left, &right).unwrap();

        let cov = &f4_matrix.covariance;
        for i in 0..cov.nrows() {
            for j in 0..cov.ncols() {
                assert!(
                    (cov[[i, j]] - cov[[j, i]]).abs() < 1e-15,
                    "Covariance not symmetric at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_covariance_diagonal_positive() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2", "Pop3"];
        let right = vec!["Pop3", "Pop4"];

        let f4_matrix = F4Matrix::from_f2_blocks(&f2_blocks, &left, &right).unwrap();

        for i in 0..f4_matrix.covariance.nrows() {
            assert!(
                f4_matrix.covariance[[i, i]] >= 0.0,
                "Variance at {} is negative: {}",
                i,
                f4_matrix.covariance[[i, i]]
            );
        }
    }
}
