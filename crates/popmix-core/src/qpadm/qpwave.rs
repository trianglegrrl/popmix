//! qpWave rank testing implementation.
//!
//! qpWave tests whether a set of "left" populations can be modeled as
//! mixtures of a common set of ancestral sources, using "right" populations
//! as outgroups.
//!
//! The test works by computing an F4 matrix and testing whether its rank
//! is less than expected under a null hypothesis of shared ancestry. Rank
//! tests follow admixtools2 `qpwave` / `drop_ranks`; see [`super::fit`].
//!
//! The estimated rank r (r + 1 streams of ancestry) is a lower bound: the
//! f4 matrix has rank at most min(left, right) − 1, right populations that
//! do not differ in their relation to the left ones hide streams, and a
//! test without the power to reject reads as "not rejected". Not rejecting
//! rank 0 does not show that the left populations form a clade.

use super::{fit, F4Matrix};
use crate::fstats::F2Blocks;
use crate::genotypes::check_populations;
use crate::Result;

/// Configuration options for qpWave analysis.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpWaveOptions {
    /// Regularization factor for numerical stability (default: 1e-4).
    pub fudge: f64,
    /// Highest rank to test. Default: number of left populations - 2, as
    /// admixtools2 `drop_ranks`.
    pub max_rank: Option<usize>,
    /// Significance level of the rank tests: the estimated rank is the
    /// lowest one whose p exceeds it (default 0.05).
    pub alpha: f64,
}

impl Default for QpWaveOptions {
    fn default() -> Self {
        Self {
            fudge: 1e-4,
            max_rank: None,
            alpha: 0.05,
        }
    }
}

/// Result of a single rank test.
#[derive(Debug, Clone)]
pub struct RankTest {
    /// The rank being tested (null hypothesis: rank <= this value).
    pub rank: usize,
    /// Chi-squared statistic.
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
    /// P-value for the test.
    pub p_value: f64,
}

/// Result of qpWave analysis.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpWaveResult {
    /// Series of rank tests from 0 to max possible.
    pub rank_tests: Vec<RankTest>,
    /// The lowest rank not rejected at `alpha` (one above the highest rank
    /// tested if every test is rejected): at least this many plus one
    /// streams of ancestry, as seen from the right populations.
    pub estimated_rank: usize,
    /// The F4 matrix used for analysis.
    pub f4_matrix: F4Matrix,
}

impl QpWaveResult {
    /// Get the full rank test result.
    pub fn full_rank(&self) -> Option<&RankTest> {
        self.rank_tests.last()
    }
}

/// Run qpWave analysis.
///
/// # Arguments
/// * `f2_blocks` - Precomputed F2 statistics
/// * `left` - Left population names
/// * `right` - Right population names (outgroups)
/// * `options` - Optional analysis parameters
///
/// # Returns
/// QpWaveResult containing rank test results
pub fn qpwave(
    f2_blocks: &F2Blocks,
    left: &[&str],
    right: &[&str],
    options: &QpWaveOptions,
) -> Result<QpWaveResult> {
    let opts = options;
    // as admixtools2 ("Duplicated pops"); qpAdm's left is target + sources
    check_populations(left.iter().chain(right).copied(), |p| {
        f2_blocks.pop_index(p).is_some()
    })?;

    qpwave_on(F4Matrix::from_f2_blocks(f2_blocks, left, right)?, opts)
}

/// qpWave on an f4 matrix already computed (from f2 blocks, or from
/// genotypes as admixtools2's `qpwave(prefix)`).
pub fn qpwave_on(f4_matrix: F4Matrix, opts: &QpWaveOptions) -> Result<QpWaveResult> {
    // admixtools2 drop_ranks tests ranks nrow(f4) - 1 down to 0
    let max_rank = opts.max_rank.unwrap_or(f4_matrix.n_rows() - 1);

    // Perform rank tests
    let rank_tests = compute_rank_tests(&f4_matrix, max_rank, opts.fudge)?;

    let estimated_rank = estimate_rank(&rank_tests, max_rank, opts.alpha);

    Ok(QpWaveResult {
        rank_tests,
        estimated_rank,
        f4_matrix,
    })
}

/// The lowest rank whose test is not rejected (p > alpha); with every rank
/// rejected, the rank is above `max_rank`.
fn estimate_rank(rank_tests: &[RankTest], max_rank: usize, alpha: f64) -> usize {
    rank_tests
        .iter()
        .find(|r| r.p_value > alpha)
        .map_or(max_rank + 1, |r| r.rank)
}

/// Rank tests 0..=max_rank, as admixtools2 `drop_ranks` / `qpadm_fit`.
fn compute_rank_tests(f4_matrix: &F4Matrix, max_rank: usize, fudge: f64) -> Result<Vec<RankTest>> {
    if f4_matrix.values.is_empty() {
        return Ok(vec![]);
    }
    let xmat = fit::to_dmatrix(&f4_matrix.values);
    let qinv = fit::fudged_inverse(&f4_matrix.covariance, fudge)?;
    (0..=max_rank)
        .map(|rank| {
            let rank_fit = fit::fit_rank(&xmat, &qinv, rank, fudge)?;
            let (chi_sq, df) = fit::evaluate_fit(&xmat, &qinv, &rank_fit, rank);
            Ok(RankTest {
                rank,
                chi_sq,
                df,
                p_value: fit::chisq_p_value(chi_sq, df),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use crate::genotypes::{Chromosome, PopulationId};
    use ndarray::Array2;

    fn create_test_f2_blocks() -> F2Blocks {
        // Create synthetic data with known structure
        let pops = ["Pop1", "Pop2", "Pop3", "Pop4", "Pop5"];
        let n_snps = 1000;

        let mut rng_state: u64 = 54321;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 5));
        let counts = Array2::from_elem((n_snps, 5), 50u32);

        for i in 0..n_snps {
            // Create populations with structured relationships
            // Pop1, Pop2, Pop3 form a clade
            // Pop4, Pop5 are outgroups
            let base1 = rand() * 0.8 + 0.1;
            let base2 = rand() * 0.8 + 0.1;

            freqs[[i, 0]] = (base1 + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 1]] = (base1 + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 2]] = (base1 + rand() * 0.05).clamp(0.0, 1.0);
            freqs[[i, 3]] = (base2 + rand() * 0.1).clamp(0.0, 1.0);
            freqs[[i, 4]] = rand();
        }

        let mut chromosomes = Vec::with_capacity(n_snps);
        let mut positions = Vec::with_capacity(n_snps);

        for i in 0..n_snps {
            chromosomes.push(Chromosome::Autosome(1));
            positions.push(i as f64 * 0.0005);
        }

        let af = AlleleFrequencies::from_positions(
            freqs,
            counts,
            pops.iter().map(|s| PopulationId::from(*s)).collect(),
            chromosomes,
            positions,
        );

        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_qpwave_basic() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2", "Pop3"];
        let right = vec!["Pop4", "Pop5"];

        let result = qpwave(&f2_blocks, &left, &right, &Default::default()).unwrap();

        // Should have rank tests from 0 to min(n_left-1, n_right-1) = min(2, 1) = 1
        assert!(!result.rank_tests.is_empty());

        // Each rank test should have valid values; one with 0 degrees of
        // freedom tests nothing and has no p-value
        for test in &result.rank_tests {
            assert!(test.chi_sq >= 0.0);
            if test.df == 0 {
                assert!(test.p_value.is_nan());
            } else {
                assert!(test.p_value >= 0.0 && test.p_value <= 1.0);
            }
        }
    }

    fn rank_test(rank: usize, p_value: f64) -> RankTest {
        RankTest {
            rank,
            chi_sq: 1.0,
            df: 1,
            p_value,
        }
    }

    #[test]
    fn test_estimated_rank_is_lowest_rank_not_rejected() {
        let tests = [rank_test(0, 0.01), rank_test(1, 0.3), rank_test(2, 0.9)];
        assert_eq!(estimate_rank(&tests, 2, 0.05), 1);
    }

    #[test]
    fn alpha_is_the_level_of_the_rank_tests() {
        let tests = [rank_test(0, 0.01), rank_test(1, 0.067), rank_test(2, 0.9)];
        assert_eq!(estimate_rank(&tests, 2, 0.05), 1);
        assert_eq!(estimate_rank(&tests, 2, 0.1), 2);
        assert_eq!(estimate_rank(&tests, 2, 0.001), 0);
    }

    #[test]
    fn test_estimated_rank_with_every_rank_rejected_is_above_max_rank() {
        // --left French_1,French_2: rank 0, the only rank tested, has p = 0.039
        assert_eq!(estimate_rank(&[rank_test(0, 0.039)], 0, 0.05), 1);
        let tests = [rank_test(0, 1e-9), rank_test(1, 0.01)];
        assert_eq!(estimate_rank(&tests, 1, 0.05), 2);
    }

    #[test]
    fn test_qpwave_requires_min_populations() {
        let f2_blocks = create_test_f2_blocks();

        // Should fail with only 1 left population
        let result = qpwave(
            &f2_blocks,
            &["Pop1"],
            &["Pop4", "Pop5"],
            &Default::default(),
        );
        assert!(result.is_err());

        // Should fail with only 1 right population
        let result = qpwave(
            &f2_blocks,
            &["Pop1", "Pop2"],
            &["Pop4"],
            &Default::default(),
        );
        assert!(result.is_err());
    }
}
