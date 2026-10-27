//! qpAdm admixture weight estimation.
//!
//! qpAdm estimates the proportions of ancestry in a target population
//! from a set of source (left) populations, using outgroup (right)
//! populations to constrain the model.
//!
//! Follows admixtools2 `qpadm`: the f4 matrix f4(source, target; R, R0) is fit
//! at rank (sources - 1) (see [`super::fit`]); the weights come from that fit,
//! their standard errors from refitting on every leave-one-block-out f4
//! matrix, and the model's chi-squared is the qpWave rank test at that rank.
//! Weights can be constrained to [0,1] or left unconstrained.

use super::compat::Admixtools2Compat;
use super::fit;
use super::qpwave::{qpwave, QpWaveOptions, QpWaveResult};
use crate::fstats::F2Blocks;
use crate::genotypes::PopulationId;
use crate::Result;

/// Configuration options for qpAdm analysis.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpAdmOptions {
    /// Regularization factor (default: 1e-4).
    pub fudge: f64,
    /// Whether to constrain weights to [0, 1] (default: false).
    pub constrained: bool,
    /// admixtools2 quirks to reproduce (default: none).
    pub compat: Admixtools2Compat,
}

impl Default for QpAdmOptions {
    fn default() -> Self {
        Self {
            fudge: 1e-4,
            constrained: false,
            compat: Admixtools2Compat::CORRECTED,
        }
    }
}

/// Result of qpAdm analysis.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpAdmResult {
    /// Target population.
    pub target: PopulationId,
    /// Source populations.
    pub sources: Vec<PopulationId>,
    /// Estimated admixture weights (sum to 1).
    pub weights: Vec<f64>,
    /// Jackknife standard errors of the weights.
    pub weight_std_errors: Vec<f64>,
    /// Chi-squared of the model (qpWave rank test at rank sources - 1).
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
    /// P-value of the model fit.
    pub p_value: f64,
    /// Whether the weights were constrained to [0, 1].
    pub constrained: bool,
    /// qpWave result for the left populations (target first, then sources).
    pub qpwave_result: QpWaveResult,
    /// Fit of every subset of the sources (admixtools2 `popdrop`).
    pub popdrop: Vec<fit::PopDropRow>,
}

impl QpAdmResult {
    /// Whether the model fits: p-value above the threshold.
    pub fn is_good_fit(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }

    /// Whether all weights lie in [0, 1].
    pub fn is_feasible(&self) -> bool {
        self.weights.iter().all(|&w| (0.0..=1.0).contains(&w))
    }
}

/// Run qpAdm: estimate the weights of `sources` in `target`.
///
/// # Arguments
/// * `f2_blocks` - Precomputed F2 statistics
/// * `target` - Target population
/// * `sources` - Source populations
/// * `right` - Right (outgroup) populations; the first is the reference
/// * `options` - Optional analysis parameters
pub fn qpadm(
    f2_blocks: &F2Blocks,
    target: &str,
    sources: &[&str],
    right: &[&str],
    options: &QpAdmOptions,
) -> Result<QpAdmResult> {
    let opts = options;

    if sources.is_empty() {
        return Err(crate::PopmixError::too_few(
            "qpAdm",
            "source populations",
            1,
            0,
        ));
    }
    if right.len() < 2 {
        return Err(crate::PopmixError::too_few(
            "qpAdm",
            "right populations",
            2,
            right.len(),
        ));
    }

    // admixtools2: left = c(target, sources)
    let mut left: Vec<&str> = Vec::with_capacity(sources.len() + 1);
    left.push(target);
    left.extend_from_slice(sources);

    let qpwave_opts = QpWaveOptions {
        fudge: opts.fudge,
        max_rank: None,
    };
    let qpwave_result = qpwave(f2_blocks, &left, right, &qpwave_opts)?;
    let f4 = &qpwave_result.f4_matrix;

    let rnk = sources.len() - 1;
    let xmat = fit::to_dmatrix(&f4.values);
    let qinv = fit::fudged_inverse(&f4.covariance, opts.fudge)?;
    let weights = fit::qpadm_weights(&xmat, &qinv, rnk, opts.fudge, opts.constrained)?;
    let weight_std_errors =
        fit::qpadm_weight_std_errors(&f4.leave_one_out, &qinv, rnk, opts.fudge, opts.constrained)?;
    let popdrop = fit::drop_pops(&xmat, &qinv, opts.fudge, opts.constrained, opts.compat)?;

    let (chi_sq, df, p_value) = qpwave_result
        .rank_tests
        .iter()
        .find(|r| r.rank == rnk)
        .map(|r| (r.chi_sq, r.df, r.p_value))
        .ok_or_else(|| {
            crate::PopmixError::Internal(format!("no qpWave rank test at rank {rnk}"))
        })?;

    Ok(QpAdmResult {
        target: PopulationId::from(target),
        sources: sources.iter().map(|s| PopulationId::from(*s)).collect(),
        weights,
        weight_std_errors,
        chi_sq,
        df,
        p_value,
        constrained: opts.constrained,
        qpwave_result,
        popdrop,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use crate::genotypes::{Chromosome, PopulationId};
    use ndarray::Array2;

    fn create_admixed_test_data() -> F2Blocks {
        // Create data where Pop3 is a 50/50 mix of Pop1 and Pop2
        let pops = ["Pop1", "Pop2", "Pop3", "Out1", "Out2"];
        let n_snps = 2000;

        let mut rng_state: u64 = 98765;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 5));
        let counts = Array2::from_elem((n_snps, 5), 100u32);

        for i in 0..n_snps {
            // Pop1 and Pop2 are divergent sources
            let freq1 = rand() * 0.8 + 0.1;
            let freq2 = rand() * 0.8 + 0.1;

            freqs[[i, 0]] = freq1; // Pop1
            freqs[[i, 1]] = freq2; // Pop2

            // Pop3 is 50/50 admixture + small noise
            freqs[[i, 2]] = (0.5 * freq1 + 0.5 * freq2 + rand() * 0.02 - 0.01).clamp(0.0, 1.0);

            // Outgroups
            freqs[[i, 3]] = rand();
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
    fn test_qpadm_basic() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1", "Out2"],
            &Default::default(),
        )
        .unwrap();

        println!("qpAdm result:");
        println!("  Weights: {:?}", result.weights);
        println!(
            "  Chi-sq: {}, df: {}, p: {}",
            result.chi_sq, result.df, result.p_value
        );

        // Weights should sum to 1
        let sum: f64 = result.weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 0.01,
            "Weights should sum to 1, got {sum}"
        );

        // Both weights should be positive (since Pop3 is mixture of both)
        // Note: unconstrained can give negative weights for noisy data
    }

    #[test]
    fn test_qpadm_constrained() {
        let f2_blocks = create_admixed_test_data();

        let opts = QpAdmOptions {
            constrained: true,
            ..Default::default()
        };

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1", "Out2"],
            &opts,
        )
        .unwrap();

        println!("Constrained qpAdm:");
        println!("  Weights: {:?}", result.weights);

        // All weights should be in [0, 1]
        for w in &result.weights {
            assert!(*w >= 0.0 && *w <= 1.0, "Weight {w} out of bounds");
        }

        // Weights should sum to 1
        let sum: f64 = result.weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 0.01,
            "Weights should sum to 1, got {sum}"
        );
    }

    #[test]
    fn test_qpadm_requires_sources() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &[],
            &["Out1", "Out2"],
            &Default::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_qpadm_requires_outgroups() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1"],
            &Default::default(),
        );
        assert!(result.is_err());
    }
}
