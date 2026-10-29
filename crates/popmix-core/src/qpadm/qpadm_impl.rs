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
use super::f4_matrix::F4Matrix;
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
        ..QpWaveOptions::default()
    };
    let qpwave_result = qpwave(f2_blocks, &left, right, &qpwave_opts)?;
    qpadm_from_qpwave(target, sources, qpwave_result, opts)
}

/// qpAdm on an f4 matrix already computed, its left populations the target
/// then `sources` (from genotypes, as admixtools2's `qpadm(prefix)`).
pub fn qpadm_on(
    target: &str,
    sources: &[&str],
    f4_matrix: F4Matrix,
    options: &QpAdmOptions,
) -> Result<QpAdmResult> {
    if sources.is_empty() || f4_matrix.n_rows() != sources.len() || f4_matrix.n_cols() < 1 {
        return Err(crate::PopmixError::DimensionMismatch {
            context: "qpAdm f4 matrix".to_string(),
            expected: format!("{} rows, one per source, and a column", sources.len()),
            actual: format!("{} x {}", f4_matrix.n_rows(), f4_matrix.n_cols()),
        });
    }
    let qpwave_opts = QpWaveOptions {
        fudge: options.fudge,
        ..QpWaveOptions::default()
    };
    let qpwave_result = super::qpwave::qpwave_on(f4_matrix, &qpwave_opts)?;
    qpadm_from_qpwave(target, sources, qpwave_result, options)
}

/// Weights, SEs, fit and popdrop from the qpWave result of the target and
/// sources.
fn qpadm_from_qpwave(
    target: &str,
    sources: &[&str],
    qpwave_result: QpWaveResult,
    opts: &QpAdmOptions,
) -> Result<QpAdmResult> {
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

            // Outgroups, each drifted from one source, so that they tell the
            // sources apart and the weights are identifiable
            freqs[[i, 3]] = (freq1 + rand() * 0.2 - 0.1).clamp(0.0, 1.0);
            freqs[[i, 4]] = (freq2 + rand() * 0.2 - 0.1).clamp(0.0, 1.0);
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
    fn test_qpadm_on_an_f4_matrix_is_qpadm_with_its_fudge() {
        let f2_blocks = create_admixed_test_data();
        let opts = QpAdmOptions {
            fudge: 0.1,
            ..Default::default()
        };
        let right = ["Out1", "Out2", "Pop2"];
        let f4 = F4Matrix::from_f2_blocks(&f2_blocks, &["Pop3", "Pop1"], &right).unwrap();
        let on = qpadm_on("Pop3", &["Pop1"], f4, &opts).unwrap();
        let full = qpadm(&f2_blocks, "Pop3", &["Pop1"], &right, &opts).unwrap();
        assert_eq!(on.chi_sq, full.chi_sq);
    }

    #[test]
    fn test_qpadm_on_needs_a_row_per_source() {
        let f2_blocks = create_admixed_test_data();
        let left = ["Pop3", "Pop1", "Pop2"];
        let f4 = F4Matrix::from_f2_blocks(&f2_blocks, &left, &["Out1", "Out2"]).unwrap();
        let err = qpadm_on("Pop3", &["Pop1"], f4, &Default::default()).unwrap_err();
        assert!(
            matches!(&err, crate::PopmixError::DimensionMismatch { context, .. } if context == "qpAdm f4 matrix"),
            "{err:?}"
        );
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

        // Pop3 is a 50/50 mixture of the sources
        for w in &result.weights {
            assert!((w - 0.5).abs() < 0.05, "weights {:?}", result.weights);
        }
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

        // inside [0, 1] already, so the constraint changes nothing
        let free = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1", "Out2"],
            &Default::default(),
        )
        .unwrap();
        for (c, f) in result.weights.iter().zip(&free.weights) {
            assert!(
                (c - f).abs() < 1e-9,
                "constrained {:?}, free {:?}",
                result.weights,
                free.weights
            );
        }
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
