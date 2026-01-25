//! qpF4ratio: Admixture proportion estimation via f4 ratios.
//!
//! This module implements the f4-ratio method for estimating admixture proportions.
//! Given 5 populations arranged as:
//!
//! ```text
//!          O (outgroup)
//!          |
//!     +---------+
//!     |         |
//!     A         |
//!     |    +----+----+
//!     |    |         |
//!     X   B         C (reference)
//! ```
//!
//! The admixture proportion α (fraction of X from the A branch) is:
//!
//! ```text
//! α = f4(A, O; X, C) / f4(A, O; B, C)
//! ```
//!
//! Where X is admixed between A and B branches, and C is a reference population
//! on the B branch.

use crate::fstats::F2Blocks;
use crate::genotypes::PopulationId;
use crate::{PopmixError, Result};

/// Result of qpF4ratio analysis.
#[derive(Debug, Clone)]
pub struct F4RatioResult {
    /// Test (admixed) population.
    pub test: PopulationId,
    /// Source population A (α proportion comes from this branch).
    pub source_a: PopulationId,
    /// Source population B ((1-α) proportion comes from this branch).
    pub source_b: PopulationId,
    /// Outgroup population.
    pub outgroup: PopulationId,
    /// Reference population (on B branch).
    pub reference: PopulationId,
    /// Estimated admixture proportion (α).
    pub alpha: f64,
    /// Standard error of α (from block jackknife).
    pub std_err: f64,
    /// Z-score (alpha / std_err).
    pub z_score: f64,
    /// Number of blocks used.
    pub n_blocks: usize,
    /// Numerator f4 value: f4(A, O; X, C).
    pub f4_numerator: f64,
    /// Denominator f4 value: f4(A, O; B, C).
    pub f4_denominator: f64,
}

/// Options for f4ratio computation.
#[derive(Debug, Clone)]
pub struct F4RatioOptions {
    /// Minimum absolute value for denominator f4 (default: 1e-6).
    /// If |f4_denominator| < min_denominator, returns an error.
    pub min_denominator: f64,
}

impl Default for F4RatioOptions {
    fn default() -> Self {
        Self {
            min_denominator: 1e-6,
        }
    }
}

/// Compute f4 ratio for admixture proportion estimation.
///
/// Formula: α = f4(A, O; X, C) / f4(A, O; B, C)
///
/// # Arguments
/// * `f2_blocks` - Pre-computed F2 blocks
/// * `test` - Test (admixed) population X
/// * `source_a` - Source population A (α proportion)
/// * `source_b` - Source population B ((1-α) proportion)
/// * `outgroup` - Outgroup population O
/// * `reference` - Reference population C (on B branch)
/// * `options` - Optional computation settings
///
/// # Returns
/// F4RatioResult with estimated α and standard error
///
/// # Errors
/// * If any population is not found in f2_blocks
/// * If denominator f4 is too close to zero
/// * If insufficient blocks with valid data
pub fn f4ratio(
    f2_blocks: &F2Blocks,
    test: &str,
    source_a: &str,
    source_b: &str,
    outgroup: &str,
    reference: &str,
    options: Option<F4RatioOptions>,
) -> Result<F4RatioResult> {
    let opts = options.unwrap_or_default();

    // Validate populations exist
    for (pop, name) in [
        (test, "test"),
        (source_a, "source_a"),
        (source_b, "source_b"),
        (outgroup, "outgroup"),
        (reference, "reference"),
    ] {
        if f2_blocks.pop_index(pop).is_none() {
            return Err(PopmixError::PopulationNotFound(format!(
                "Population '{}' ({}) not found in F2 blocks",
                pop, name
            )));
        }
    }

    // Compute block-level f4 values for numerator and denominator
    // f4(A, O; X, C) = (f2(A,C) + f2(O,X) - f2(A,X) - f2(O,C)) / 2
    // f4(A, O; B, C) = (f2(A,C) + f2(O,B) - f2(A,B) - f2(O,C)) / 2

    let f2_ac = f2_blocks
        .get_block_f2_values(source_a, reference)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(A,C)".to_string()))?;
    let f2_ox = f2_blocks
        .get_block_f2_values(outgroup, test)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(O,X)".to_string()))?;
    let f2_ax = f2_blocks
        .get_block_f2_values(source_a, test)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(A,X)".to_string()))?;
    let f2_oc = f2_blocks
        .get_block_f2_values(outgroup, reference)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(O,C)".to_string()))?;
    let f2_ob = f2_blocks
        .get_block_f2_values(outgroup, source_b)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(O,B)".to_string()))?;
    let f2_ab = f2_blocks
        .get_block_f2_values(source_a, source_b)
        .ok_or_else(|| PopmixError::Numerical("Failed to get f2(A,B)".to_string()))?;

    let n_blocks = f2_blocks.n_blocks;

    // Compute block-level f4 values
    let mut block_f4_num = Vec::with_capacity(n_blocks);
    let mut block_f4_den = Vec::with_capacity(n_blocks);

    for i in 0..n_blocks {
        let f4_num = (f2_ac[i] + f2_ox[i] - f2_ax[i] - f2_oc[i]) / 2.0;
        let f4_den = (f2_ac[i] + f2_ob[i] - f2_ab[i] - f2_oc[i]) / 2.0;

        // Skip blocks with NaN or zero denominator
        if f4_num.is_nan() || f4_den.is_nan() || f4_den.abs() < 1e-15 {
            continue;
        }

        block_f4_num.push(f4_num);
        block_f4_den.push(f4_den);
    }

    let valid_blocks = block_f4_num.len();
    if valid_blocks < 2 {
        return Err(PopmixError::Numerical(
            "Insufficient valid blocks for f4ratio computation".to_string(),
        ));
    }

    // Compute overall f4 values (mean of block estimates)
    let f4_numerator: f64 = block_f4_num.iter().sum::<f64>() / valid_blocks as f64;
    let f4_denominator: f64 = block_f4_den.iter().sum::<f64>() / valid_blocks as f64;

    // Check denominator magnitude
    if f4_denominator.abs() < opts.min_denominator {
        return Err(PopmixError::Numerical(format!(
            "f4 denominator too small: {} (threshold: {})",
            f4_denominator, opts.min_denominator
        )));
    }

    // Compute overall ratio and alpha
    let ratio = f4_numerator / f4_denominator;
    // α = 1 - ratio (derived from: f4(A,O;X,C) = (1-α)*f4(A,O;B,C))
    // α represents the proportion of X from source_a
    let alpha = 1.0 - ratio;

    // Compute jackknife standard error for the ratio
    // SE(alpha) = SE(1-ratio) = SE(ratio)
    let std_err = jackknife_std_err_ratio(&block_f4_num, &block_f4_den, ratio);

    let z_score = if std_err > 0.0 {
        alpha / std_err
    } else {
        0.0
    };

    Ok(F4RatioResult {
        test: PopulationId::new(test),
        source_a: PopulationId::new(source_a),
        source_b: PopulationId::new(source_b),
        outgroup: PopulationId::new(outgroup),
        reference: PopulationId::new(reference),
        alpha,
        std_err,
        z_score,
        n_blocks: valid_blocks,
        f4_numerator,
        f4_denominator,
    })
}

/// Compute jackknife standard error for a ratio of sums.
///
/// For ratio = sum(numerators) / sum(denominators), compute SE using
/// delete-one jackknife where each LOO estimate is:
///   (total_num - num_i) / (total_den - den_i)
fn jackknife_std_err_ratio(
    block_numerators: &[f64],
    block_denominators: &[f64],
    full_estimate: f64,
) -> f64 {
    let n = block_numerators.len();
    if n < 2 {
        return 0.0;
    }

    let total_num: f64 = block_numerators.iter().sum();
    let total_den: f64 = block_denominators.iter().sum();

    let mut sum_sq_diff = 0.0;
    for i in 0..n {
        // Leave-one-out estimate using ratio of sums
        let loo_num = total_num - block_numerators[i];
        let loo_den = total_den - block_denominators[i];
        let loo_estimate = loo_num / loo_den;
        let diff = loo_estimate - full_estimate;
        sum_sq_diff += diff * diff;
    }

    let variance = ((n - 1) as f64 / n as f64) * sum_sq_diff;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::blocks::BlockDefinition;
    use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    /// Create test data with known admixture proportions.
    fn create_test_f2_blocks() -> F2Blocks {
        // Create 5 populations: A, B, X (admixed), O (outgroup), C (reference)
        // X = 0.6 * A + 0.4 * B (true alpha = 0.6)

        let mut samples = Vec::new();
        for i in 0..20 {
            samples.push(Sample::new(format!("A{}", i), PopulationId::new("PopA")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("B{}", i), PopulationId::new("PopB")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("X{}", i), PopulationId::new("PopX")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("O{}", i), PopulationId::new("PopO")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("C{}", i), PopulationId::new("PopC")));
        }

        // Create SNPs across multiple blocks
        // 200 SNPs, 0.1 cM spacing = 20 cM total, with 5 cM block size = 4 blocks
        let snps: Vec<Snp> = (0..200)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.1, // 0.1 cM spacing
                    'A',
                    'G',
                )
            })
            .collect();

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        // Use a simple random-ish pattern based on indices
        for snp_idx in 0..n_snps {
            // Base frequencies vary by SNP
            let base = (snp_idx * 7 % 100) as f64 / 100.0;

            // PopA: low frequency (0.1 - 0.4)
            let freq_a = 0.1 + base * 0.3;
            // PopB: high frequency (0.6 - 0.9)
            let freq_b = 0.6 + base * 0.3;
            // PopX: admixed (0.6*A + 0.4*B)
            let freq_x = 0.6 * freq_a + 0.4 * freq_b;
            // PopO: distant outgroup (0.05 - 0.15)
            let freq_o = 0.05 + base * 0.1;
            // PopC: similar to B (0.55 - 0.85)
            let freq_c = 0.55 + base * 0.3;

            for sample_idx in 0..n_samples {
                let pop_idx = sample_idx / 20;
                let freq = match pop_idx {
                    0 => freq_a,
                    1 => freq_b,
                    2 => freq_x,
                    3 => freq_o,
                    4 => freq_c,
                    _ => 0.5,
                };

                // Assign genotype based on frequency and sample position
                // This creates variation within populations
                let threshold = freq * 2.0; // Scale to [0, 2] for genotype assignment
                let noise = ((sample_idx * 13 + snp_idx * 17) % 100) as f64 / 100.0 * 2.0;
                data[[snp_idx, sample_idx]] = if noise < threshold { 2 } else { 0 };
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);

        let block_def = BlockDefinition::new(0.05); // 5 cM blocks
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_f4ratio_basic() {
        let f2_blocks = create_test_f2_blocks();

        let result = f4ratio(
            &f2_blocks,
            "PopX",  // test (admixed)
            "PopA",  // source_a
            "PopB",  // source_b
            "PopO",  // outgroup
            "PopC",  // reference
            None,
        )
        .unwrap();

        // Verify the result structure is valid
        assert_eq!(result.test.as_str(), "PopX");
        assert_eq!(result.source_a.as_str(), "PopA");
        assert_eq!(result.source_b.as_str(), "PopB");
        assert_eq!(result.outgroup.as_str(), "PopO");
        assert_eq!(result.reference.as_str(), "PopC");

        // Alpha should be finite
        assert!(result.alpha.is_finite(), "Alpha should be finite");

        // Standard error should be positive (unless all blocks identical)
        assert!(result.std_err >= 0.0, "Standard error should be non-negative");

        // Should have at least 2 valid blocks
        assert!(result.n_blocks >= 2, "Should have at least 2 valid blocks");

        // f4 values should be non-zero (populations are distinct)
        assert!(
            result.f4_denominator.abs() > 1e-10,
            "f4_denominator should be non-zero"
        );
    }

    #[test]
    fn test_f4ratio_population_not_found() {
        let f2_blocks = create_test_f2_blocks();

        let result = f4ratio(
            &f2_blocks,
            "NonExistent",
            "PopA",
            "PopB",
            "PopO",
            "PopC",
            None,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_f4ratio_options() {
        let f2_blocks = create_test_f2_blocks();

        let opts = F4RatioOptions {
            min_denominator: 1e-10,
        };

        let result = f4ratio(
            &f2_blocks, "PopX", "PopA", "PopB", "PopO", "PopC",
            Some(opts),
        );

        assert!(result.is_ok());
    }
}
