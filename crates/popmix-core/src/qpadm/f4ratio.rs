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
//! The admixture proportion α (fraction of X from the A branch) is, as in
//! admixtools2 `qpf4ratio(pops = c(A, O, X, B, C))`:
//!
//! ```text
//! α = f4(A, O; X, B) / f4(A, O; C, B)
//! ```
//!
//! Where X is admixed between A and B branches, and C is a reference population
//! on the B branch. (By f4 additivity this equals 1 - f4(A, O; X, C) / f4(A, O; B, C).)
//!
//! The ratio is formed per block and jackknifed as admixtools2 does: blocks
//! whose denominator is below the threshold are left out, the estimate is the
//! ratio of the length-weighted mean numerator and denominator, and the
//! leave-one-out ratios go through `jack_mat_stats` with that total.

use crate::fstats::{jackknife_with_total, leave_one_out, F2Blocks};
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
    /// Blocks whose denominator f4 is smaller than this in absolute value are
    /// left out (admixtools2 `thresh`, default: 1e-6).
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
/// Formula: α = f4(A, O; X, B) / f4(A, O; C, B), as admixtools2 `qpf4ratio`
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

    let f2 = |a: &str, b: &str| {
        f2_blocks
            .get_block_f2_values(a, b)
            .ok_or_else(|| PopmixError::Numerical(format!("Failed to get f2({a}, {b})")))
    };
    let (a, o, x, b, c) = (source_a, outgroup, test, source_b, reference);
    // f2_f4(f2_14, f2_23, f2_13, f2_24) for f4(A, O; X, B) and f4(A, O; C, B)
    let (f2_ab, f2_ox, f2_ax, f2_ob) = (f2(a, b)?, f2(o, x)?, f2(a, x)?, f2(o, b)?);
    let (f2_oc, f2_ac) = (f2(o, c)?, f2(a, c)?);

    let n_blocks = f2_blocks.n_blocks;
    let mut num = Vec::with_capacity(n_blocks);
    let mut den = Vec::with_capacity(n_blocks);
    for k in 0..n_blocks {
        let nk = (f2_ab[k] + f2_ox[k] - f2_ax[k] - f2_ob[k]) / 2.0;
        let dk = (f2_ab[k] + f2_oc[k] - f2_ac[k] - f2_ob[k]) / 2.0;
        // setmiss = abs(f4_den) < thresh
        if dk.is_nan() || dk.abs() < opts.min_denominator || !nk.is_finite() {
            num.push(f64::NAN);
            den.push(f64::NAN);
        } else {
            num.push(nk);
            den.push(dk);
        }
    }

    let valid_blocks = den.iter().filter(|v| v.is_finite()).count();
    if valid_blocks < 2 {
        return Err(PopmixError::Numerical(
            "Insufficient valid blocks for f4ratio computation".to_string(),
        ));
    }

    let lengths = f2_blocks.block_lengths_f64();
    let weighted_mean = |v: &[f64]| {
        let (s, w) = v
            .iter()
            .zip(&lengths)
            .filter(|(x, _)| x.is_finite())
            .fold((0.0, 0.0), |(s, w), (x, l)| (s + x * l, w + l));
        s / w
    };
    let f4_numerator = weighted_mean(&num);
    let f4_denominator = weighted_mean(&den);
    let tot = f4_numerator / f4_denominator;

    let loo_ratio: Vec<f64> = leave_one_out(&num, &lengths)
        .into_iter()
        .zip(leave_one_out(&den, &lengths))
        .map(|(n, d)| n / d)
        .collect();
    let (alpha, std_err) = jackknife_with_total(&loo_ratio, &lengths, tot).ok_or_else(|| {
        PopmixError::Numerical("Insufficient valid blocks for f4ratio computation".to_string())
    })?;

    let z_score = if std_err > 0.0 { alpha / std_err } else { 0.0 };

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
            &f2_blocks, "PopX", // test (admixed)
            "PopA", // source_a
            "PopB", // source_b
            "PopO", // outgroup
            "PopC", // reference
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
        assert!(
            result.std_err >= 0.0,
            "Standard error should be non-negative"
        );

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
            &f2_blocks,
            "PopX",
            "PopA",
            "PopB",
            "PopO",
            "PopC",
            Some(opts),
        );

        assert!(result.is_ok());
    }
}
