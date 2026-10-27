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
//! The ratio is jackknifed over blocks: blocks whose denominator is below the
//! threshold in absolute value (or with a non-finite value) are left out, the
//! estimate is the ratio of the length-weighted mean numerator and
//! denominator, and the leave-one-out ratios (leave-one-out numerator over
//! leave-one-out denominator) go through the weighted block jackknife with
//! that total ([`jackknife_with_total`]).
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/small-helpers.md` (section 1), after Reich et al. 2009
//! (Nature 461:489, supplement) and Patterson et al. 2012 (Genetics
//! 192:1065), and checked against admixtools2's black-box outputs.

use crate::fstats::jackknife::z_score;
use crate::fstats::{jackknife_with_total, leave_one_out, F2Blocks};
use crate::genotypes::PopulationId;
use crate::{PopmixError, Result};

/// Result of qpF4ratio analysis.
#[derive(Debug, Clone)]
#[non_exhaustive]
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
    /// Numerator f4 value: f4(A, O; X, B), the mean over the blocks used.
    pub f4_numerator: f64,
    /// Denominator f4 value: f4(A, O; C, B), the mean over the blocks used.
    pub f4_denominator: f64,
}

/// Options for f4ratio computation.
#[derive(Debug, Clone)]
#[non_exhaustive]
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
    options: &F4RatioOptions,
) -> Result<F4RatioResult> {
    let opts = options;
    for pop in [test, source_a, source_b, outgroup, reference] {
        if f2_blocks.pop_index(pop).is_none() {
            return Err(PopmixError::PopulationNotFound(pop.to_string()));
        }
    }
    let f2 = |p: &str, q: &str| {
        f2_blocks
            .get_block_f2_values(p, q)
            .ok_or_else(|| PopmixError::PopulationNotFound(format!("{p} or {q}")))
    };
    // f4(A, B; C, D) = (f2(A, D) + f2(B, C) - f2(A, C) - f2(B, D)) / 2
    let f4 = |a: &str, b: &str, c: &str, d: &str| -> Result<Vec<f64>> {
        Ok(((f2(a, d)? + f2(b, c)? - f2(a, c)? - f2(b, d)?) / 2.0).to_vec())
    };
    let num_all = f4(source_a, outgroup, test, source_b)?;
    let den_all = f4(source_a, outgroup, reference, source_b)?;

    // Blocks with a usable denominator and finite values; the rest are absent.
    let (mut num, mut den, mut lengths) = (Vec::new(), Vec::new(), Vec::new());
    for ((n, d), m) in num_all
        .iter()
        .zip(&den_all)
        .zip(f2_blocks.block_lengths_f64())
    {
        if n.is_finite() && d.is_finite() && d.abs() >= opts.min_denominator {
            num.push(*n);
            den.push(*d);
            lengths.push(m);
        }
    }
    if num.len() < 2 {
        return Err(PopmixError::too_few(
            "the f4ratio jackknife",
            format!(
                "blocks with a denominator of at least {} in absolute value",
                opts.min_denominator
            ),
            2,
            num.len(),
        ));
    }

    let total: f64 = lengths.iter().sum();
    let mean = |x: &[f64]| x.iter().zip(&lengths).map(|(v, m)| v * m).sum::<f64>() / total;
    let f4_numerator = mean(&num);
    let f4_denominator = mean(&den);
    let tot = f4_numerator / f4_denominator;

    let loo: Vec<f64> = leave_one_out(&num, &lengths)
        .iter()
        .zip(leave_one_out(&den, &lengths))
        .map(|(n, d)| n / d)
        .collect();
    let (alpha, std_err) = jackknife_with_total(&loo, &lengths, tot)
        .ok_or_else(|| PopmixError::data("f4ratio: no block has a finite leave-one-out ratio"))?;

    Ok(F4RatioResult {
        test: PopulationId::new(test),
        source_a: PopulationId::new(source_a),
        source_b: PopulationId::new(source_b),
        outgroup: PopulationId::new(outgroup),
        reference: PopulationId::new(reference),
        alpha,
        std_err,
        z_score: z_score(alpha, std_err),
        n_blocks: num.len(),
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
            samples.push(Sample::new(format!("A{i}"), PopulationId::new("PopA")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("B{i}"), PopulationId::new("PopB")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("X{i}"), PopulationId::new("PopX")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("O{i}"), PopulationId::new("PopO")));
        }
        for i in 0..20 {
            samples.push(Sample::new(format!("C{i}"), PopulationId::new("PopC")));
        }

        // Create SNPs across multiple blocks
        // 200 SNPs, 0.1 cM spacing = 20 cM total, with 5 cM block size = 4 blocks
        let snps: Vec<Snp> = (0..200)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    f64::from(i) * 0.1, // 0.1 cM spacing
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
            "PopX", // test (admixed)
            "PopA", // source_a
            "PopB", // source_b
            "PopO", // outgroup
            "PopC", // reference
            &Default::default(),
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
            &Default::default(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_f4ratio_options() {
        let f2_blocks = create_test_f2_blocks();

        let opts = F4RatioOptions {
            min_denominator: 1e-10,
        };

        let result = f4ratio(&f2_blocks, "PopX", "PopA", "PopB", "PopO", "PopC", &opts);

        assert!(result.is_ok());
    }
}
