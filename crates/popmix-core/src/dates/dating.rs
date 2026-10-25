//! Main DATES algorithm implementation.
//!
//! Estimates admixture timing from LD decay patterns.

use super::compat::DatesCompat;
use super::exponential_fit::{fit_exponential, weighted_jackknife, DecayFit, FitOptions};
use super::ld_decay::{compute_ld_decay_by_chromosome, LdBin, LdOptions};
use super::weights::{compute_weights, WeightResult};
use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::PopulationId;
use crate::Result;
use rayon::prelude::*;

/// Options for DATES algorithm.
#[derive(Debug, Clone)]
pub struct DatesOptions {
    /// Options for LD computation.
    pub ld_options: LdOptions,
    /// Options for exponential fitting.
    pub fit_options: FitOptions,
    /// Generation time in years (for converting to years BP). Default: 29.0
    pub generation_time: f64,
    /// Minimum observed allele count in each source population for a SNP to
    /// be used. Default: 1, so one called sample suffices, as in DATES v753
    /// (which reads `mincount` but ignores it: "mincount stripped").
    pub min_source_count: u32,
    /// Run jackknife for standard error estimation. Default: true
    pub run_jackknife: bool,
    /// DATES v753 quirks to reproduce. Applies to every stage, overriding
    /// `ld_options.compat`. Default: none, all corrected; use
    /// [`DatesCompat::V753`] to compare against DATES output.
    pub compat: DatesCompat,
}

impl Default for DatesOptions {
    fn default() -> Self {
        DatesOptions {
            ld_options: LdOptions::default(),
            fit_options: FitOptions::default(),
            generation_time: 29.0,
            min_source_count: 1,
            run_jackknife: true,
            compat: DatesCompat::CORRECTED,
        }
    }
}

/// Round to `decimals` places as printf would print it.
fn round_printed(x: f64, decimals: usize) -> f64 {
    format!("{x:.decimals$}").parse().unwrap_or(x)
}

/// Fit LD bins, first rounding correlations as DATES prints them if asked.
fn fit_bins(bins: &[LdBin], opts: &DatesOptions) -> Result<DecayFit> {
    if opts.compat.fit_printed_values {
        let rounded: Vec<LdBin> = bins
            .iter()
            .map(|b| LdBin {
                correlation: round_printed(b.correlation, 6),
                ..b.clone()
            })
            .collect();
        fit_exponential(&rounded, Some(opts.fit_options.clone()))
    } else {
        fit_exponential(bins, Some(opts.fit_options.clone()))
    }
}

/// Result of DATES analysis.
#[derive(Debug, Clone)]
pub struct DatesResult {
    /// Test (admixed) population.
    pub test_pop: PopulationId,
    /// First source population.
    pub source1: PopulationId,
    /// Second source population.
    pub source2: PopulationId,
    /// Estimated number of generations since admixture: the weighted
    /// jackknife estimate if the jackknife ran (as DATES reports), otherwise
    /// the full-data fit.
    pub generations: f64,
    /// Standard error of generations (from jackknife).
    pub generations_se: f64,
    /// Estimated years before present.
    pub years_bp: f64,
    /// Standard error of years BP.
    pub years_bp_se: f64,
    /// LD decay bins (all chromosomes combined).
    pub ld_bins: Vec<LdBin>,
    /// Exponential decay fit result.
    pub fit: DecayFit,
    /// Number of SNPs used.
    pub n_snps: usize,
    /// Weight statistics.
    pub weights: WeightResult,
    /// Number of chromosomes used.
    pub n_chromosomes: usize,
    /// Per-chromosome jackknife blocks, empty if the jackknife did not run.
    pub jackknife_by_chromosome: Vec<JackknifeBlock>,
}

/// One leave-one-chromosome-out jackknife block (a row of the DATES `.jin`).
#[derive(Debug, Clone)]
pub struct JackknifeBlock {
    /// Chromosome left out (1..=22).
    pub chromosome: u8,
    /// Block weight: number of SNPs on the chromosome in the input.
    pub weight: f64,
    /// Generations estimated without this chromosome.
    pub generations: f64,
}

/// Estimate admixture time using the DATES algorithm.
///
/// DATES estimates when admixture occurred by measuring the decay of
/// linkage disequilibrium with genetic distance. LD from recent
/// admixture decays exponentially, with the decay rate proportional
/// to the number of generations since admixture.
///
/// # Arguments
/// * `geno` - Genotype matrix with all relevant populations
/// * `af` - Allele frequencies (or will be computed if None)
/// * `test_pop` - Name of the test (admixed) population
/// * `source1` - Name of first source population
/// * `source2` - Name of second source population
/// * `options` - Algorithm options
///
/// # Returns
/// DatesResult with estimated admixture timing
///
/// # Example
/// ```ignore
/// use popmix_core::dates::{dates, DatesOptions};
///
/// let result = dates(
///     &geno,
///     Some(&af),
///     "Admixed_Pop",
///     "Source1",
///     "Source2",
///     None,
/// )?;
///
/// println!("Admixture time: {:.0} ± {:.0} generations",
///          result.generations, result.generations_se);
/// println!("Years BP: {:.0} ± {:.0}",
///          result.years_bp, result.years_bp_se);
/// ```
pub fn dates(
    geno: &GenotypeMatrix,
    af: Option<&AlleleFrequencies>,
    test_pop: &str,
    source1: &str,
    source2: &str,
    options: Option<DatesOptions>,
) -> Result<DatesResult> {
    let opts = options.unwrap_or_default();
    crate::options::check_positive("generation_time", opts.generation_time)?;

    // Compute allele frequencies if not provided
    let owned_af;
    let af = match af {
        Some(af) => af,
        None => {
            let pops = vec![test_pop, source1, source2];
            owned_af = AlleleFrequencies::from_genotypes_for_populations(geno, Some(&pops));
            &owned_af
        }
    };

    // Step 1: Compute weights (p_source1 - p_source2)
    let weights = compute_weights(af, source1, source2, Some(opts.min_source_count))?;

    if weights.n_valid == 0 {
        return Err(crate::error::PopmixError::InsufficientData {
            context: "No SNPs with valid weights".to_string(),
            required: 100,
            actual: 0,
        });
    }

    // Get source population frequency arrays
    let source1_idx = af
        .pop_index(source1)
        .ok_or_else(|| crate::error::PopmixError::PopulationNotFound(source1.to_string()))?;
    let source2_idx = af
        .pop_index(source2)
        .ok_or_else(|| crate::error::PopmixError::PopulationNotFound(source2.to_string()))?;

    let source1_freqs: Vec<f64> = (0..geno.n_snps())
        .map(|i| af.frequencies[[i, source1_idx]])
        .collect();
    let source2_freqs: Vec<f64> = (0..geno.n_snps())
        .map(|i| af.frequencies[[i, source2_idx]])
        .collect();

    // Step 2: Compute LD decay, for all chromosomes and leaving each one out
    let ld = compute_ld_decay_by_chromosome(
        geno,
        af,
        test_pop,
        &weights.weights,
        &source1_freqs,
        &source2_freqs,
        Some(LdOptions {
            compat: opts.compat,
            ..opts.ld_options.clone()
        }),
    )?;

    if ld.bins.is_empty() {
        return Err(crate::error::PopmixError::InsufficientData {
            context: "No LD bins computed".to_string(),
            required: 10,
            actual: 0,
        });
    }

    // Step 3: Fit exponential decay
    let mut fit = fit_bins(&ld.bins, &opts)?;

    // Step 4: Weighted jackknife over chromosomes (DATES dates_wtjack + dowtjack)
    let (generations, generations_se, jackknife_by_chromosome) = if opts.run_jackknife {
        let blocks = ld
            .chromosomes
            .par_iter()
            .map(|c| {
                let loo = fit_bins(&c.leave_one_out, &opts)?;
                Ok(JackknifeBlock {
                    chromosome: c.chromosome,
                    weight: c.n_snps as f64,
                    generations: loo.decay_rate,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let printed = |x: f64| {
            if opts.compat.fit_printed_values {
                round_printed(x, 3)
            } else {
                x
            }
        };
        let estimates: Vec<(f64, f64)> = blocks
            .iter()
            .map(|b| (b.weight, printed(b.generations)))
            .collect();
        let (est, se) = weighted_jackknife(printed(fit.decay_rate), &estimates)?;
        fit.decay_rate_se = se;
        (est, se, blocks)
    } else {
        (fit.decay_rate, f64::NAN, Vec::new())
    };

    // Convert to years
    let years_bp = generations * opts.generation_time;
    let years_bp_se = if generations_se.is_finite() {
        generations_se * opts.generation_time
    } else {
        f64::NAN
    };

    // Count chromosomes used
    let n_chromosomes = geno
        .snps()
        .iter()
        .filter(|s| s.is_autosomal())
        .map(|s| &s.chromosome)
        .collect::<std::collections::HashSet<_>>()
        .len();

    Ok(DatesResult {
        test_pop: PopulationId::new(test_pop),
        source1: PopulationId::new(source1),
        source2: PopulationId::new(source2),
        generations,
        generations_se,
        years_bp,
        years_bp_se,
        ld_bins: ld.bins,
        fit,
        n_snps: weights.n_valid,
        weights,
        n_chromosomes,
        jackknife_by_chromosome,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    fn create_test_genotypes() -> GenotypeMatrix {
        // Create test data with 3 populations
        let mut samples = Vec::new();
        for i in 0..20 {
            samples.push(Sample::new(
                format!("test_{}", i),
                PopulationId::new("test"),
            ));
        }
        for i in 0..10 {
            samples.push(Sample::new(
                format!("src1_{}", i),
                PopulationId::new("source1"),
            ));
        }
        for i in 0..10 {
            samples.push(Sample::new(
                format!("src2_{}", i),
                PopulationId::new("source2"),
            ));
        }

        // Create SNPs across multiple chromosomes
        let mut snps = Vec::new();
        for chrom in 1..=3 {
            for i in 0..100 {
                snps.push(Snp::new(
                    format!("rs{}_{}", chrom, i),
                    Chromosome::Autosome(chrom),
                    (i * 10000) as u32,
                    i as f64 * 0.003, // 0.3 cM = 0.003 M spacing
                    'A',
                    'G',
                ));
            }
        }

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        // Create genotypes with ancestry-specific patterns
        for snp_idx in 0..n_snps {
            // Source1 pattern (higher frequency)
            for sample_idx in 20..30 {
                data[[snp_idx, sample_idx]] = if snp_idx % 3 == 0 { 2 } else { 1 };
            }
            // Source2 pattern (lower frequency)
            for sample_idx in 30..40 {
                data[[snp_idx, sample_idx]] = if snp_idx % 3 == 0 { 0 } else { 1 };
            }
            // Test population (admixed - mix of patterns with some LD)
            for sample_idx in 0..20 {
                let base = if (snp_idx + sample_idx) % 2 == 0 {
                    data[[snp_idx, 20 + sample_idx % 10]] // Like source1
                } else {
                    data[[snp_idx, 30 + sample_idx % 10]] // Like source2
                };
                data[[snp_idx, sample_idx]] = base;
            }
        }

        GenotypeMatrix::new(data, samples, snps).unwrap()
    }

    #[test]
    fn test_dates_basic() {
        let geno = create_test_genotypes();

        let mut opts = DatesOptions::default();
        opts.ld_options.min_pairs_per_bin = 10; // Lower threshold for test
        opts.run_jackknife = false; // Skip jackknife for faster test

        let result = dates(&geno, None, "test", "source1", "source2", Some(opts));

        // Should complete without error
        assert!(result.is_ok(), "DATES should complete: {:?}", result);

        let result = result.unwrap();

        // Check basic properties
        assert!(result.n_snps > 0);
        assert!(result.generations > 0.0);
        assert!(result.years_bp > 0.0);
        assert!(!result.ld_bins.is_empty());
    }

    #[test]
    fn test_dates_missing_population() {
        let geno = create_test_genotypes();

        let result = dates(&geno, None, "nonexistent", "source1", "source2", None);
        assert!(result.is_err());
    }
}
