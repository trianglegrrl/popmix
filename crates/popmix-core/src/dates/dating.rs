//! Main DATES algorithm implementation.
//!
//! Estimates admixture timing from LD decay patterns.

use super::exponential_fit::{fit_exponential, jackknife_decay_rate, DecayFit, FitOptions};
use super::ld_decay::{compute_ld_decay, LdBin, LdOptions};
use super::weights::{compute_weights, WeightResult};
use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::{Chromosome, PopulationId};
use crate::Result;
use std::collections::HashMap;

/// Options for DATES algorithm.
#[derive(Debug, Clone)]
pub struct DatesOptions {
    /// Options for LD computation.
    pub ld_options: LdOptions,
    /// Options for exponential fitting.
    pub fit_options: FitOptions,
    /// Generation time in years (for converting to years BP). Default: 29.0
    pub generation_time: f64,
    /// Minimum allele count in source populations for weight computation.
    pub min_source_count: u32,
    /// Run jackknife for standard error estimation. Default: true
    pub run_jackknife: bool,
}

impl Default for DatesOptions {
    fn default() -> Self {
        DatesOptions {
            ld_options: LdOptions::default(),
            fit_options: FitOptions::default(),
            generation_time: 29.0,
            min_source_count: 4,
            run_jackknife: true,
        }
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
    /// Estimated number of generations since admixture.
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

    // Compute allele frequencies if not provided
    let owned_af;
    let af = match af {
        Some(af) => af,
        None => {
            let pops = vec![test_pop, source1, source2];
            owned_af =
                AlleleFrequencies::from_genotypes_for_populations(geno, Some(&pops));
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
    let source1_idx = af.pop_index(source1).ok_or_else(|| {
        crate::error::PopmixError::PopulationNotFound(source1.to_string())
    })?;
    let source2_idx = af.pop_index(source2).ok_or_else(|| {
        crate::error::PopmixError::PopulationNotFound(source2.to_string())
    })?;

    let source1_freqs: Vec<f64> = (0..geno.n_snps())
        .map(|i| af.frequencies[[i, source1_idx]])
        .collect();
    let source2_freqs: Vec<f64> = (0..geno.n_snps())
        .map(|i| af.frequencies[[i, source2_idx]])
        .collect();

    // Step 2: Compute LD decay
    let ld_bins = compute_ld_decay(
        geno,
        af,
        test_pop,
        &weights.weights,
        &source1_freqs,
        &source2_freqs,
        Some(opts.ld_options.clone()),
    )?;

    if ld_bins.is_empty() {
        return Err(crate::error::PopmixError::InsufficientData {
            context: "No LD bins computed".to_string(),
            required: 10,
            actual: 0,
        });
    }

    // Step 3: Fit exponential decay
    let mut fit = fit_exponential(&ld_bins, Some(opts.fit_options.clone()))?;

    // Step 4: Jackknife for standard error (if enabled)
    let (generations, generations_se) = if opts.run_jackknife {
        // Group LD bins by chromosome
        let ld_bins_by_chrom = compute_ld_by_chromosome(
            geno,
            af,
            test_pop,
            &weights.weights,
            Some(opts.ld_options.clone()),
        )?;

        match jackknife_decay_rate(&ld_bins_by_chrom, Some(opts.fit_options.clone())) {
            Ok((mean, se)) => {
                fit.decay_rate_se = se;
                (mean, se)
            }
            Err(_) => (fit.decay_rate, f64::NAN),
        }
    } else {
        (fit.decay_rate, f64::NAN)
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
        ld_bins,
        fit,
        n_snps: weights.n_valid,
        weights,
        n_chromosomes,
    })
}

/// Compute LD decay separately for each chromosome.
fn compute_ld_by_chromosome(
    geno: &GenotypeMatrix,
    _af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    options: Option<LdOptions>,
) -> Result<Vec<Vec<LdBin>>> {
    let opts = options.unwrap_or_default();

    // Get indices of test population samples
    let test_pop_id = PopulationId::new(test_pop);
    let test_indices = geno.population_indices(&test_pop_id);

    if test_indices.is_empty() {
        return Err(crate::error::PopmixError::PopulationNotFound(
            test_pop.to_string(),
        ));
    }

    // Group SNPs by chromosome
    let snps = geno.snps();
    let mut snps_by_chrom: HashMap<Chromosome, Vec<usize>> = HashMap::new();
    for (idx, snp) in snps.iter().enumerate() {
        if snp.is_autosomal() && weights[idx].is_finite() {
            snps_by_chrom
                .entry(snp.chromosome.clone())
                .or_default()
                .push(idx);
        }
    }

    // Sort chromosomes
    let mut chroms: Vec<_> = snps_by_chrom.keys().cloned().collect();
    chroms.sort();

    // Compute LD for each chromosome
    let mut results = Vec::with_capacity(chroms.len());

    for chrom in chroms {
        let chrom_snps = &snps_by_chrom[&chrom];
        if chrom_snps.len() < 10 {
            continue;
        }

        // Create a filtered genotype view for this chromosome
        // For now, just compute all LD and filter by chromosome
        // This is less efficient but simpler
        let n_bins = (opts.max_distance_cm / opts.bin_width_cm).ceil() as usize;
        let chrom_bins = compute_chromosome_ld(
            geno,
            snps,
            chrom_snps,
            &test_indices,
            weights,
            &opts,
            n_bins,
        );

        if !chrom_bins.is_empty() {
            results.push(chrom_bins);
        }
    }

    Ok(results)
}

/// Compute LD for a single chromosome and return as LdBin vector.
fn compute_chromosome_ld(
    geno: &GenotypeMatrix,
    snps: &[crate::genotypes::Snp],
    chrom_snps: &[usize],
    test_indices: &[usize],
    weights: &[f64],
    opts: &LdOptions,
    n_bins: usize,
) -> Vec<LdBin> {
    // Initialize bins
    let mut bin_sums = vec![0.0_f64; n_bins];
    let mut bin_sum_sq = vec![0.0_f64; n_bins];
    let mut bin_weight_sums = vec![0.0_f64; n_bins];
    let mut bin_counts = vec![0_usize; n_bins];

    // Pre-compute genotype means
    let mut snp_data: HashMap<usize, (f64, Vec<f64>)> = HashMap::new();
    for &snp_idx in chrom_snps {
        let (mean, genos) = compute_snp_stats_for_chrom(geno, snp_idx, test_indices);
        if mean.is_finite() {
            snp_data.insert(snp_idx, (mean, genos));
        }
    }

    // Sort SNPs by position
    let mut sorted_snps: Vec<usize> = chrom_snps
        .iter()
        .filter(|&&idx| snp_data.contains_key(&idx))
        .copied()
        .collect();
    sorted_snps.sort_by(|&a, &b| {
        snps[a]
            .position_morgans
            .partial_cmp(&snps[b].position_morgans)
            .unwrap()
    });

    // Convert options to Morgans to match DATES exactly
    let binsize_morgans = opts.bin_width_cm / 100.0;
    let max_distance_morgans = opts.max_distance_cm / 100.0;

    // Compute pairwise LD
    for (i, &snp_i) in sorted_snps.iter().enumerate() {
        let pos_i = snps[snp_i].position_morgans;
        let (mean_i, genos_i) = &snp_data[&snp_i];
        let w_i = weights[snp_i];

        for &snp_j in sorted_snps.iter().skip(i + 1) {
            let pos_j = snps[snp_j].position_morgans;
            let dis = pos_j - pos_i;

            if dis >= max_distance_morgans {
                break;
            }

            let bin_idx = (dis / binsize_morgans) as usize;
            if bin_idx >= n_bins {
                continue;
            }

            let (mean_j, genos_j) = &snp_data[&snp_j];
            let w_j = weights[snp_j];

            let cov = compute_cov_for_chrom(genos_i, genos_j, *mean_i, *mean_j);
            let weighted_ld = cov * w_i * w_j;
            let weight_product = w_i * w_j;

            bin_sums[bin_idx] += weighted_ld;
            bin_sum_sq[bin_idx] += weighted_ld * weighted_ld;
            bin_weight_sums[bin_idx] += weight_product.abs();
            bin_counts[bin_idx] += 1;
        }
    }

    // Convert to LdBin
    let mut results = Vec::new();
    for bin_idx in 0..n_bins {
        let count = bin_counts[bin_idx];
        if count >= opts.min_pairs_per_bin {
            let distance_cm = (bin_idx as f64 + 0.5) * opts.bin_width_cm;
            let mean_ld = bin_sums[bin_idx] / count as f64;
            let variance = (bin_sum_sq[bin_idx] / count as f64) - mean_ld * mean_ld;
            let std_err = if count > 1 {
                (variance / count as f64).sqrt()
            } else {
                f64::NAN
            };
            let z_score = if std_err > 0.0 && std_err.is_finite() {
                mean_ld / std_err
            } else {
                f64::NAN
            };
            let correlation = if bin_weight_sums[bin_idx].abs() > 1e-15 {
                bin_sums[bin_idx] / bin_weight_sums[bin_idx]
            } else {
                f64::NAN
            };

            results.push(LdBin {
                distance_cm,
                weighted_ld: mean_ld,
                z_score,
                correlation,
                n_pairs: count,
                std_err,
            });
        }
    }

    results
}

fn compute_snp_stats_for_chrom(
    geno: &GenotypeMatrix,
    snp_idx: usize,
    sample_indices: &[usize],
) -> (f64, Vec<f64>) {
    let mut sum = 0.0;
    let mut count = 0;
    let mut genos = Vec::with_capacity(sample_indices.len());

    for &sample_idx in sample_indices {
        let g = geno.snp_raw(snp_idx)[sample_idx];
        if g <= 2 {
            genos.push(g as f64);
            sum += g as f64;
            count += 1;
        } else {
            genos.push(f64::NAN);
        }
    }

    let mean = if count > 0 { sum / count as f64 } else { f64::NAN };
    (mean, genos)
}

fn compute_cov_for_chrom(genos_i: &[f64], genos_j: &[f64], mean_i: f64, mean_j: f64) -> f64 {
    let mut sum = 0.0;
    let mut count = 0;

    for (gi, gj) in genos_i.iter().zip(genos_j.iter()) {
        if gi.is_finite() && gj.is_finite() {
            sum += (gi - mean_i) * (gj - mean_j);
            count += 1;
        }
    }

    if count > 1 {
        sum / (count - 1) as f64
    } else {
        f64::NAN
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{Sample, Snp};
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
