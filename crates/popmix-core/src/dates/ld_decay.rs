//! LD decay computation for DATES algorithm.
//!
//! Computes weighted linkage disequilibrium covariance binned by genetic distance.
//!
//! The DATES algorithm computes LD as follows:
//! 1. For each individual in the test population:
//!    - Compute residuals: res[snp] = genotype - predicted
//!    - predicted = α × p1 + (1-α) × p2 where α is estimated from regression
//!    - Compute weighted residuals: wr[snp] = res[snp] × (p1[snp] - p2[snp])
//! 2. For each SNP pair (x, y) at distance d:
//!    - Accumulate wr[x] × wr[y] into distance bins
//! 3. Compute covariance = E[y1×y2] - E[y1]×E[y2] over all pairs at distance d

use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::Chromosome;
use crate::Result;
use rayon::prelude::*;
use std::collections::HashMap;

/// Options for LD decay computation.
#[derive(Debug, Clone)]
pub struct LdOptions {
    /// Bin width in cM (default: 0.1)
    pub bin_width_cm: f64,
    /// Maximum genetic distance to compute (default: 50.0 cM)
    pub max_distance_cm: f64,
    /// Minimum number of SNP pairs per bin (default: 100)
    pub min_pairs_per_bin: usize,
    /// Minimum allele count in test population (default: 1)
    pub min_test_count: u32,
}

impl Default for LdOptions {
    fn default() -> Self {
        LdOptions {
            bin_width_cm: 0.1,
            max_distance_cm: 100.0, // DATES default: maxdis = 1.0 Morgan = 100 cM
            min_pairs_per_bin: 100,
            min_test_count: 1,
        }
    }
}

/// Result for a single genetic distance bin.
#[derive(Debug, Clone)]
pub struct LdBin {
    /// Center of the distance bin in cM.
    pub distance_cm: f64,
    /// Weighted LD covariance.
    pub weighted_ld: f64,
    /// Z-score for the weighted LD.
    pub z_score: f64,
    /// Correlation coefficient.
    pub correlation: f64,
    /// Number of SNP pairs in this bin.
    pub n_pairs: usize,
    /// Standard error of weighted LD.
    pub std_err: f64,
}

/// Correlation accumulator (matching DATES CORR structure).
#[derive(Debug, Clone, Default)]
struct CorrBin {
    s0: f64,  // count
    s1: f64,  // sum of y1
    s2: f64,  // sum of y2
    s11: f64, // sum of y1²
    s12: f64, // sum of y1×y2
    s22: f64, // sum of y2²
}

impl CorrBin {
    fn add(&mut self, y1: f64, y2: f64) {
        self.s0 += 1.0;
        self.s1 += y1;
        self.s2 += y2;
        self.s11 += y1 * y1;
        self.s12 += y1 * y2;
        self.s22 += y2 * y2;
    }

    fn merge(&mut self, other: &CorrBin) {
        self.s0 += other.s0;
        self.s1 += other.s1;
        self.s2 += other.s2;
        self.s11 += other.s11;
        self.s12 += other.s12;
        self.s22 += other.s22;
    }

    /// Compute covariance: E[y1×y2] - E[y1]×E[y2]
    fn covariance(&self) -> f64 {
        if self.s0 < 1.0 {
            return f64::NAN;
        }
        let n = self.s0;
        let m1 = self.s1 / n;
        let m2 = self.s2 / n;
        (self.s12 / n) - m1 * m2
    }

    /// Compute v12 without mean subtraction (DATES mode=1 behavior).
    /// This is just S12/n, matching DATES's calccorr with mode=1.
    fn v12_no_mean(&self) -> f64 {
        if self.s0 < 1.0 {
            return f64::NAN;
        }
        self.s12 / self.s0
    }

    /// Compute variances for z-score calculation
    fn variances(&self) -> (f64, f64) {
        if self.s0 < 1.0 {
            return (f64::NAN, f64::NAN);
        }
        let n = self.s0;
        let m1 = self.s1 / n;
        let m2 = self.s2 / n;
        let v1 = (self.s11 / n) - m1 * m1;
        let v2 = (self.s22 / n) - m2 * m2;
        (v1, v2)
    }

    /// Compute correlation coefficient (with mean subtraction)
    fn correlation(&self) -> f64 {
        let (v1, v2) = self.variances();
        let cov = self.covariance();
        if v1 > 0.0 && v2 > 0.0 {
            cov / (v1 * v2).sqrt()
        } else {
            f64::NAN
        }
    }

    /// Compute correlation WITHOUT mean subtraction (DATES mode=1).
    /// This is the formula DATES uses for output column 4.
    /// corr = (S12/n) / sqrt((S11/n) * (S22/n))
    fn correlation_no_mean(&self) -> f64 {
        if self.s0 < 1.0 {
            return f64::NAN;
        }
        let n = self.s0;
        let v11 = self.s11 / n;
        let v12 = self.s12 / n;
        let v22 = self.s22 / n;
        if v11 > 0.0 && v22 > 0.0 {
            v12 / (v11 * v22).sqrt()
        } else {
            f64::NAN
        }
    }
}

/// Compute weighted LD decay by genetic distance using the DATES algorithm.
///
/// # Arguments
/// * `geno` - Genotype matrix (needed for individual genotypes)
/// * `af` - Allele frequencies (for source population frequencies)
/// * `test_pop` - Test (admixed) population
/// * `weights` - Pre-computed weights for each SNP (p_source1 - p_source2)
/// * `source1_freqs` - Allele frequencies in source1 population
/// * `source2_freqs` - Allele frequencies in source2 population
/// * `options` - Computation options
///
/// # Returns
/// Vector of LdBin results, one per distance bin
pub fn compute_ld_decay(
    geno: &GenotypeMatrix,
    _af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    options: Option<LdOptions>,
) -> Result<Vec<LdBin>> {
    let opts = options.unwrap_or_default();

    // Get indices of test population samples
    let test_pop_id = crate::genotypes::PopulationId::new(test_pop);
    let test_indices = geno.population_indices(&test_pop_id);

    if test_indices.is_empty() {
        return Err(crate::error::PopmixError::PopulationNotFound(
            test_pop.to_string(),
        ));
    }

    // Group SNPs by chromosome
    // DATES filters out SNPs where |weight| < 0.001 (frequency difference < 0.1%)
    // Also filters SNPs where mean frequency variance is too low
    let snps = geno.snps();
    let mut snps_by_chrom: HashMap<Chromosome, Vec<usize>> = HashMap::new();
    for (idx, snp) in snps.iter().enumerate() {
        let w = weights[idx];
        // DATES filtering: |freq_diff| >= 0.001
        if !snp.is_autosomal() || !w.is_finite() || w.abs() < 0.001 {
            continue;
        }
        // DATES also filters by mean frequency variance: p*(1-p) >= 0.001
        let p1 = source1_freqs[idx];
        let p2 = source2_freqs[idx];
        let mean_freq = (p1 + p2) * 0.5;
        let variance = mean_freq * (1.0 - mean_freq);
        if variance < 0.001 {
            continue;
        }
        snps_by_chrom
            .entry(snp.chromosome.clone())
            .or_default()
            .push(idx);
    }

    // Number of bins
    let n_bins = (opts.max_distance_cm / opts.bin_width_cm).ceil() as usize;

    // Process each individual and aggregate results
    let all_bins: Vec<Vec<CorrBin>> = test_indices
        .par_iter()
        .map(|&sample_idx| {
            compute_individual_ld(
                geno,
                snps,
                &snps_by_chrom,
                sample_idx,
                weights,
                source1_freqs,
                source2_freqs,
                &opts,
                n_bins,
            )
        })
        .collect();

    // Merge results from all individuals
    let mut merged_bins: Vec<CorrBin> = vec![CorrBin::default(); n_bins];
    for ind_bins in all_bins {
        for (i, bin) in ind_bins.into_iter().enumerate() {
            merged_bins[i].merge(&bin);
        }
    }

    // Compute final statistics per bin
    // DATES uses bin center = (k+1)*binsize, so bin 0 is centered at 1*binsize
    let mut results = Vec::with_capacity(n_bins);
    for (bin_idx, bin) in merged_bins.iter().enumerate() {
        let count = bin.s0 as usize;
        if count >= opts.min_pairs_per_bin {
            // DATES convention: bin center = (bin_idx + 1) * bin_width
            let distance_cm = (bin_idx as f64 + 1.0) * opts.bin_width_cm;
            // Use v12 without mean subtraction to match DATES mode=1
            let weighted_ld = bin.v12_no_mean();
            // DATES mode=1: correlation without mean subtraction
            let correlation = bin.correlation_no_mean();

            // Z-score approximation using DATES formula
            // DATES: Z = sqrt(n) * correlation
            let z_score = if count > 0 && correlation.is_finite() {
                (count as f64).sqrt() * correlation
            } else {
                f64::NAN
            };

            // Regression coefficient: S12 / S11
            let std_err = if bin.s11 > 1e-20 {
                bin.s12 / bin.s11
            } else {
                f64::NAN
            };

            results.push(LdBin {
                distance_cm,
                weighted_ld,
                z_score,
                correlation,
                n_pairs: count,
                std_err,
            });
        }
    }

    Ok(results)
}

/// Compute LD for a single individual across all chromosomes.
fn compute_individual_ld(
    geno: &GenotypeMatrix,
    snps: &[crate::genotypes::Snp],
    snps_by_chrom: &HashMap<Chromosome, Vec<usize>>,
    sample_idx: usize,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    opts: &LdOptions,
    n_bins: usize,
) -> Vec<CorrBin> {
    let mut bins: Vec<CorrBin> = vec![CorrBin::default(); n_bins];

    // First pass: collect valid genotypes and compute regression coefficient α
    // α = dot(g - p2, p1 - p2) / dot(p1 - p2, p1 - p2)
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut valid_snps: Vec<(usize, f64)> = Vec::new(); // (snp_idx, genotype/2)

    for chrom_snps in snps_by_chrom.values() {
        for &snp_idx in chrom_snps {
            let g = geno.snp_raw(snp_idx)[sample_idx];
            if g <= 2 {
                let w = weights[snp_idx]; // p1 - p2
                // Skip SNPs with invalid weights (filtered out by freq diff or MAF checks)
                if !w.is_finite() {
                    continue;
                }

                let geno_val = g as f64 / 2.0; // Scale to [0, 1]
                let p1 = source1_freqs[snp_idx];
                let p2 = source2_freqs[snp_idx];

                // For regression: (geno - p2) regressed on (p1 - p2)
                numerator += (geno_val - p2) * w;
                denominator += w * w;
                valid_snps.push((snp_idx, geno_val));
            }
        }
    }

    // Compute mixture coefficient
    let alpha = if denominator > 1e-15 {
        numerator / denominator
    } else {
        0.5 // Default to 50% if no valid weights
    };

    // Second pass: compute weighted residuals
    // res = genotype - (α * p1 + (1-α) * p2)
    // weighted_res = res * weight
    let mut weighted_res: HashMap<usize, f64> = HashMap::new();

    for (snp_idx, geno_val) in valid_snps {
        let p1 = source1_freqs[snp_idx];
        let p2 = source2_freqs[snp_idx];
        let prediction = alpha * p1 + (1.0 - alpha) * p2;
        let residual = geno_val - prediction;
        let wr = residual * weights[snp_idx];
        weighted_res.insert(snp_idx, wr);
    }

    // Compute pairwise LD for each chromosome
    for (_, chrom_snps) in snps_by_chrom {
        // Sort SNPs by position
        let mut sorted_snps: Vec<usize> = chrom_snps
            .iter()
            .filter(|&&idx| weighted_res.contains_key(&idx))
            .copied()
            .collect();
        sorted_snps.sort_by(|&a, &b| {
            snps[a]
                .position_morgans
                .partial_cmp(&snps[b].position_morgans)
                .unwrap()
        });

        // Compute pairwise weighted LD
        // Convert options to Morgans to match DATES exactly
        // DATES works in Morgans: binsize = 0.001 M (= 0.1 cM), maxdis = 1.0 M (= 100 cM)
        let binsize_morgans = opts.bin_width_cm / 100.0;
        let max_distance_morgans = opts.max_distance_cm / 100.0;

        for (i, &snp_i) in sorted_snps.iter().enumerate() {
            let pos_i = snps[snp_i].position_morgans;
            let wr_i = weighted_res[&snp_i];

            for &snp_j in sorted_snps.iter().skip(i + 1) {
                let pos_j = snps[snp_j].position_morgans;
                let dis = pos_j - pos_i;

                // DATES: if(dis >= maxdis) break
                if dis >= max_distance_morgans {
                    break;
                }

                // DATES binning: s = (int)(dis/binsize)
                // Positions are now stored in Morgans, matching DATES exactly
                let bin_idx = (dis / binsize_morgans) as usize;
                if bin_idx >= n_bins {
                    continue;
                }

                let wr_j = weighted_res[&snp_j];
                bins[bin_idx].add(wr_i, wr_j);
            }
        }
    }

    bins
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{PopulationId, Sample, Snp};
    use ndarray::Array2;

    fn create_test_data() -> (GenotypeMatrix, AlleleFrequencies, Vec<f64>) {
        // Create 3 populations: test, source1, source2
        // with 50 SNPs on chromosome 1
        let mut samples = Vec::new();
        for i in 0..10 {
            samples.push(Sample::new(
                format!("test_{}", i),
                PopulationId::new("test"),
            ));
        }
        for i in 0..5 {
            samples.push(Sample::new(
                format!("src1_{}", i),
                PopulationId::new("source1"),
            ));
        }
        for i in 0..5 {
            samples.push(Sample::new(
                format!("src2_{}", i),
                PopulationId::new("source2"),
            ));
        }

        let snps: Vec<Snp> = (0..50)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 10000) as u32,
                    i as f64 * 0.005, // 0.005 M spacing (= 0.5 cM)
                    'A',
                    'G',
                )
            })
            .collect();

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        // Create genotypes with some LD structure
        for snp_idx in 0..n_snps {
            for sample_idx in 0..n_samples {
                let base = ((snp_idx + sample_idx) % 3) as u8;
                data[[snp_idx, sample_idx]] = base;
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);

        // Create weights (source1 - source2 frequencies)
        let idx1 = af.pop_index("source1").unwrap();
        let idx2 = af.pop_index("source2").unwrap();
        let weights: Vec<f64> = (0..n_snps)
            .map(|snp_idx| af.frequencies[[snp_idx, idx1]] - af.frequencies[[snp_idx, idx2]])
            .collect();

        (geno, af, weights)
    }

    #[test]
    fn test_ld_decay_basic() {
        let (geno, af, weights) = create_test_data();

        // Get source frequencies
        let idx1 = af.pop_index("source1").unwrap();
        let idx2 = af.pop_index("source2").unwrap();
        let n_snps = geno.n_snps();
        let source1_freqs: Vec<f64> = (0..n_snps)
            .map(|i| af.frequencies[[i, idx1]])
            .collect();
        let source2_freqs: Vec<f64> = (0..n_snps)
            .map(|i| af.frequencies[[i, idx2]])
            .collect();

        let opts = LdOptions {
            bin_width_cm: 1.0,
            max_distance_cm: 25.0,
            min_pairs_per_bin: 1,
            min_test_count: 1,
        };

        let result = compute_ld_decay(
            &geno, &af, "test", &weights, &source1_freqs, &source2_freqs, Some(opts)
        ).unwrap();

        // Should have some bins
        assert!(!result.is_empty(), "Should have LD decay bins");

        // All bins should have valid data
        for bin in &result {
            assert!(bin.distance_cm > 0.0);
            assert!(bin.n_pairs > 0);
        }
    }

    #[test]
    fn test_ld_decay_missing_population() {
        let (geno, af, weights) = create_test_data();

        let n_snps = geno.n_snps();
        let source1_freqs: Vec<f64> = vec![0.5; n_snps];
        let source2_freqs: Vec<f64> = vec![0.5; n_snps];

        let result = compute_ld_decay(
            &geno, &af, "nonexistent", &weights, &source1_freqs, &source2_freqs, None
        );
        assert!(result.is_err());
    }
}
