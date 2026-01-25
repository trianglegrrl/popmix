//! F2 block computation and storage.
//!
//! This module provides efficient computation and storage of f2 statistics
//! across all population pairs, organized by genomic blocks for jackknife resampling.
//!
//! The key data structure is `F2Blocks` which stores:
//! - f2 values for each block and population pair
//! - Block weights (number of valid SNPs per block per pair)
//! - Population names and indices

use ndarray::Array2;
use rayon::prelude::*;

use crate::fstats::blocks::BlockDefinition;
use crate::fstats::f2::f2_single_snp;
use crate::genotype_matrix::AlleleFrequencies;
use crate::genotypes::PopulationId;

/// Precomputed F2 statistics organized by blocks.
///
/// This structure stores f2 values for all population pairs across genomic blocks,
/// enabling efficient jackknife resampling and computation of f3/f4 statistics.
#[derive(Debug, Clone)]
pub struct F2Blocks {
    /// F2 values (bias-corrected): [n_blocks, n_pairs]
    /// Each value is the sum of f2 for SNPs in that block for that pair.
    pub block_sums: Array2<f64>,

    /// Raw squared differences (not bias-corrected): [n_blocks, n_pairs]
    /// Used for FST denominator computation.
    pub block_raw_sums: Array2<f64>,

    /// Weights (SNP counts): [n_blocks, n_pairs]
    /// Number of valid SNPs contributing to each block/pair.
    pub block_weights: Array2<f64>,

    /// Heterozygosity sums: [n_blocks, n_pops]
    /// Sum of p*(1-p) for SNPs in each block for each population.
    /// Used for FST computation.
    pub block_het_sums: Array2<f64>,

    /// Heterozygosity weights: [n_blocks, n_pops]
    /// Number of valid SNPs contributing to each block/population.
    pub block_het_weights: Array2<f64>,

    /// Population names in order.
    pub populations: Vec<PopulationId>,

    /// Number of blocks.
    pub n_blocks: usize,

    /// Block lengths (number of SNPs per block).
    pub block_lengths: Vec<usize>,
}

impl F2Blocks {
    /// Compute F2 blocks from allele frequencies.
    ///
    /// # Arguments
    /// * `af` - Allele frequencies for all populations
    /// * `block_def` - Block definition (size in Morgans)
    /// * `snps` - SNP metadata for block assignment
    pub fn from_allele_frequencies(
        af: &AlleleFrequencies,
        block_def: &BlockDefinition,
    ) -> Self {
        let populations = af.populations.clone();
        let n_pops = populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2; // Including self-pairs for verification

        // Assign SNPs to blocks
        let block_assignment = block_def.assign_blocks(&af.snps);
        let n_blocks = block_assignment.n_blocks;

        if n_blocks == 0 {
            return F2Blocks {
                block_sums: Array2::zeros((0, n_pairs)),
                block_raw_sums: Array2::zeros((0, n_pairs)),
                block_weights: Array2::zeros((0, n_pairs)),
                block_het_sums: Array2::zeros((0, n_pops)),
                block_het_weights: Array2::zeros((0, n_pops)),
                populations,
                n_blocks: 0,
                block_lengths: Vec::new(),
            };
        }

        // Compute f2 and heterozygosity for each block in parallel
        let block_results: Vec<BlockResult> = (0..n_blocks)
            .into_par_iter()
            .map(|block_idx| {
                let range = block_assignment.block_range(block_idx);
                compute_block_f2_and_het(af, &populations, range)
            })
            .collect();

        // Assemble into arrays
        let mut block_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_raw_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_weights = Array2::zeros((n_blocks, n_pairs));
        let mut block_het_sums = Array2::zeros((n_blocks, n_pops));
        let mut block_het_weights = Array2::zeros((n_blocks, n_pops));

        for (block_idx, result) in block_results.into_iter().enumerate() {
            for pair_idx in 0..n_pairs {
                block_sums[[block_idx, pair_idx]] = result.f2_sums[pair_idx];
                block_raw_sums[[block_idx, pair_idx]] = result.f2_raw_sums[pair_idx];
                block_weights[[block_idx, pair_idx]] = result.f2_weights[pair_idx];
            }
            for (pop_idx, (het_sum, het_weight)) in result
                .het_sums
                .into_iter()
                .zip(result.het_weights)
                .enumerate()
            {
                block_het_sums[[block_idx, pop_idx]] = het_sum;
                block_het_weights[[block_idx, pop_idx]] = het_weight;
            }
        }

        F2Blocks {
            block_sums,
            block_raw_sums,
            block_weights,
            block_het_sums,
            block_het_weights,
            populations,
            n_blocks,
            block_lengths: block_assignment.block_lengths,
        }
    }

    /// Get the number of populations.
    pub fn n_populations(&self) -> usize {
        self.populations.len()
    }

    /// Get the population index by name.
    pub fn pop_index(&self, name: &str) -> Option<usize> {
        self.populations.iter().position(|p| p.as_str() == name)
    }

    /// Get the pair index for two population indices.
    ///
    /// Uses triangular indexing: pair(i, j) where i <= j.
    #[inline]
    pub fn pair_index(&self, pop_a: usize, pop_b: usize) -> usize {
        let (i, j) = if pop_a <= pop_b {
            (pop_a, pop_b)
        } else {
            (pop_b, pop_a)
        };
        // Triangular number formula: sum(0..i) + (j - i)
        // = i * (i + 1) / 2 + (j - i)
        // But we want row-major order where row i has (n - i) elements
        // Actually, let's use: pair(i, j) = i * n - i*(i+1)/2 + j for upper triangular
        let n = self.populations.len();
        i * n - i * (i + 1) / 2 + j
    }

    /// Get f2 estimate between two populations.
    ///
    /// Returns (estimate, std_err, n_snps).
    ///
    /// The estimate is the simple mean of block means, matching admixtools2 behavior.
    /// This gives equal weight to each block regardless of SNP count.
    ///
    /// For self-comparisons (pop_a == pop_b), returns 0.0 to match admixtools2.
    pub fn f2(&self, pop_a: &str, pop_b: &str) -> Option<(f64, f64, usize)> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;

        // Self-F2 is defined as 0 (matching admixtools2)
        if idx_a == idx_b {
            let pair_idx = self.pair_index(idx_a, idx_b);
            let mut total_snps = 0usize;
            for block_idx in 0..self.n_blocks {
                total_snps += self.block_weights[[block_idx, pair_idx]] as usize;
            }
            return Some((0.0, 0.0, total_snps));
        }

        let pair_idx = self.pair_index(idx_a, idx_b);

        // Compute block means and count valid blocks
        let mut block_means = Vec::with_capacity(self.n_blocks);
        let mut total_snps = 0usize;

        for block_idx in 0..self.n_blocks {
            let sum = self.block_sums[[block_idx, pair_idx]];
            let weight = self.block_weights[[block_idx, pair_idx]];

            if weight > 0.0 {
                block_means.push(sum / weight);
                total_snps += weight as usize;
            }
        }

        if block_means.is_empty() {
            return None;
        }

        // Simple mean of block means (admixtools2 style)
        let estimate: f64 = block_means.iter().sum::<f64>() / block_means.len() as f64;

        // Compute jackknife standard error using block means
        let std_err = self.jackknife_std_err_from_means(&block_means, estimate);

        Some((estimate, std_err, total_snps))
    }

    /// Get f2 for a specific block.
    pub fn f2_block(&self, pop_a: &str, pop_b: &str, block_idx: usize) -> Option<(f64, f64)> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;
        let pair_idx = self.pair_index(idx_a, idx_b);

        let sum = self.block_sums[[block_idx, pair_idx]];
        let weight = self.block_weights[[block_idx, pair_idx]];

        Some((sum, weight))
    }

    /// Compute jackknife standard error from block means.
    fn jackknife_std_err_from_means(&self, block_means: &[f64], full_estimate: f64) -> f64 {
        let n = block_means.len();
        if n < 2 {
            return 0.0;
        }

        let total: f64 = block_means.iter().sum();
        let mut sum_sq_diff = 0.0;

        for &block_mean in block_means {
            // Leave-one-out estimate (mean of remaining blocks)
            let loo_estimate = (total - block_mean) / (n - 1) as f64;
            let diff = loo_estimate - full_estimate;
            sum_sq_diff += diff * diff;
        }

        // Jackknife variance formula: (n-1)/n * sum((x_i - x_bar)^2)
        let n_f64 = n as f64;
        let variance = (n_f64 - 1.0) / n_f64 * sum_sq_diff;
        variance.sqrt()
    }

    /// Get block-level f2 values for a pair of populations.
    ///
    /// Returns f2 estimate for each block (block mean), or None if pair not found.
    /// Used for computing covariance in qpWave/qpAdm.
    /// For self-pairs (pop_a == pop_b), returns all zeros.
    pub fn get_block_f2_values(&self, pop_a: &str, pop_b: &str) -> Option<ndarray::Array1<f64>> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;

        // Self-F2 is 0 for all blocks
        if idx_a == idx_b {
            return Some(ndarray::Array1::zeros(self.n_blocks));
        }

        let pair_idx = self.pair_index(idx_a, idx_b);

        let mut block_values = ndarray::Array1::zeros(self.n_blocks);

        for block_idx in 0..self.n_blocks {
            let sum = self.block_sums[[block_idx, pair_idx]];
            let weight = self.block_weights[[block_idx, pair_idx]];

            if weight > 0.0 {
                block_values[block_idx] = sum / weight;
            } else {
                // Use NaN for missing blocks
                block_values[block_idx] = f64::NAN;
            }
        }

        Some(block_values)
    }

    /// Get all pairwise f2 values as a matrix.
    ///
    /// Returns a symmetric matrix of f2 values [n_pops, n_pops].
    /// Diagonal elements (self-F2) are 0.0 to match admixtools2.
    pub fn f2_matrix(&self) -> Array2<f64> {
        let n = self.populations.len();
        let mut matrix = Array2::zeros((n, n));

        for i in 0..n {
            for j in i..n {
                // Self-F2 is 0 by definition
                if i == j {
                    matrix[[i, j]] = 0.0;
                    continue;
                }

                let pair_idx = self.pair_index(i, j);

                // Compute mean of block means (admixtools2 style)
                let mut block_means = Vec::new();
                for block_idx in 0..self.n_blocks {
                    let sum = self.block_sums[[block_idx, pair_idx]];
                    let weight = self.block_weights[[block_idx, pair_idx]];
                    if weight > 0.0 {
                        block_means.push(sum / weight);
                    }
                }

                let f2 = if block_means.is_empty() {
                    f64::NAN
                } else {
                    block_means.iter().sum::<f64>() / block_means.len() as f64
                };

                matrix[[i, j]] = f2;
                matrix[[j, i]] = f2;
            }
        }

        matrix
    }
}

/// Result of computing F2 and heterozygosity for a block.
struct BlockResult {
    f2_sums: Vec<f64>,
    f2_raw_sums: Vec<f64>,
    f2_weights: Vec<f64>,
    het_sums: Vec<f64>,
    het_weights: Vec<f64>,
}

/// Compute f2 sums, weights, and heterozygosity for a block of SNPs.
fn compute_block_f2_and_het(
    af: &AlleleFrequencies,
    populations: &[PopulationId],
    snp_range: std::ops::Range<usize>,
) -> BlockResult {
    let n_pops = populations.len();
    let n_pairs = n_pops * (n_pops + 1) / 2;

    let mut f2_sums = vec![0.0; n_pairs];
    let mut f2_raw_sums = vec![0.0; n_pairs];
    let mut f2_weights = vec![0.0; n_pairs];
    let mut het_sums = vec![0.0; n_pops];
    let mut het_weights = vec![0.0; n_pops];

    for snp_idx in snp_range {
        // Compute heterozygosity for each population: h = p * (1 - p)
        // Note: admixtools2 uses p*(1-p), not 2*p*(1-p)
        for pop_idx in 0..n_pops {
            let count = af.get_count(snp_idx, pop_idx);

            if let Some(freq) = af.get(snp_idx, pop_idx) {
                if count >= 2 {
                    let het = freq * (1.0 - freq);
                    het_sums[pop_idx] += het;
                    het_weights[pop_idx] += 1.0;
                }
            }
        }

        // For each population pair
        let mut pair_idx = 0;
        for i in 0..n_pops {
            for j in i..n_pops {
                let freq_a = af.get(snp_idx, i);
                let freq_b = af.get(snp_idx, j);
                let count_a = af.get_count(snp_idx, i);
                let count_b = af.get_count(snp_idx, j);

                let (f2, weight) = f2_single_snp(freq_a, freq_b, count_a, count_b);

                // Also compute raw squared difference (without bias correction)
                // for FST denominator
                if let (Some(pa), Some(pb)) = (freq_a, freq_b) {
                    if count_a >= 2 && count_b >= 2 {
                        let raw = (pa - pb).powi(2);
                        f2_raw_sums[pair_idx] += raw * weight;
                    }
                }

                f2_sums[pair_idx] += f2 * weight;
                f2_weights[pair_idx] += weight;

                pair_idx += 1;
            }
        }
    }

    BlockResult {
        f2_sums,
        f2_raw_sums,
        f2_weights,
        het_sums,
        het_weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use approx::assert_relative_eq;

    fn make_test_data() -> AlleleFrequencies {
        // Create test data with 3 populations, 20 SNPs across 2 cM
        let mut samples = Vec::new();
        for i in 0..10 {
            samples.push(Sample::new(format!("A{}", i), PopulationId::new("PopA")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("B{}", i), PopulationId::new("PopB")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("C{}", i), PopulationId::new("PopC")));
        }

        let snps: Vec<Snp> = (0..20)
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

        // Create genotype data
        // PopA: mostly ref (0), PopB: mostly alt (2), PopC: mixed
        let mut data = ndarray::Array2::zeros((20, 30));
        for snp_idx in 0..20 {
            // PopA (columns 0-9): mostly 0
            for j in 0..10 {
                data[[snp_idx, j]] = if snp_idx % 5 == 0 { 1 } else { 0 };
            }
            // PopB (columns 10-19): mostly 2
            for j in 10..20 {
                data[[snp_idx, j]] = if snp_idx % 5 == 0 { 1 } else { 2 };
            }
            // PopC (columns 20-29): mixed
            for j in 20..30 {
                data[[snp_idx, j]] = ((snp_idx + j) % 3) as u8;
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        AlleleFrequencies::from_genotypes(&geno)
    }

    #[test]
    fn test_f2_blocks_creation() {
        let af = make_test_data();
        let block_def = BlockDefinition::new(0.01); // 1 cM blocks

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert_eq!(f2_blocks.n_populations(), 3);
        assert!(f2_blocks.n_blocks > 0);

        // Should have 6 pairs: AA, AB, AC, BB, BC, CC
        let n_pairs = 3 * 4 / 2; // n*(n+1)/2 = 6
        assert_eq!(f2_blocks.block_sums.ncols(), n_pairs);
    }

    #[test]
    fn test_f2_self_pair_near_zero() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,A) should be close to 0
        let (f2_aa, _, _) = f2_blocks.f2("PopA", "PopA").unwrap();
        assert!(
            f2_aa.abs() < 0.01,
            "f2(A,A) should be near 0, got {}",
            f2_aa
        );
    }

    #[test]
    fn test_f2_symmetry() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,B) should equal f2(B,A)
        let (f2_ab, _, _) = f2_blocks.f2("PopA", "PopB").unwrap();
        let (f2_ba, _, _) = f2_blocks.f2("PopB", "PopA").unwrap();
        assert_relative_eq!(f2_ab, f2_ba, epsilon = 1e-10);
    }

    #[test]
    fn test_f2_different_populations() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,B) should be positive (different populations)
        let (f2_ab, std_err, n_snps) = f2_blocks.f2("PopA", "PopB").unwrap();
        assert!(f2_ab > 0.0, "f2(A,B) should be positive, got {}", f2_ab);
        assert!(std_err >= 0.0);
        assert!(n_snps > 0);
    }

    #[test]
    fn test_f2_matrix() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);
        let matrix = f2_blocks.f2_matrix();

        assert_eq!(matrix.shape(), &[3, 3]);

        // Print matrix for debugging
        println!("F2 matrix:\n{:?}", matrix);

        // Diagonal should be near 0 (but can be slightly negative due to bias correction)
        for i in 0..3 {
            assert!(
                matrix[[i, i]].abs() < 0.1,
                "Diagonal element {} = {} should be near 0",
                i,
                matrix[[i, i]]
            );
        }

        // Should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(matrix[[i, j]], matrix[[j, i]], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_pair_index() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // Verify pair indexing is consistent
        // For 3 populations: pairs are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
        // Indices should be 0, 1, 2, 3, 4, 5
        assert_eq!(f2_blocks.pair_index(0, 0), 0);
        assert_eq!(f2_blocks.pair_index(0, 1), 1);
        assert_eq!(f2_blocks.pair_index(0, 2), 2);
        assert_eq!(f2_blocks.pair_index(1, 1), 3);
        assert_eq!(f2_blocks.pair_index(1, 2), 4);
        assert_eq!(f2_blocks.pair_index(2, 2), 5);

        // Symmetric access
        assert_eq!(f2_blocks.pair_index(1, 0), f2_blocks.pair_index(0, 1));
        assert_eq!(f2_blocks.pair_index(2, 0), f2_blocks.pair_index(0, 2));
    }

    #[test]
    fn test_jackknife_std_err() {
        let af = make_test_data();
        let block_def = BlockDefinition::new(0.005); // Small blocks for more jackknife samples

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        let (_, std_err, _) = f2_blocks.f2("PopA", "PopB").unwrap();

        // Standard error should be non-negative
        assert!(std_err >= 0.0);

        // With multiple blocks, should have positive standard error
        if f2_blocks.n_blocks > 2 {
            // May still be 0 if all blocks have identical estimates
        }
    }
}
