//! qpWave pairs: Pairwise cladality testing.
//!
//! This module implements pairwise qpWave testing to screen for cladality
//! (shared ancestry) between all pairs of populations.

use crate::fstats::F2Blocks;
use crate::genotypes::PopulationId;
use crate::qpadm::qpwave::{qpwave, QpWaveOptions};
use crate::Result;
use itertools::Itertools;

/// Result of qpWave test for a single population pair.
#[derive(Debug, Clone)]
pub struct QpWavePairResult {
    /// First population.
    pub pop_a: PopulationId,
    /// Second population.
    pub pop_b: PopulationId,
    /// P-value for rank 0 (cladality) test.
    /// High p-value suggests the populations form a clade.
    pub p_value: f64,
    /// Chi-squared statistic for rank 0 test.
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
}

/// Result of pairwise qpWave analysis.
#[derive(Debug, Clone)]
pub struct QpWavePairsResult {
    /// Populations tested.
    pub populations: Vec<PopulationId>,
    /// Outgroup populations.
    pub outgroups: Vec<PopulationId>,
    /// Results for all population pairs.
    pub pairs: Vec<QpWavePairResult>,
    /// Number of pairs tested.
    pub n_pairs: usize,
}

impl QpWavePairsResult {
    /// Get pairs that pass the cladality test (p > threshold).
    pub fn clades(&self, p_threshold: f64) -> Vec<&QpWavePairResult> {
        self.pairs
            .iter()
            .filter(|p| p.p_value > p_threshold)
            .collect()
    }

    /// Get the pair with highest p-value (most likely clade); pairs with no
    /// p-value (0 degrees of freedom) are not candidates.
    pub fn best_clade(&self) -> Option<&QpWavePairResult> {
        self.pairs
            .iter()
            .filter(|p| !p.p_value.is_nan())
            .max_by(|a, b| a.p_value.total_cmp(&b.p_value))
    }

    /// Get p-value matrix as a hashmap for easy lookup.
    pub fn p_value_matrix(&self) -> std::collections::HashMap<(String, String), f64> {
        self.pairs
            .iter()
            .map(|p| ((p.pop_a.to_string(), p.pop_b.to_string()), p.p_value))
            .collect()
    }
}

/// Run qpWave for all population pairs.
///
/// Tests cladality (rank 0) for all C(n, 2) pairs of populations.
/// A high p-value indicates the pair forms a clade relative to the outgroups.
///
/// # Arguments
/// * `f2_blocks` - Pre-computed F2 blocks
/// * `populations` - Populations to test pairwise
/// * `outgroups` - Outgroup populations (right populations)
/// * `options` - Optional qpWave settings
///
/// # Returns
/// QpWavePairsResult with all pairwise test results
pub fn qpwave_pairs(
    f2_blocks: &F2Blocks,
    populations: &[&str],
    outgroups: &[&str],
    options: Option<QpWaveOptions>,
) -> Result<QpWavePairsResult> {
    let opts = options.unwrap_or_default();

    let mut pairs = Vec::new();

    // Test all population pairs
    for combo in populations.iter().combinations(2) {
        let pop_a = combo[0];
        let pop_b = combo[1];

        // Run qpWave with left = [pop_a, pop_b]
        let left = vec![*pop_a, *pop_b];

        match qpwave(f2_blocks, &left, outgroups, Some(opts.clone())) {
            Ok(result) => {
                // Extract rank 0 test (cladality)
                if let Some(rank0) = result.rank_tests.first() {
                    pairs.push(QpWavePairResult {
                        pop_a: PopulationId::new(*pop_a),
                        pop_b: PopulationId::new(*pop_b),
                        p_value: rank0.p_value,
                        chi_sq: rank0.chi_sq,
                        df: rank0.df,
                    });
                }
            }
            Err(_) => {
                // Skip failed pairs (e.g., singular matrix)
                continue;
            }
        }
    }

    let n_pairs = pairs.len();

    Ok(QpWavePairsResult {
        populations: populations.iter().map(|&s| PopulationId::new(s)).collect(),
        outgroups: outgroups.iter().map(|&s| PopulationId::new(s)).collect(),
        pairs,
        n_pairs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::blocks::BlockDefinition;
    use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    /// Create test F2 blocks.
    fn create_test_f2_blocks() -> F2Blocks {
        let mut samples = Vec::new();
        for (pop_name, count) in [
            ("PopA", 10),
            ("PopB", 10),
            ("PopC", 10),
            ("Out1", 10),
            ("Out2", 10),
        ] {
            for i in 0..count {
                samples.push(Sample::new(
                    format!("{}_{}", pop_name, i),
                    PopulationId::new(pop_name),
                ));
            }
        }

        let snps: Vec<Snp> = (0..100)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.2,
                    'A',
                    'G',
                )
            })
            .collect();

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        for snp_idx in 0..n_snps {
            let base = (snp_idx * 7 % 100) as f64 / 100.0;
            for sample_idx in 0..n_samples {
                let pop_idx = sample_idx / 10;
                let freq = match pop_idx {
                    0 => 0.3 + base * 0.2,  // PopA
                    1 => 0.35 + base * 0.2, // PopB (similar to A)
                    2 => 0.6 + base * 0.2,  // PopC (different)
                    3 => 0.1 + base * 0.1,  // Out1
                    4 => 0.15 + base * 0.1, // Out2
                    _ => 0.5,
                };
                let threshold = freq * 2.0;
                let noise = ((sample_idx * 13 + snp_idx * 17) % 100) as f64 / 100.0 * 2.0;
                data[[snp_idx, sample_idx]] = if noise < threshold { 2 } else { 0 };
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);
        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_qpwave_pairs_basic() {
        let f2_blocks = create_test_f2_blocks();

        let result = qpwave_pairs(
            &f2_blocks,
            &["PopA", "PopB", "PopC"],
            &["Out1", "Out2"],
            None,
        )
        .unwrap();

        // Should test C(3,2) = 3 pairs
        assert_eq!(result.n_pairs, 3);
        assert_eq!(result.pairs.len(), 3);

        // Check all pairs have valid p-values
        for pair in &result.pairs {
            assert!(pair.p_value >= 0.0 && pair.p_value <= 1.0);
            assert!(pair.chi_sq >= 0.0);
        }
    }

    #[test]
    fn test_qpwave_pairs_helper_methods() {
        let f2_blocks = create_test_f2_blocks();

        let result = qpwave_pairs(
            &f2_blocks,
            &["PopA", "PopB", "PopC"],
            &["Out1", "Out2"],
            None,
        )
        .unwrap();

        // Test clades method
        let clades = result.clades(0.05);
        for clade in clades {
            assert!(clade.p_value > 0.05);
        }

        // Test best_clade method
        if let Some(best) = result.best_clade() {
            for pair in &result.pairs {
                assert!(best.p_value >= pair.p_value);
            }
        }

        // Test p_value_matrix
        let matrix = result.p_value_matrix();
        assert_eq!(matrix.len(), result.n_pairs);
    }
}
