//! Weight computation for DATES algorithm.
//!
//! Computes SNP weights based on allele frequency differences between source populations.

use crate::genotype_matrix::AlleleFrequencies;
use crate::error::PopmixError;
use crate::Result;

/// Result of weight computation for a single SNP.
#[derive(Debug, Clone)]
pub struct WeightResult {
    /// Weights for each SNP (p_source1 - p_source2).
    pub weights: Vec<f64>,
    /// Number of SNPs with valid weights.
    pub n_valid: usize,
    /// Indices of SNPs with valid weights.
    pub valid_indices: Vec<usize>,
}

/// Compute weights for each SNP based on allele frequency difference between sources.
///
/// The weight for SNP x is: w(x) = p_source1(x) - p_source2(x)
///
/// # Arguments
/// * `af` - Allele frequencies for all populations
/// * `source1` - First source population name
/// * `source2` - Second source population name
/// * `min_count` - Minimum allele count per source population (default: 4, i.e., 2 diploid individuals)
///
/// # Returns
/// WeightResult containing the weights for each SNP
pub fn compute_weights(
    af: &AlleleFrequencies,
    source1: &str,
    source2: &str,
    min_count: Option<u32>,
) -> Result<WeightResult> {
    let min_n = min_count.unwrap_or(4); // Default: at least 2 diploid individuals

    // Get population indices
    let idx1 = af.pop_index(source1).ok_or_else(|| {
        PopmixError::PopulationNotFound(source1.to_string())
    })?;
    let idx2 = af.pop_index(source2).ok_or_else(|| {
        PopmixError::PopulationNotFound(source2.to_string())
    })?;

    let n_snps = af.n_snps();
    let mut weights = vec![f64::NAN; n_snps];
    let mut valid_indices = Vec::with_capacity(n_snps);

    for snp_idx in 0..n_snps {
        // Check if both populations have sufficient data
        let count1 = af.counts[[snp_idx, idx1]];
        let count2 = af.counts[[snp_idx, idx2]];

        if count1 >= min_n && count2 >= min_n {
            let p1 = af.frequencies[[snp_idx, idx1]];
            let p2 = af.frequencies[[snp_idx, idx2]];

            // Both frequencies must be valid
            if p1.is_finite() && p2.is_finite() {
                let freq_diff = (p1 - p2).abs();

                // DATES filter: frequency difference must be >= 0.001
                if freq_diff < 0.001 {
                    continue;
                }

                // DATES filter: MAF check using mean frequency
                // y = (p1 + p2) / 2, then check y*(1-y) >= 0.001
                let mean_freq = (p1 + p2) * 0.5;
                let maf_term = mean_freq * (1.0 - mean_freq);
                if maf_term < 0.001 {
                    continue;
                }

                weights[snp_idx] = p1 - p2;
                valid_indices.push(snp_idx);
            }
        }
    }

    let n_valid = valid_indices.len();

    Ok(WeightResult {
        weights,
        n_valid,
        valid_indices,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, PopulationId, Sample, Snp};
    use ndarray::Array2;

    fn create_test_allele_frequencies() -> AlleleFrequencies {
        // Create simple test data with 2 populations and 10 SNPs
        let samples = vec![
            Sample::new("s1", PopulationId::new("pop1")),
            Sample::new("s2", PopulationId::new("pop1")),
            Sample::new("s3", PopulationId::new("pop1")),
            Sample::new("s4", PopulationId::new("pop2")),
            Sample::new("s5", PopulationId::new("pop2")),
            Sample::new("s6", PopulationId::new("pop2")),
        ];

        let snps: Vec<Snp> = (0..10)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.1,
                    'A',
                    'G',
                )
            })
            .collect();

        // Create genotype matrix with controlled values
        // pop1 has higher frequencies, pop2 has lower
        let mut data = Array2::zeros((10, 6));
        for snp_idx in 0..10 {
            // pop1 (samples 0-2): higher frequency
            data[[snp_idx, 0]] = 2;
            data[[snp_idx, 1]] = 2;
            data[[snp_idx, 2]] = 1;
            // pop2 (samples 3-5): lower frequency
            data[[snp_idx, 3]] = 0;
            data[[snp_idx, 4]] = 1;
            data[[snp_idx, 5]] = 0;
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        AlleleFrequencies::from_genotypes(&geno)
    }

    #[test]
    fn test_weight_computation() {
        let af = create_test_allele_frequencies();

        // Use min_count=2 so we need at least 1 diploid individual
        let result = compute_weights(&af, "pop1", "pop2", Some(2)).unwrap();

        // All 10 SNPs should have valid weights
        assert_eq!(result.n_valid, 10);
        assert_eq!(result.valid_indices.len(), 10);

        // Check that weights are p_pop1 - p_pop2
        // pop1: (2+2+1)/(3*2) = 5/6 = 0.833...
        // pop2: (0+1+0)/(3*2) = 1/6 = 0.166...
        // weight = 5/6 - 1/6 = 4/6 = 0.666...
        for &idx in &result.valid_indices {
            let w = result.weights[idx];
            assert!(w.is_finite());
            // pop1 has higher freq than pop2, so weights should be positive
            assert!(w > 0.0, "Weight at {} should be positive, got {}", idx, w);
            // Check approximate value
            assert!((w - 4.0/6.0).abs() < 0.01, "Expected ~0.666, got {}", w);
        }
    }

    #[test]
    fn test_weight_signs() {
        let af = create_test_allele_frequencies();

        // Swap source order - weights should be negated
        let result1 = compute_weights(&af, "pop1", "pop2", Some(2)).unwrap();
        let result2 = compute_weights(&af, "pop2", "pop1", Some(2)).unwrap();

        for idx in 0..result1.weights.len() {
            if result1.weights[idx].is_finite() {
                assert!(
                    (result1.weights[idx] + result2.weights[idx]).abs() < 1e-10,
                    "Swapping sources should negate weights"
                );
            }
        }
    }

    #[test]
    fn test_missing_population() {
        let af = create_test_allele_frequencies();

        let result = compute_weights(&af, "pop1", "nonexistent", Some(2));
        assert!(result.is_err());
    }
}
