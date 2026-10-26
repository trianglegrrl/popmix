//! SNP weights for DATES: the difference of the two source allele
//! frequencies, at the SNPs informative enough to use.
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/dates.md` (section 3), Narasimhan et al. 2019 (Science,
//! supplementary methods of DATES), DATES's parameter-file documentation and
//! DATES v753 run as a black box.

use crate::error::{PopmixError, Result};
use crate::genotype_matrix::AlleleFrequencies;

/// Smallest |p1 − p2| for a SNP to be used.
const MIN_WEIGHT: f64 = 0.001;
/// Smallest m(1 − m), m the mean source frequency, for a SNP to be used.
const MIN_MEAN_HET: f64 = 0.001;

/// SNP weights w = p1 − p2.
#[derive(Debug, Clone)]
pub struct WeightResult {
    /// One per SNP: p1 − p2 where the SNP is valid, NaN elsewhere.
    pub weights: Vec<f64>,
    /// Number of valid SNPs.
    pub n_valid: usize,
    /// Indices of the valid SNPs, increasing.
    pub valid_indices: Vec<usize>,
}

/// Whether a SNP with weight `w` and source frequencies `p1`, `p2` is
/// informative enough: `w` finite with |w| ≥ 0.001, and m(1 − m) ≥ 0.001
/// for the mean frequency m = (p1 + p2) / 2.
pub(crate) fn informative(w: f64, p1: f64, p2: f64) -> bool {
    if !w.is_finite() || w.abs() < MIN_WEIGHT {
        return false;
    }
    let m = (p1 + p2) * 0.5;
    m * (1.0 - m) >= MIN_MEAN_HET
}

/// Weights p1 − p2 from the frequencies of `source1` and `source2` in `af`.
///
/// A SNP is valid when each source has at least `min_count` observed alleles
/// (default 4), both frequencies are finite, |p1 − p2| ≥ 0.001 and the mean
/// frequency m has m(1 − m) ≥ 0.001.
pub fn compute_weights(
    af: &AlleleFrequencies,
    source1: &str,
    source2: &str,
    min_count: Option<u32>,
) -> Result<WeightResult> {
    let i1 = af
        .pop_index(source1)
        .ok_or_else(|| PopmixError::PopulationNotFound(source1.to_string()))?;
    let i2 = af
        .pop_index(source2)
        .ok_or_else(|| PopmixError::PopulationNotFound(source2.to_string()))?;
    let min_count = min_count.unwrap_or(4);

    let mut weights = vec![f64::NAN; af.n_snps()];
    let mut valid_indices = Vec::new();
    for (i, w_out) in weights.iter_mut().enumerate() {
        if af.counts[[i, i1]] < min_count || af.counts[[i, i2]] < min_count {
            continue;
        }
        let p1 = af.frequencies[[i, i1]];
        let p2 = af.frequencies[[i, i2]];
        if !p1.is_finite() || !p2.is_finite() {
            continue;
        }
        let w = p1 - p2;
        if informative(w, p1, p2) {
            *w_out = w;
            valid_indices.push(i);
        }
    }
    Ok(WeightResult {
        weights,
        n_valid: valid_indices.len(),
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
            assert!((w - 4.0 / 6.0).abs() < 0.01, "Expected ~0.666, got {}", w);
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
