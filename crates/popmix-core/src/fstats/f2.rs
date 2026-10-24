//! F2 statistic computation.
//!
//! f2(A,B) measures the genetic drift between two populations A and B.
//!
//! Formula: f2(A,B) = E[(pA - pB)²] - bias_correction
//!
//! Where:
//! - pA, pB are allele frequencies in populations A and B
//! - bias_correction = pA(1-pA)/(nA-1) + pB(1-pB)/(nB-1)
//! - nA, nB are the number of observed alleles (not individuals)
//!
//! The bias correction accounts for sampling variance and ensures
//! f2(A,A) = 0 for any population A.

use crate::genotype_matrix::AlleleFrequencies;

/// Result of an f2 computation.
#[derive(Debug, Clone)]
pub struct F2Result {
    /// The f2 estimate.
    pub estimate: f64,
    /// Standard error (from jackknife, if available).
    pub std_err: Option<f64>,
    /// Number of SNPs used in the computation.
    pub n_snps: usize,
}

/// Compute f2 for a single SNP between two populations.
///
/// Returns (f2_value, weight) where weight is 1.0 if both populations have data,
/// 0.0 otherwise (for missing data).
///
/// # Arguments
/// * `freq_a` - Allele frequency in population A (None if missing)
/// * `freq_b` - Allele frequency in population B (None if missing)
/// * `count_a` - Number of observed alleles in population A
/// * `count_b` - Number of observed alleles in population B
#[inline]
pub fn f2_single_snp(
    freq_a: Option<f64>,
    freq_b: Option<f64>,
    count_a: u32,
    count_b: u32,
) -> (f64, f64) {
    // Both populations need data
    let (pa, pb) = match (freq_a, freq_b) {
        (Some(a), Some(b)) => (a, b),
        _ => return (0.0, 0.0),
    };

    // Need at least 1 allele in each population
    if count_a < 1 || count_b < 1 {
        return (0.0, 0.0);
    }

    // Use max(1, count-1) for denominator, matching admixtools2 behavior
    // This allows single-sample populations (count=1) to work
    let denom_a = (count_a as f64 - 1.0).max(1.0);
    let denom_b = (count_b as f64 - 1.0).max(1.0);

    // Raw squared difference
    let diff = pa - pb;
    let raw_f2 = diff * diff;

    // Bias correction for sampling variance
    // This ensures f2(A,A) = 0
    let bias_a = pa * (1.0 - pa) / denom_a;
    let bias_b = pb * (1.0 - pb) / denom_b;

    let f2 = raw_f2 - bias_a - bias_b;

    (f2, 1.0)
}

/// Compute the f2 statistic between two populations across all SNPs.
///
/// Returns the weighted average f2 value.
///
/// # Arguments
/// * `af` - Allele frequencies for all populations
/// * `pop_a` - Name of population A
/// * `pop_b` - Name of population B
pub fn f2_estimate(af: &AlleleFrequencies, pop_a: &str, pop_b: &str) -> Option<F2Result> {
    let idx_a = af.pop_index(pop_a)?;
    let idx_b = af.pop_index(pop_b)?;

    let mut sum_f2 = 0.0;
    let mut sum_weight = 0.0;
    let mut n_snps_used = 0usize;

    for snp_idx in 0..af.n_snps() {
        let freq_a = af.get(snp_idx, idx_a);
        let freq_b = af.get(snp_idx, idx_b);
        let count_a = af.get_count(snp_idx, idx_a);
        let count_b = af.get_count(snp_idx, idx_b);

        let (f2, weight) = f2_single_snp(freq_a, freq_b, count_a, count_b);

        if weight > 0.0 {
            sum_f2 += f2 * weight;
            sum_weight += weight;
            n_snps_used += 1;
        }
    }

    if sum_weight == 0.0 {
        return None;
    }

    Some(F2Result {
        estimate: sum_f2 / sum_weight,
        std_err: None, // Will be computed with jackknife
        n_snps: n_snps_used,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_f2_identical_frequencies() {
        // When frequencies are identical, f2 should be close to 0
        // (exactly 0 after bias correction if same population)
        let (f2, weight) = f2_single_snp(Some(0.5), Some(0.5), 100, 100);
        assert_eq!(weight, 1.0);
        // Raw f2 = 0, bias corrections cancel out symmetrically
        // bias = 0.5 * 0.5 / 99 = 0.00252525... for each
        // f2 = 0 - 0.00252525 - 0.00252525 = -0.00505...
        // This is expected - f2 can be slightly negative due to bias correction
        assert!(f2.abs() < 0.01);
    }

    #[test]
    fn test_f2_different_frequencies() {
        // pA = 0.8, pB = 0.2, each with 100 alleles
        let (f2, weight) = f2_single_snp(Some(0.8), Some(0.2), 100, 100);
        assert_eq!(weight, 1.0);

        // Raw f2 = (0.8 - 0.2)² = 0.36
        // bias_a = 0.8 * 0.2 / 99 = 0.001616...
        // bias_b = 0.2 * 0.8 / 99 = 0.001616...
        // f2 = 0.36 - 0.001616 - 0.001616 ≈ 0.3568
        assert_relative_eq!(f2, 0.36 - 2.0 * 0.16 / 99.0, epsilon = 1e-10);
    }

    #[test]
    fn test_f2_missing_data() {
        // Missing frequency in population A
        let (f2, weight) = f2_single_snp(None, Some(0.5), 0, 100);
        assert_eq!(weight, 0.0);
        assert_eq!(f2, 0.0);

        // Missing frequency in population B
        let (f2, weight) = f2_single_snp(Some(0.5), None, 100, 0);
        assert_eq!(weight, 0.0);
        assert_eq!(f2, 0.0);
    }

    #[test]
    fn test_f2_single_sample() {
        // Single-sample populations should work (matching admixtools2)
        // Uses max(1, count-1) = max(1, 0) = 1 for denominator
        let (f2, weight) = f2_single_snp(Some(0.5), Some(0.5), 1, 100);
        assert_eq!(weight, 1.0);
        // bias_a = 0.5 * 0.5 / 1 = 0.25
        // bias_b = 0.5 * 0.5 / 99 ≈ 0.00252525
        // f2 = 0 - 0.25 - 0.00252525 ≈ -0.2525
        assert!(f2 < 0.0); // Expected to be negative due to bias correction

        let (f2, weight) = f2_single_snp(Some(0.5), Some(0.5), 100, 1);
        assert_eq!(weight, 1.0);
        assert!(f2 < 0.0);
    }

    #[test]
    fn test_f2_fixed_alleles() {
        // Population A fixed for ref (0), Population B fixed for alt (1)
        let (f2, weight) = f2_single_snp(Some(0.0), Some(1.0), 100, 100);
        assert_eq!(weight, 1.0);

        // Raw f2 = 1.0
        // bias_a = 0 * 1 / 99 = 0
        // bias_b = 1 * 0 / 99 = 0
        // f2 = 1.0
        assert_relative_eq!(f2, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_f2_symmetry() {
        // f2(A,B) should equal f2(B,A)
        let (f2_ab, _) = f2_single_snp(Some(0.3), Some(0.7), 50, 80);
        let (f2_ba, _) = f2_single_snp(Some(0.7), Some(0.3), 80, 50);
        assert_relative_eq!(f2_ab, f2_ba, epsilon = 1e-10);
    }

    #[test]
    fn test_f2_estimate_with_allele_frequencies() {
        use crate::genotype_matrix::GenotypeMatrix;
        use crate::genotypes::{Chromosome, PopulationId, Sample, Snp};
        use ndarray::array;

        // Create a simple test case with 2 populations and 4 SNPs
        let samples = vec![
            Sample::new("IND1", PopulationId::new("PopA")),
            Sample::new("IND2", PopulationId::new("PopA")),
            Sample::new("IND3", PopulationId::new("PopB")),
            Sample::new("IND4", PopulationId::new("PopB")),
        ];

        let snps = vec![
            Snp::new("rs1", Chromosome::Autosome(1), 100, 0.01, 'A', 'G'),
            Snp::new("rs2", Chromosome::Autosome(1), 200, 0.02, 'A', 'G'),
            Snp::new("rs3", Chromosome::Autosome(1), 300, 0.03, 'A', 'G'),
            Snp::new("rs4", Chromosome::Autosome(1), 400, 0.04, 'A', 'G'),
        ];

        // Genotypes: 0=HomRef, 1=Het, 2=HomAlt, 3=Missing
        // PopA individuals: columns 0,1
        // PopB individuals: columns 2,3
        let data = array![
            [0, 0, 2, 2], // SNP1: PopA=0, PopB=1
            [2, 2, 0, 0], // SNP2: PopA=1, PopB=0
            [1, 1, 1, 1], // SNP3: both=0.5
            [0, 2, 0, 2], // SNP4: both=0.5
        ];

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);

        let result = f2_estimate(&af, "PopA", "PopB").unwrap();

        // With these genotypes:
        // SNP1: pA=0, pB=1, f2 = 1.0 - 0 - 0 = 1.0
        // SNP2: pA=1, pB=0, f2 = 1.0 - 0 - 0 = 1.0
        // SNP3: pA=0.5, pB=0.5, f2 ≈ 0 (with small bias correction)
        // SNP4: pA=0.5, pB=0.5, f2 ≈ 0

        // Average should be around 0.5 (two SNPs with f2≈1, two with f2≈0)
        assert!(result.estimate > 0.4 && result.estimate < 0.6);
        assert_eq!(result.n_snps, 4);
    }

    #[test]
    fn test_f2_estimate_population_not_found() {
        use crate::genotype_matrix::GenotypeMatrix;
        use crate::genotypes::{Chromosome, PopulationId, Sample, Snp};
        use ndarray::array;

        let samples = vec![Sample::new("IND1", PopulationId::new("PopA"))];
        let snps = vec![Snp::new(
            "rs1",
            Chromosome::Autosome(1),
            100,
            0.01,
            'A',
            'G',
        )];
        let data = array![[0]];

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);

        // Unknown population should return None
        assert!(f2_estimate(&af, "PopA", "PopX").is_none());
        assert!(f2_estimate(&af, "PopX", "PopA").is_none());
    }
}
