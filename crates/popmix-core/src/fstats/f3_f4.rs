//! F3, F4, and FST statistics.
//!
//! These statistics are derived from F2 blocks:
//! - f3(A; B, C) = (f2(A,B) + f2(A,C) - f2(B,C)) / 2
//!   Tests if A is admixed between B and C (negative f3 indicates admixture)
//!
//! - f4(A,B; C,D) = (f2(A,D) + f2(B,C) - f2(A,C) - f2(B,D)) / 2
//!   Tests if there is gene flow between (A,B) and (C,D) branches
//!   Also known as D-statistic or ABBA-BABA test
//!
//! - FST(A,B) = f2(A,B) / (f2(A,B) + expected_het)
//!   Fixation index measuring population differentiation

use crate::fstats::f2_blocks::F2Blocks;
use crate::fstats::jackknife::{block_jackknife, z_score};

/// Result of an f3 computation.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct F3Result {
    /// The f3 estimate.
    pub estimate: f64,
    /// Standard error from jackknife.
    pub std_err: f64,
    /// Z-score (estimate / std_err), NaN when the SE is 0.
    pub z_score: f64,
    /// Number of blocks used.
    pub n_blocks: usize,
}

/// Result of an f4 computation.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct F4Result {
    /// The f4 estimate.
    pub estimate: f64,
    /// Standard error from jackknife.
    pub std_err: f64,
    /// Z-score (estimate / std_err), NaN when the SE is 0.
    pub z_score: f64,
    /// Number of blocks used.
    pub n_blocks: usize,
}

/// Result of an FST computation.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct FstResult {
    /// The FST estimate (Hudson's FST).
    pub estimate: f64,
    /// Standard error from jackknife.
    pub std_err: f64,
    /// Number of SNPs used.
    pub n_snps: usize,
}

impl F2Blocks {
    /// Compute f3(A; B, C) statistic.
    ///
    /// f3(A; B, C) = (f2(A,B) + f2(A,C) - f2(B,C)) / 2
    ///
    /// A negative f3 indicates that A is admixed between B and C.
    ///
    /// # Arguments
    /// * `pop_a` - Target population (potentially admixed)
    /// * `pop_b` - Source population 1
    /// * `pop_c` - Source population 2
    pub fn f3(&self, pop_a: &str, pop_b: &str, pop_c: &str) -> Option<F3Result> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;
        let idx_c = self.pop_index(pop_c)?;

        // f3 for each block; NA where a pair has no SNPs in the block
        let block_f3: Vec<f64> = (0..self.n_blocks)
            .map(|block| {
                let mean = |i, j| self.block_mean(block, i, j);
                (mean(idx_a, idx_b) + mean(idx_a, idx_c) - mean(idx_b, idx_c)) / 2.0
            })
            .collect();

        // Weighted block jackknife, as admixtools2 f3()
        let (estimate, std_err) = block_jackknife(&block_f3, &self.block_lengths_f64())?;
        let n_blocks = block_f3.iter().filter(|v| v.is_finite()).count();
        let z_score = z_score(estimate, std_err);

        Some(F3Result {
            estimate,
            std_err,
            z_score,
            n_blocks,
        })
    }

    /// Compute f4(A, B; C, D) statistic.
    ///
    /// f4(A,B; C,D) = (f2(A,D) + f2(B,C) - f2(A,C) - f2(B,D)) / 2
    ///
    /// A non-zero f4 indicates gene flow between the (A,B) and (C,D) branches.
    ///
    /// # Arguments
    /// * `pop_a`, `pop_b` - First pair of populations
    /// * `pop_c`, `pop_d` - Second pair of populations
    pub fn f4(&self, pop_a: &str, pop_b: &str, pop_c: &str, pop_d: &str) -> Option<F4Result> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;
        let idx_c = self.pop_index(pop_c)?;
        let idx_d = self.pop_index(pop_d)?;

        // f4 for each block; NA where a pair has no SNPs in the block. With
        // A = B or C = D it is 0 by construction: exactly 0, not rounding
        // noise whose "SE" would give a z in the thousands (ADMIXTOOLS2-0008).
        let zero = idx_a == idx_b || idx_c == idx_d;
        let block_f4: Vec<f64> = (0..self.n_blocks)
            .map(|block| {
                let mean = |i, j| self.block_mean(block, i, j);
                let f4 = (mean(idx_a, idx_d) + mean(idx_b, idx_c)
                    - mean(idx_a, idx_c)
                    - mean(idx_b, idx_d))
                    / 2.0;
                if zero && f4.is_finite() {
                    0.0
                } else {
                    f4
                }
            })
            .collect();

        // Weighted block jackknife, as admixtools2 f4()
        let (estimate, std_err) = block_jackknife(&block_f4, &self.block_lengths_f64())?;
        let n_blocks = block_f4.iter().filter(|v| v.is_finite()).count();
        let z_score = z_score(estimate, std_err);

        Some(F4Result {
            estimate,
            std_err,
            z_score,
            n_blocks,
        })
    }

    /// f4(A, B; C, D), under its old name.
    ///
    /// This is f4, not Patterson's D, which is normalized per SNP and needs
    /// the genotypes ([`crate::f4_from_genotypes`] with [`crate::F4Mode::D`]).
    pub fn d_stat(&self, pop_a: &str, pop_b: &str, pop_c: &str, pop_d: &str) -> Option<F4Result> {
        self.f4(pop_a, pop_b, pop_c, pop_d)
    }

    /// Compute FST (fixation index) between two populations.
    ///
    /// FST measures population differentiation on a scale of 0 to 1:
    /// - FST = 0: No differentiation (populations are identical)
    /// - FST = 1: Complete fixation (populations share no alleles)
    ///
    /// Uses the Reich lab FST estimator (matching admixtools2):
    /// For each block: FST_block = mean(f2) / mean(f2 + het_A + het_B)
    /// Overall: FST = mean(FST_block values)
    /// This is mean-of-ratios, matching admixtools2's approach.
    ///
    /// # Arguments
    /// * `pop_a` - First population
    /// * `pop_b` - Second population
    pub fn fst(&self, pop_a: &str, pop_b: &str) -> Option<FstResult> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;

        // Self-FST is 0 by definition
        if idx_a == idx_b {
            // Count total SNPs from heterozygosity weights
            let mut total_snps = 0usize;
            for block_idx in 0..self.n_blocks {
                total_snps += self.block_het_weights[[block_idx, idx_a]] as usize;
            }
            return Some(FstResult {
                estimate: 0.0,
                std_err: 0.0,
                n_snps: total_snps,
            });
        }

        let pair_idx = self.pair_index(idx_a, idx_b);

        // admixtools2 fst(): weighted block jackknife of per-block FST on its
        // own SNP set and blocks
        if let Some(fst_blocks) = &self.fst {
            let values = fst_blocks.block_values(pair_idx);
            let lengths: Vec<f64> = fst_blocks.block_lengths.iter().map(|&l| l as f64).collect();
            let (estimate, std_err) = block_jackknife(&values, &lengths)?;
            let n_snps = values
                .iter()
                .zip(&fst_blocks.block_lengths)
                .filter(|(v, _)| v.is_finite())
                .map(|(_, &l)| l)
                .sum();
            return Some(FstResult {
                estimate,
                std_err,
                n_snps,
            });
        }

        // F2 directories from storage versions 1-3 have no FST blocks: an
        // approximation from the f2 blocks (polymorphic SNPs only).
        // Compute per-block FST values (mean-of-ratios approach, matching admixtools2)
        // For each block: FST_block = f2 / (f2 + het_a + het_b)
        // where f2 is bias-corrected and the denominator simplifies to
        // raw_squared_diff + het_a + het_b (see admixtools2 mats_to_fstarr)
        let mut block_fst_values = Vec::with_capacity(self.n_blocks);
        let mut total_snps = 0usize;

        for block_idx in 0..self.n_blocks {
            let f2_sum = self.block_sums[[block_idx, pair_idx]];
            let f2_raw_sum = self.block_raw_sums[[block_idx, pair_idx]];
            let f2_weight = self.block_weights[[block_idx, pair_idx]];

            let het_a_sum = self.block_het_sums[[block_idx, idx_a]];
            let het_a_weight = self.block_het_weights[[block_idx, idx_a]];
            let het_b_sum = self.block_het_sums[[block_idx, idx_b]];
            let het_b_weight = self.block_het_weights[[block_idx, idx_b]];

            // Need valid data for f2 and both heterozygosities
            if f2_weight > 0.0 && het_a_weight > 0.0 && het_b_weight > 0.0 {
                let f2 = f2_sum / f2_weight;
                let f2_raw = f2_raw_sum / f2_weight;
                let het_a = het_a_sum / het_a_weight;
                let het_b = het_b_sum / het_b_weight;

                // Denominator: raw squared difference + het_a + het_b
                // (this is equivalent to f2 + het_ub_a + het_ub_b after simplification)
                let denom = f2_raw + het_a + het_b;

                if denom > 0.0 {
                    // Per-block FST ratio
                    block_fst_values.push(f2 / denom);
                    total_snps += f2_weight as usize;
                }
            }
        }

        if block_fst_values.is_empty() {
            return None;
        }

        // Mean-of-ratios: FST = mean(FST_block values)
        let n_blocks = block_fst_values.len();
        let estimate: f64 = block_fst_values.iter().sum::<f64>() / n_blocks as f64;

        // Jackknife standard error
        let std_err = jackknife_std_err(&block_fst_values, estimate);

        Some(FstResult {
            estimate,
            std_err,
            n_snps: total_snps,
        })
    }
}

/// Compute jackknife standard error from block estimates.
fn jackknife_std_err(block_values: &[f64], full_estimate: f64) -> f64 {
    let n = block_values.len();
    if n < 2 {
        return 0.0;
    }

    let total: f64 = block_values.iter().sum();

    let mut sum_sq_diff = 0.0;
    for value in block_values {
        // Leave-one-out estimate
        let loo_estimate = (total - value) / (n - 1) as f64;
        let diff = loo_estimate - full_estimate;
        sum_sq_diff += diff * diff;
    }

    let variance = ((n - 1) as f64 / n as f64) * sum_sq_diff;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::blocks::BlockDefinition;
    use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
    use crate::genotypes::{Chromosome, PopulationId, Sample, Snp};

    fn make_test_data_4pops() -> AlleleFrequencies {
        // Create 4 populations with distinct frequency patterns
        let mut samples = Vec::new();
        for i in 0..10 {
            samples.push(Sample::new(format!("A{i}"), PopulationId::new("PopA")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("B{i}"), PopulationId::new("PopB")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("C{i}"), PopulationId::new("PopC")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("D{i}"), PopulationId::new("PopD")));
        }

        let snps: Vec<Snp> = (0..40)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    f64::from(i) * 0.2, // 0.2 cM spacing
                    'A',
                    'G',
                )
            })
            .collect();

        // Create distinct genotype patterns
        // PopA: freq ~0.2, PopB: freq ~0.8, PopC: freq ~0.5, PopD: freq ~0.5
        let mut data = ndarray::Array2::zeros((40, 40));
        for snp_idx in 0..40 {
            // PopA (columns 0-9): mostly ref
            for j in 0..10 {
                data[[snp_idx, j]] = if (snp_idx + j) % 5 == 0 { 2 } else { 0 };
            }
            // PopB (columns 10-19): mostly alt
            for j in 10..20 {
                data[[snp_idx, j]] = if (snp_idx + j) % 5 == 0 { 0 } else { 2 };
            }
            // PopC (columns 20-29): mixed, similar to D
            for j in 20..30 {
                data[[snp_idx, j]] = ((snp_idx + j) % 3) as u8;
            }
            // PopD (columns 30-39): mixed, similar to C
            for j in 30..40 {
                data[[snp_idx, j]] = ((snp_idx + j + 1) % 3) as u8;
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        AlleleFrequencies::from_genotypes(&geno)
    }

    #[test]
    fn a_population_has_no_f2_with_itself() {
        // admixtools2 sets the f2 of a population with itself to 0, so
        // f3(A; A, B) = 0, f3(A; B, B) = f2(A, B), f4(A, B; A, C) = f3(A; B, C)
        // and f4(A, B; A, B) = f2(A, B)
        let af = make_test_data_4pops();
        let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::default());
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * b.abs().max(1e-12);
        assert_eq!(f2.f3("PopA", "PopA", "PopB").unwrap().estimate, 0.0);
        assert_eq!(f2.f3("PopA", "PopB", "PopA").unwrap().estimate, 0.0);
        let f2_ab = f2.f2("PopA", "PopB").unwrap().0;
        assert!(close(
            f2.f3("PopA", "PopB", "PopB").unwrap().estimate,
            f2_ab
        ));
        let f3 = f2.f3("PopA", "PopB", "PopC").unwrap().estimate;
        assert!(close(
            f2.f4("PopA", "PopB", "PopA", "PopC").unwrap().estimate,
            f3
        ));
        assert!(close(
            f2.f4("PopA", "PopB", "PopA", "PopB").unwrap().estimate,
            f2_ab
        ));
    }

    #[test]
    fn z_is_nan_without_an_se() {
        // f3(A; A, B) and f4(A, A; B, C) are 0 in every block: SE 0, and z
        // is undefined, not 0 (TASKS.md G8)
        let af = make_test_data_4pops();
        let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::default());
        let f3 = f2.f3("PopA", "PopA", "PopB").unwrap();
        assert_eq!(f3.std_err, 0.0);
        assert!(f3.z_score.is_nan());
        for (a, b, c, d) in [
            ("PopA", "PopA", "PopB", "PopC"),
            ("PopA", "PopB", "PopC", "PopC"),
        ] {
            let f4 = f2.f4(a, b, c, d).unwrap();
            assert_eq!(f4.estimate, 0.0);
            assert_eq!(f4.std_err, 0.0);
            assert!(f4.z_score.is_nan());
        }
    }

    #[test]
    fn test_f3_computation() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        let result = f2_blocks.f3("PopC", "PopA", "PopB").unwrap();

        // f3 should have a value and standard error
        assert!(result.n_blocks > 0);
        assert!(result.std_err >= 0.0);

        println!(
            "f3(C; A, B) = {:.6} ± {:.6} (z = {:.2})",
            result.estimate, result.std_err, result.z_score
        );
    }

    #[test]
    fn test_f3_symmetry_in_sources() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f3(A; B, C) should equal f3(A; C, B)
        let f3_abc = f2_blocks.f3("PopA", "PopB", "PopC").unwrap();
        let f3_acb = f2_blocks.f3("PopA", "PopC", "PopB").unwrap();

        assert!(
            (f3_abc.estimate - f3_acb.estimate).abs() < 1e-10,
            "f3 should be symmetric in source populations"
        );
    }

    #[test]
    fn test_f4_computation() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        let result = f2_blocks.f4("PopA", "PopB", "PopC", "PopD").unwrap();

        assert!(result.n_blocks > 0);
        assert!(result.std_err >= 0.0);

        println!(
            "f4(A, B; C, D) = {:.6} ± {:.6} (z = {:.2})",
            result.estimate, result.std_err, result.z_score
        );
    }

    #[test]
    fn test_f4_symmetries() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f4(A,B; C,D) = f4(B,A; D,C) (swap both pairs)
        let f4_abcd = f2_blocks.f4("PopA", "PopB", "PopC", "PopD").unwrap();
        let f4_badc = f2_blocks.f4("PopB", "PopA", "PopD", "PopC").unwrap();
        assert!(
            (f4_abcd.estimate - f4_badc.estimate).abs() < 1e-10,
            "f4(A,B;C,D) should equal f4(B,A;D,C)"
        );

        // f4(A,B; C,D) = -f4(B,A; C,D) (swap first pair only)
        let f4_bacd = f2_blocks.f4("PopB", "PopA", "PopC", "PopD").unwrap();
        assert!(
            (f4_abcd.estimate + f4_bacd.estimate).abs() < 1e-10,
            "f4(A,B;C,D) should equal -f4(B,A;C,D)"
        );

        // f4(A,B; C,D) = f4(C,D; A,B) (swap the pairs)
        let f4_cdab = f2_blocks.f4("PopC", "PopD", "PopA", "PopB").unwrap();
        assert!(
            (f4_abcd.estimate - f4_cdab.estimate).abs() < 1e-10,
            "f4(A,B;C,D) should equal f4(C,D;A,B)"
        );
    }

    #[test]
    fn test_f4_same_population_is_zero() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f4(A,A; C,D) should be 0 (or very close)
        let result = f2_blocks.f4("PopA", "PopA", "PopC", "PopD").unwrap();
        assert!(
            result.estimate.abs() < 0.01,
            "f4(A,A;C,D) should be ~0, got {}",
            result.estimate
        );
    }

    #[test]
    fn test_d_stat_equals_f4() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        let f4 = f2_blocks.f4("PopA", "PopB", "PopC", "PopD").unwrap();
        let d = f2_blocks.d_stat("PopA", "PopB", "PopC", "PopD").unwrap();

        assert!((f4.estimate - d.estimate).abs() < 1e-10);
        assert!((f4.std_err - d.std_err).abs() < 1e-10);
    }

    #[test]
    fn test_f3_unknown_population() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert!(f2_blocks.f3("PopA", "PopB", "Unknown").is_none());
        assert!(f2_blocks.f3("Unknown", "PopB", "PopC").is_none());
    }

    #[test]
    fn test_f4_unknown_population() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert!(f2_blocks.f4("PopA", "PopB", "PopC", "Unknown").is_none());
    }

    // =========================================================================
    // FST tests
    // =========================================================================

    #[test]
    fn test_fst_computation() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // FST between divergent populations should be positive
        let result = f2_blocks.fst("PopA", "PopB").unwrap();
        assert!(
            result.estimate > 0.0,
            "FST should be positive for divergent populations"
        );
        assert!(result.estimate <= 1.0, "FST should be <= 1.0");
        assert!(
            result.std_err >= 0.0,
            "Standard error should be non-negative"
        );

        println!("FST(A, B) = {:.4} ± {:.4}", result.estimate, result.std_err);
    }

    #[test]
    fn test_fst_symmetry() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // FST(A,B) should equal FST(B,A)
        let fst_ab = f2_blocks.fst("PopA", "PopB").unwrap();
        let fst_ba = f2_blocks.fst("PopB", "PopA").unwrap();

        assert!(
            (fst_ab.estimate - fst_ba.estimate).abs() < 1e-10,
            "FST should be symmetric"
        );
    }

    #[test]
    fn test_fst_same_population_is_zero() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // FST(A,A) should be 0 (or very close)
        let result = f2_blocks.fst("PopA", "PopA").unwrap();
        assert!(
            result.estimate.abs() < 0.01,
            "FST(A,A) should be ~0, got {}",
            result.estimate
        );
    }

    #[test]
    fn test_fst_divergent_vs_similar() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // PopA and PopB are very different (0.2 vs 0.8 frequency)
        // PopC and PopD are similar (both ~0.5)
        let fst_ab = f2_blocks.fst("PopA", "PopB").unwrap();
        let fst_cd = f2_blocks.fst("PopC", "PopD").unwrap();

        assert!(
            fst_ab.estimate > fst_cd.estimate,
            "FST(A,B) should be > FST(C,D) since A,B are more divergent"
        );
    }

    /// F2 blocks of A and B from storage versions 1-3 (no FST blocks):
    /// three blocks with everything FST needs, then one each without SNPs
    /// for the pair, for A's heterozygosity and for B's, and one whose
    /// denominator is 0. Per usable block, FST = f2 / (raw f2 + het A +
    /// het B): 0.1 / 0.65, 0.05 / 0.65 and 0.1 / 0.55.
    fn f2_blocks_without_fst() -> F2Blocks {
        // columns: pairs (A, A), (A, B), (B, B); A, B
        let pair = |v: [f64; 7]| {
            ndarray::Array2::from_shape_fn((7, 3), |(b, p)| if p == 1 { v[b] } else { 0.0 })
        };
        let pop = |a: [f64; 7], b: [f64; 7]| {
            ndarray::Array2::from_shape_fn((7, 2), |(k, p)| [a, b][p][k])
        };
        F2Blocks {
            block_sums: pair([0.2, 0.15, 0.4, 0.1, 0.1, 0.1, -0.1]),
            block_raw_sums: pair([0.3, 0.3, 0.8, 0.2, 0.2, 0.2, 0.0]),
            block_weights: pair([2.0, 3.0, 4.0, 0.0, 1.0, 1.0, 1.0]),
            block_het_sums: pop(
                [0.6, 0.75, 1.0, 0.3, 0.3, 0.3, 0.0],
                [0.8, 1.2, 0.3, 0.3, 0.3, 0.3, 0.0],
            ),
            block_het_weights: pop(
                [2.0, 3.0, 5.0, 1.0, 0.0, 1.0, 1.0],
                [4.0, 4.0, 2.0, 1.0, 1.0, 0.0, 1.0],
            ),
            populations: vec![PopulationId::new("A"), PopulationId::new("B")],
            n_blocks: 7,
            block_lengths: vec![10; 7],
            fst: None,
            provenance: None,
        }
    }

    /// The mean of `x` and its delete-one jackknife standard error.
    fn mean_and_jackknife_se(x: &[f64]) -> (f64, f64) {
        let g = x.len() as f64;
        let mean = x.iter().sum::<f64>() / g;
        let loo = x.iter().map(|v| (mean * g - v) / (g - 1.0));
        let var = (g - 1.0) / g * loo.map(|l| (l - mean).powi(2)).sum::<f64>();
        (mean, var.sqrt())
    }

    #[test]
    fn test_fst_without_fst_blocks_is_the_mean_of_the_block_ratios() {
        let fst = f2_blocks_without_fst().fst("A", "B").unwrap();
        let (mean, se) = mean_and_jackknife_se(&[0.1 / 0.65, 0.05 / 0.65, 0.1 / 0.55]);
        assert!((fst.estimate - mean).abs() < 1e-15, "{}", fst.estimate);
        assert!((fst.std_err - se).abs() < 1e-15, "{}", fst.std_err);
        // the SNPs of the blocks used
        assert_eq!(fst.n_snps, 9);
    }

    #[test]
    fn test_fst_without_fst_blocks_of_two_blocks() {
        let fst = f2_blocks_without_fst()
            .select_blocks(&[0, 1])
            .unwrap()
            .fst("A", "B")
            .unwrap();
        let (_, se) = mean_and_jackknife_se(&[0.1 / 0.65, 0.05 / 0.65]);
        assert!((fst.std_err - se).abs() < 1e-15, "{}", fst.std_err);
    }

    #[test]
    fn test_fst_without_fst_blocks_of_one_block_has_no_se() {
        let fst = f2_blocks_without_fst()
            .select_blocks(&[0])
            .unwrap()
            .fst("A", "B")
            .unwrap();
        assert_eq!(fst.std_err, 0.0);
    }

    #[test]
    fn test_self_fst_counts_the_snps_with_a_heterozygosity() {
        let fst = f2_blocks_without_fst().fst("A", "A").unwrap();
        assert_eq!(fst.n_snps, 13);
    }

    #[test]
    fn test_fst_unknown_population() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert!(f2_blocks.fst("PopA", "Unknown").is_none());
        assert!(f2_blocks.fst("Unknown", "PopA").is_none());
    }
}
