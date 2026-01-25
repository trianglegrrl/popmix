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

/// Result of an f3 computation.
#[derive(Debug, Clone)]
pub struct F3Result {
    /// The f3 estimate.
    pub estimate: f64,
    /// Standard error from jackknife.
    pub std_err: f64,
    /// Z-score (estimate / std_err).
    pub z_score: f64,
    /// Number of blocks used.
    pub n_blocks: usize,
}

/// Result of an f4 computation.
#[derive(Debug, Clone)]
pub struct F4Result {
    /// The f4 estimate.
    pub estimate: f64,
    /// Standard error from jackknife.
    pub std_err: f64,
    /// Z-score (estimate / std_err).
    pub z_score: f64,
    /// Number of blocks used.
    pub n_blocks: usize,
}

/// Result of an FST computation.
#[derive(Debug, Clone)]
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

        let pair_ab = self.pair_index(idx_a, idx_b);
        let pair_ac = self.pair_index(idx_a, idx_c);
        let pair_bc = self.pair_index(idx_b, idx_c);

        // Compute f3 for each block
        let mut block_f3 = Vec::with_capacity(self.n_blocks);
        let mut block_weights = Vec::with_capacity(self.n_blocks);

        for block_idx in 0..self.n_blocks {
            let (sum_ab, weight_ab) = (
                self.block_sums[[block_idx, pair_ab]],
                self.block_weights[[block_idx, pair_ab]],
            );
            let (sum_ac, weight_ac) = (
                self.block_sums[[block_idx, pair_ac]],
                self.block_weights[[block_idx, pair_ac]],
            );
            let (sum_bc, weight_bc) = (
                self.block_sums[[block_idx, pair_bc]],
                self.block_weights[[block_idx, pair_bc]],
            );

            // Need all pairs to have data
            let min_weight = weight_ab.min(weight_ac).min(weight_bc);
            if min_weight > 0.0 {
                let f2_ab = sum_ab / weight_ab;
                let f2_ac = sum_ac / weight_ac;
                let f2_bc = sum_bc / weight_bc;

                let f3 = (f2_ab + f2_ac - f2_bc) / 2.0;
                block_f3.push(f3);
                block_weights.push(min_weight);
            }
        }

        if block_f3.is_empty() {
            return None;
        }

        // Compute simple mean of block f3 values (admixtools2 style)
        // This gives equal weight to each block regardless of SNP count
        let n_blocks = block_f3.len();
        let estimate: f64 = block_f3.iter().sum::<f64>() / n_blocks as f64;

        // Jackknife standard error using simple means
        let std_err = jackknife_std_err(&block_f3, estimate);
        let z_score = if std_err > 0.0 {
            estimate / std_err
        } else {
            0.0
        };

        Some(F3Result {
            estimate,
            std_err,
            z_score,
            n_blocks: block_f3.len(),
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

        let pair_ad = self.pair_index(idx_a, idx_d);
        let pair_bc = self.pair_index(idx_b, idx_c);
        let pair_ac = self.pair_index(idx_a, idx_c);
        let pair_bd = self.pair_index(idx_b, idx_d);

        // Compute f4 for each block
        let mut block_f4 = Vec::with_capacity(self.n_blocks);
        let mut block_weights = Vec::with_capacity(self.n_blocks);

        for block_idx in 0..self.n_blocks {
            let (sum_ad, weight_ad) = (
                self.block_sums[[block_idx, pair_ad]],
                self.block_weights[[block_idx, pair_ad]],
            );
            let (sum_bc, weight_bc) = (
                self.block_sums[[block_idx, pair_bc]],
                self.block_weights[[block_idx, pair_bc]],
            );
            let (sum_ac, weight_ac) = (
                self.block_sums[[block_idx, pair_ac]],
                self.block_weights[[block_idx, pair_ac]],
            );
            let (sum_bd, weight_bd) = (
                self.block_sums[[block_idx, pair_bd]],
                self.block_weights[[block_idx, pair_bd]],
            );

            // Need all pairs to have data
            let min_weight = weight_ad.min(weight_bc).min(weight_ac).min(weight_bd);
            if min_weight > 0.0 {
                let f2_ad = sum_ad / weight_ad;
                let f2_bc = sum_bc / weight_bc;
                let f2_ac = sum_ac / weight_ac;
                let f2_bd = sum_bd / weight_bd;

                let f4 = (f2_ad + f2_bc - f2_ac - f2_bd) / 2.0;
                block_f4.push(f4);
                block_weights.push(min_weight);
            }
        }

        if block_f4.is_empty() {
            return None;
        }

        // Compute simple mean of block f4 values (admixtools2 style)
        // This gives equal weight to each block regardless of SNP count
        let n_blocks = block_f4.len();
        let estimate: f64 = block_f4.iter().sum::<f64>() / n_blocks as f64;

        // Jackknife standard error using simple means
        let std_err = jackknife_std_err(&block_f4, estimate);
        let z_score = if std_err > 0.0 {
            estimate / std_err
        } else {
            0.0
        };

        Some(F4Result {
            estimate,
            std_err,
            z_score,
            n_blocks: block_f4.len(),
        })
    }

    /// Compute D-statistic (Patterson's D), which is equivalent to f4.
    ///
    /// D(A, B; C, D) = f4(A, B; C, D)
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

    /// Compute D-statistics (qpdstat) for multiple population combinations.
    ///
    /// This is equivalent to f4 statistics and provides a convenient way to test
    /// for gene flow between population pairs.
    ///
    /// D(A, B; C, D) tests if there's differential gene flow from C and D
    /// to A vs B. A significantly positive D indicates more gene flow
    /// from D to A (or from C to B).
    ///
    /// # Arguments
    /// * `pop1` - First population set (A in the formula)
    /// * `pop2` - Second population set (B in the formula)
    /// * `pop3` - Third population set (C in the formula)
    /// * `pop4` - Fourth population set (D in the formula)
    ///
    /// # Returns
    /// Vector of (pop1, pop2, pop3, pop4, F4Result) tuples
    pub fn qpdstat(
        &self,
        pop1: &[&str],
        pop2: &[&str],
        pop3: &[&str],
        pop4: &[&str],
    ) -> Vec<(String, String, String, String, F4Result)> {
        let mut results = Vec::new();

        for &p1 in pop1 {
            for &p2 in pop2 {
                for &p3 in pop3 {
                    for &p4 in pop4 {
                        // Skip if any pair is identical (would give 0)
                        if p1 == p2 || p3 == p4 {
                            continue;
                        }
                        if let Some(result) = self.f4(p1, p2, p3, p4) {
                            results.push((
                                p1.to_string(),
                                p2.to_string(),
                                p3.to_string(),
                                p4.to_string(),
                                result,
                            ));
                        }
                    }
                }
            }
        }

        results
    }

    /// Compute f3 outgroup statistics (qp3pop) for multiple population combinations.
    ///
    /// f3(Target; Source1, Source2) tests if Target is admixed between
    /// Source1 and Source2. A significantly negative f3 indicates admixture.
    ///
    /// # Arguments
    /// * `targets` - Target populations to test for admixture
    /// * `source1` - First source population set
    /// * `source2` - Second source population set
    ///
    /// # Returns
    /// Vector of (target, source1, source2, F3Result) tuples
    pub fn qp3pop(
        &self,
        targets: &[&str],
        source1: &[&str],
        source2: &[&str],
    ) -> Vec<(String, String, String, F3Result)> {
        let mut results = Vec::new();

        for &target in targets {
            for &s1 in source1 {
                for &s2 in source2 {
                    // Skip if sources are identical
                    if s1 == s2 {
                        continue;
                    }
                    // Skip if target is same as either source
                    if target == s1 || target == s2 {
                        continue;
                    }
                    if let Some(result) = self.f3(target, s1, s2) {
                        results.push((
                            target.to_string(),
                            s1.to_string(),
                            s2.to_string(),
                            result,
                        ));
                    }
                }
            }
        }

        results
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
            samples.push(Sample::new(format!("A{}", i), PopulationId::new("PopA")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("B{}", i), PopulationId::new("PopB")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("C{}", i), PopulationId::new("PopC")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("D{}", i), PopulationId::new("PopD")));
        }

        let snps: Vec<Snp> = (0..40)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.2, // 0.2 cM spacing
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
        assert!(result.estimate > 0.0, "FST should be positive for divergent populations");
        assert!(result.estimate <= 1.0, "FST should be <= 1.0");
        assert!(result.std_err >= 0.0, "Standard error should be non-negative");

        println!(
            "FST(A, B) = {:.4} ± {:.4}",
            result.estimate, result.std_err
        );
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

    #[test]
    fn test_fst_unknown_population() {
        let af = make_test_data_4pops();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert!(f2_blocks.fst("PopA", "Unknown").is_none());
        assert!(f2_blocks.fst("Unknown", "PopA").is_none());
    }
}
