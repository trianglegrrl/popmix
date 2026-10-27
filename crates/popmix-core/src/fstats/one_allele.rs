//! Populations with at most one allele per SNP (TASKS.md G7).
//!
//! f2(A, B) subtracts each population's sampling bias, p(1 − p)/(n − 1)
//! from its n alleles. With one allele (a single pseudo-haploid sample) the
//! bias cannot be estimated, so f2 with that population, and f3 with it as
//! the target, keep its sampling variance and come out too large; in f4,
//! and in f3 with it as a source, the term cancels. admixtools2 computes
//! them all the same, and its `qp3pop` on genotypes gives NaN with n = 0
//! for such a target.

use super::f2_blocks::F2Blocks;

impl F2Blocks {
    /// Whether `pop` has at most one allele at every SNP, so that its
    /// sampling bias was not corrected. False when the f2 blocks do not
    /// know (those read from admixtools2's f2 directories, whose
    /// heterozygosity is not stored) or `pop` is unknown.
    pub fn has_one_allele(&self, pop: &str) -> bool {
        let Some(i) = self.pop_index(pop) else {
            return false;
        };
        let sums = self.block_het_sums().column(i);
        let known = sums.iter().all(|s| s.is_finite());
        known && self.block_het_weights().column(i).sum() == 0.0
    }
}

#[cfg(test)]
mod tests {
    use crate::fstats::blocks::BlockDefinition;
    use crate::fstats::F2Blocks;
    use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
    use crate::genotypes::{Chromosome, Ploidy, PopulationId, Sample, Snp};

    fn f2_blocks() -> F2Blocks {
        let mut samples: Vec<Sample> = (0..4)
            .map(|i| Sample::new(format!("A{i}"), PopulationId::new("A")))
            .collect();
        samples.push(Sample::with_ploidy(
            "T0",
            PopulationId::new("T"),
            Ploidy::Pseudohaploid,
        ));
        let snps: Vec<Snp> = (0..30u32)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    1000 * (i + 1),
                    f64::from(i) * 0.01,
                    'A',
                    'G',
                )
            })
            .collect();
        let data = ndarray::Array2::from_shape_fn((30, 5), |(s, j)| {
            if j == 4 {
                if s % 2 == 0 {
                    0
                } else {
                    2
                }
            } else {
                ((s + j) % 3) as u8
            }
        });
        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        F2Blocks::from_allele_frequencies(
            &AlleleFrequencies::from_genotypes(&geno),
            &BlockDefinition::default(),
        )
    }

    #[test]
    fn a_single_pseudohaploid_sample_has_one_allele() {
        let f2 = f2_blocks();
        assert!(f2.has_one_allele("T"));
        assert!(!f2.has_one_allele("A"));
        assert!(!f2.has_one_allele("nobody"));
    }
}
