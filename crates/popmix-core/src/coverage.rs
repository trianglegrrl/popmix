//! Each sample's genotype calls: how many SNPs it has a call at and how many
//! of those are heterozygous (TASKS.md H10). Ancient samples are often of
//! low coverage, with most SNPs uncalled and one allele read at the rest
//! (pseudo-haploid, coded 0 or 2); the counts tell which samples to project
//! rather than use to build PCA axes, which share too few SNPs for kinship,
//! and how much a ploidy call ("no heterozygous call") rests on.

use crate::fstats::GenotypeRows;
use crate::Result;

/// One sample's calls.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SampleCalls {
    /// SNPs with a call (genotype 0, 1 or 2).
    pub called: u64,
    /// Of them, heterozygous (1).
    pub heterozygous: u64,
}

/// The calls of every sample over `n_snps` rows of `n_samples` genotype
/// codes (3 missing).
pub fn sample_calls(
    rows: &mut dyn GenotypeRows,
    n_samples: usize,
    n_snps: usize,
) -> Result<Vec<SampleCalls>> {
    let mut out = vec![SampleCalls::default(); n_samples];
    let mut row = vec![0u8; n_samples];
    for _ in 0..n_snps {
        rows.next_row(&mut row)?;
        for (c, &g) in out.iter_mut().zip(&row) {
            if g <= 2 {
                c.called += 1;
                c.heterozygous += u64::from(g == 1);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::MatrixRows;
    use crate::{GenotypeMatrix, Sample, Snp};
    use ndarray::array;

    #[test]
    fn calls_and_heterozygous_calls_are_counted() {
        let snps: Vec<Snp> = (0..3)
            .map(|i| {
                Snp::new(
                    format!("s{i}"),
                    crate::genotypes::Chromosome::Autosome(1),
                    i,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let samples = vec![Sample::new("a", "P"), Sample::new("b", "P")];
        let m = GenotypeMatrix::new(array![[0, 3], [1, 2], [3, 1]], samples, snps).unwrap();
        let mut rows = MatrixRows::new(&m, vec![0, 1, 2], vec![0, 1]);
        let calls = sample_calls(&mut rows, 2, 3).unwrap();
        assert_eq!(
            calls[0],
            SampleCalls {
                called: 2,
                heterozygous: 1
            }
        );
        assert_eq!(
            calls[1],
            SampleCalls {
                called: 2,
                heterozygous: 1
            }
        );
    }
}
