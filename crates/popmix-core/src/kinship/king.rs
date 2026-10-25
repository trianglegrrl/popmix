//! KING-robust kinship (Manichaikul et al. 2010), as PLINK 2
//! `--make-king`/`--make-king-table`.
//!
//! For two samples, over the SNPs where both are called: the het-het count,
//! the opposite-homozygote (IBS0) count, and the het-hom counts each way.
//! The between-family estimator is
//!
//! ```text
//! kinship = 1/2 - (4 IBS0 + HET1_HOM2 + HET2_HOM1) / (4 min(het1, het2))
//! ```
//!
//! with each sample's heterozygote count taken over the same SNPs; this is
//! `(HETHET - 2 IBS0) / (2 min) + 1/2 - (het1 + het2) / (4 min)` rearranged.
//! A sample with no heterozygous call (as a pseudo-haploid one) gives
//! minus infinity, as in PLINK 2.

use rayon::prelude::*;

use crate::error::Result;
use crate::fstats::GenotypeRows;

/// Genotypes as bit sets over the SNPs, one per sample: heterozygous,
/// homozygous for allele 0 (genotype 0) and for allele 1 (genotype 2).
#[derive(Debug, Clone)]
pub struct BitGenotypes {
    n_samples: usize,
    n_snps: usize,
    words: usize,
    het: Vec<u64>,
    hom0: Vec<u64>,
    hom2: Vec<u64>,
}

impl BitGenotypes {
    /// Read `n_snps` rows of `n_samples` genotypes (0, 1, 2; 3 missing).
    pub fn read(rows: &mut dyn GenotypeRows, n_samples: usize, n_snps: usize) -> Result<Self> {
        let words = n_snps.div_ceil(64);
        let mut g = BitGenotypes {
            n_samples,
            n_snps,
            words,
            het: vec![0; n_samples * words],
            hom0: vec![0; n_samples * words],
            hom2: vec![0; n_samples * words],
        };
        let mut row = vec![0u8; n_samples];
        for snp in 0..n_snps {
            rows.next_row(&mut row)?;
            let (w, bit) = (snp / 64, 1u64 << (snp % 64));
            for (s, &x) in row.iter().enumerate() {
                let i = s * words + w;
                match x {
                    0 => g.hom0[i] |= bit,
                    1 => g.het[i] |= bit,
                    2 => g.hom2[i] |= bit,
                    _ => {}
                }
            }
        }
        Ok(g)
    }

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    pub fn n_snps(&self) -> usize {
        self.n_snps
    }

    /// One sample's (heterozygous, homozygous 0, homozygous 2) bit sets.
    pub fn masks(&self, s: usize) -> (&[u64], &[u64], &[u64]) {
        (
            self.sample(&self.het, s),
            self.sample(&self.hom0, s),
            self.sample(&self.hom2, s),
        )
    }

    fn sample<'a>(&self, v: &'a [u64], s: usize) -> &'a [u64] {
        &v[s * self.words..(s + 1) * self.words]
    }
}

/// KING-robust counts and kinship for one pair of samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KingPair {
    /// The pair, `i > j` (PLINK 2's row order: (1, 0), (2, 0), (2, 1), ...).
    pub i: usize,
    pub j: usize,
    /// SNPs with both called.
    pub nsnp: u64,
    pub hethet: u64,
    pub ibs0: u64,
    /// Het in `i` and homozygous in `j`, and the reverse.
    pub het1_hom2: u64,
    pub het2_hom1: u64,
    pub kinship: f64,
}

/// KING-robust kinship for every pair, in PLINK 2's order.
pub fn king_robust(g: &BitGenotypes) -> Vec<KingPair> {
    let n = g.n_samples;
    (1..n)
        .into_par_iter()
        .flat_map_iter(|i| (0..i).map(move |j| king_pair(g, i, j)))
        .collect()
}

fn king_pair(g: &BitGenotypes, i: usize, j: usize) -> KingPair {
    let (het_i, a0, a2) = (
        g.sample(&g.het, i),
        g.sample(&g.hom0, i),
        g.sample(&g.hom2, i),
    );
    let (het_j, b0, b2) = (
        g.sample(&g.het, j),
        g.sample(&g.hom0, j),
        g.sample(&g.hom2, j),
    );
    let (mut nsnp, mut hethet, mut ibs0, mut h1, mut h2) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for w in 0..g.words {
        let called_i = het_i[w] | a0[w] | a2[w];
        let called_j = het_j[w] | b0[w] | b2[w];
        let hom_i = a0[w] | a2[w];
        let hom_j = b0[w] | b2[w];
        nsnp += (called_i & called_j).count_ones() as u64;
        hethet += (het_i[w] & het_j[w]).count_ones() as u64;
        ibs0 += ((a0[w] & b2[w]) | (a2[w] & b0[w])).count_ones() as u64;
        h1 += (het_i[w] & hom_j).count_ones() as u64;
        h2 += (het_j[w] & hom_i).count_ones() as u64;
    }
    let min_het = (hethet + h1).min(hethet + h2);
    let kinship = 0.5 - (4 * ibs0 + h1 + h2) as f64 / (4 * min_het) as f64;
    KingPair {
        i,
        j,
        nsnp,
        hethet,
        ibs0,
        het1_hom2: h1,
        het2_hom1: h2,
        kinship,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::MatrixRows;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    fn bits(geno: Array2<u8>) -> BitGenotypes {
        let (n_snps, n) = geno.dim();
        let samples: Vec<Sample> = (0..n).map(|i| Sample::new(format!("s{i}"), "P")).collect();
        let snps: Vec<Snp> = (0..n_snps)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    i as u32,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let m = GenotypeMatrix::new(geno, samples, snps).unwrap();
        let mut rows = MatrixRows::new(&m, (0..n_snps).collect(), (0..n).collect());
        BitGenotypes::read(&mut rows, n, n_snps).unwrap()
    }

    /// Unrelated parents at frequency 1/2, and their child.
    fn trio(n_snps: usize) -> Array2<u8> {
        let mut state = 12345u64;
        let mut coin = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state & 1) as u8
        };
        let mut g = Array2::<u8>::zeros((n_snps, 3));
        for snp in 0..n_snps {
            let (m1, m2, f1, f2) = (coin(), coin(), coin(), coin());
            g[[snp, 0]] = m1 + m2;
            g[[snp, 1]] = f1 + f2;
            let from_m = if coin() == 1 { m1 } else { m2 };
            let from_f = if coin() == 1 { f1 } else { f2 };
            g[[snp, 2]] = from_m + from_f;
        }
        g
    }

    #[test]
    fn a_child_has_kinship_one_quarter_with_each_parent_and_parents_zero() {
        let pairs = king_robust(&bits(trio(20000)));
        let k = |i: usize, j: usize| pairs.iter().find(|p| p.i == i && p.j == j).unwrap().kinship;
        assert!((k(2, 0) - 0.25).abs() < 0.02, "{}", k(2, 0));
        assert!((k(2, 1) - 0.25).abs() < 0.02, "{}", k(2, 1));
        assert!(k(1, 0).abs() < 0.02, "{}", k(1, 0));
        // parent and child never share no allele
        assert_eq!(pairs.iter().find(|p| p.i == 2 && p.j == 0).unwrap().ibs0, 0);
    }

    #[test]
    fn a_sample_with_itself_has_kinship_one_half() {
        let t = trio(500);
        let mut g = Array2::<u8>::zeros((500, 2));
        g.column_mut(0).assign(&t.column(0));
        g.column_mut(1).assign(&t.column(0));
        let p = king_robust(&bits(g))[0];
        assert_eq!(p.kinship, 0.5);
        assert_eq!(p.ibs0, 0);
    }

    #[test]
    fn counts_use_only_snps_called_in_both() {
        // sample 1 is missing at SNP 0
        let g = ndarray::arr2(&[[1u8, 3], [1, 1], [0, 2], [2, 1], [1, 0]]);
        let p = king_robust(&bits(g))[0];
        assert_eq!((p.i, p.j), (1, 0));
        assert_eq!(p.nsnp, 4);
        assert_eq!(p.hethet, 1);
        assert_eq!(p.ibs0, 1);
        // het in sample 1 (i), hom in 0 (j): SNP 3; the reverse: SNP 4
        assert_eq!((p.het1_hom2, p.het2_hom1), (1, 1));
        assert_eq!(p.kinship, 0.5 - 6.0 / 8.0);
    }

    #[test]
    fn no_heterozygote_gives_minus_infinity() {
        let g = ndarray::arr2(&[[0u8, 2], [2, 2], [0, 0], [2, 1]]);
        assert_eq!(king_robust(&bits(g))[0].kinship, f64::NEG_INFINITY);
    }

    #[test]
    fn pairs_are_in_plink_order() {
        let g = Array2::<u8>::zeros((10, 4));
        let order: Vec<(usize, usize)> = king_robust(&bits(g)).iter().map(|p| (p.i, p.j)).collect();
        assert_eq!(order, vec![(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]);
    }
}
