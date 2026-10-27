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
        if n_samples == 0 {
            return Ok(g);
        }
        // a batch of whole words of SNPs at a time (rows decoded in
        // parallel), then each word of every sample from its 64 rows, read
        // in order, the words in parallel
        const BATCH_WORDS: usize = 256;
        let mut batch = vec![0u8; BATCH_WORDS * 64 * n_samples];
        for w0 in (0..words).step_by(BATCH_WORDS) {
            let first = w0 * 64;
            let n_rows = (first + BATCH_WORDS * 64).min(n_snps) - first;
            rows.next_rows(&mut batch[..n_rows * n_samples], n_samples)?;
            let batch = &batch[..n_rows * n_samples];
            let span = n_rows.div_ceil(64);
            // per word: the het, hom0 and hom2 bits of every sample
            let built: Vec<Vec<[u64; 3]>> = (0..span)
                .into_par_iter()
                .map(|w| {
                    let mut bits = vec![[0u64; 3]; n_samples];
                    let rows = w * 64..((w + 1) * 64).min(n_rows);
                    for (bit, r) in rows.enumerate() {
                        let row = &batch[r * n_samples..(r + 1) * n_samples];
                        for (b, &code) in bits.iter_mut().zip(row) {
                            b[0] |= ((code == 1) as u64) << bit;
                            b[1] |= ((code == 0) as u64) << bit;
                            b[2] |= ((code == 2) as u64) << bit;
                        }
                    }
                    bits
                })
                .collect();
            for (w, bits) in built.iter().enumerate() {
                for (s, b) in bits.iter().enumerate() {
                    let at = s * words + w0 + w;
                    g.het[at] = b[0];
                    g.hom0[at] = b[1];
                    g.hom2[at] = b[2];
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

/// Samples per side of a tile of pairs, and 64-SNP words per chunk: a
/// tile's bit sets for one chunk stay in cache while all its pairs are
/// counted.
const TILE: usize = 32;
const CHUNK: usize = 512;

/// KING-robust kinship for every pair, in PLINK 2's order.
///
/// Pairs are counted a tile at a time, the tiles in parallel; counts are
/// integers, so the result does not depend on the thread count.
pub fn king_robust(g: &BitGenotypes) -> Vec<KingPair> {
    let n = g.n_samples;
    let n_tiles = n.div_ceil(TILE);
    let tiles: Vec<(usize, usize)> = (0..n_tiles)
        .flat_map(|a| (0..=a).map(move |b| (a, b)))
        .collect();
    let counted: Vec<_> = tiles
        .into_par_iter()
        .map(|(a, b)| ((a, b), count_tile(g, a * TILE, b * TILE)))
        .collect();
    let mut pairs = vec![
        KingPair {
            i: 0,
            j: 0,
            nsnp: 0,
            hethet: 0,
            ibs0: 0,
            het1_hom2: 0,
            het2_hom1: 0,
            kinship: 0.0,
        };
        n * n.saturating_sub(1) / 2
    ];
    for ((a, b), counts) in counted {
        for (k, c) in counts.into_iter().enumerate() {
            let (i, j) = (a * TILE + k / TILE, b * TILE + k % TILE);
            if j < i && i < n {
                pairs[i * (i - 1) / 2 + j] = king_pair(i, j, c);
            }
        }
    }
    pairs
}

/// Counts for the pairs of samples `i0..i0 + TILE` with `j0..j0 + TILE`
/// (row-major, `TILE` per row; pairs with `j >= i` or past the last sample
/// are left zero), with AVX2 and POPCNT where the CPU has them.
fn count_tile(g: &BitGenotypes, i0: usize, j0: usize) -> Vec<[u64; 5]> {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("popcnt")
    {
        // SAFETY: the CPU has the features the function is compiled for.
        return unsafe { count_tile_avx2(g, i0, j0) };
    }
    count_tile_portable(g, i0, j0)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
unsafe fn count_tile_avx2(g: &BitGenotypes, i0: usize, j0: usize) -> Vec<[u64; 5]> {
    count_tile_portable(g, i0, j0)
}

#[inline(always)]
fn count_tile_portable(g: &BitGenotypes, i0: usize, j0: usize) -> Vec<[u64; 5]> {
    let n = g.n_samples;
    let mut acc = vec![[0u64; 5]; TILE * TILE];
    for lo in (0..g.words).step_by(CHUNK) {
        let hi = (lo + CHUNK).min(g.words);
        let part = |s: usize| {
            let r = s * g.words + lo..s * g.words + hi;
            (&g.het[r.clone()], &g.hom0[r.clone()], &g.hom2[r])
        };
        for i in i0..(i0 + TILE).min(n) {
            let a = part(i);
            for j in j0..(j0 + TILE).min(n).min(i) {
                let b = part(j);
                let c = counts(a, b);
                let e = &mut acc[(i - i0) * TILE + (j - j0)];
                for (x, y) in e.iter_mut().zip(c) {
                    *x += y;
                }
            }
        }
    }
    acc
}

type Masks<'a> = (&'a [u64], &'a [u64], &'a [u64]);

/// (SNPs called in both, het-het, IBS0, het in i and hom in j, the reverse)
/// over the words of two samples' (het, hom0, hom2) bit sets.
#[inline(always)]
fn counts((het_i, a0, a2): Masks, (het_j, b0, b2): Masks) -> [u64; 5] {
    // equal lengths up front, so the loop has no bounds checks and vectorizes
    let n = het_i.len();
    let (a0, a2, het_j, b0, b2) = (&a0[..n], &a2[..n], &het_j[..n], &b0[..n], &b2[..n]);
    let (mut nsnp, mut hethet, mut ibs0, mut h1, mut h2) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for w in 0..n {
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
    [nsnp, hethet, ibs0, h1, h2]
}

fn king_pair(i: usize, j: usize, [nsnp, hethet, ibs0, h1, h2]: [u64; 5]) -> KingPair {
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
    fn counts_match_a_direct_count_over_many_samples_and_snps() {
        // several tiles of samples and chunks of SNPs, with missing calls
        let (n_snps, n) = (5_000, 75);
        let mut state = 99u64;
        let geno = Array2::from_shape_fn((n_snps, n), |_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % 4) as u8
        });
        let pairs = king_robust(&bits(geno.clone()));
        assert_eq!(pairs.len(), n * (n - 1) / 2);
        let mut k = 0;
        for i in 1..n {
            for j in 0..i {
                let p = pairs[k];
                k += 1;
                assert_eq!((p.i, p.j), (i, j));
                let (mut nsnp, mut hethet, mut ibs0, mut h1, mut h2) = (0, 0, 0, 0, 0);
                for snp in 0..n_snps {
                    let (a, b) = (geno[[snp, i]], geno[[snp, j]]);
                    if a > 2 || b > 2 {
                        continue;
                    }
                    nsnp += 1;
                    hethet += (a == 1 && b == 1) as u64;
                    ibs0 += (a.abs_diff(b) == 2) as u64;
                    h1 += (a == 1 && b != 1) as u64;
                    h2 += (b == 1 && a != 1) as u64;
                }
                assert_eq!(
                    (p.nsnp, p.hethet, p.ibs0, p.het1_hom2, p.het2_hom1),
                    (nsnp, hethet, ibs0, h1, h2),
                    "pair ({i}, {j})"
                );
            }
        }
    }

    #[test]
    fn the_fast_kernel_counts_as_the_portable_one() {
        let mut state = 7u64;
        let geno = Array2::from_shape_fn((3_000, 40), |_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % 4) as u8
        });
        let g = bits(geno);
        assert_eq!(count_tile(&g, 32, 0), count_tile_portable(&g, 32, 0));
        assert_eq!(count_tile(&g, 32, 32), count_tile_portable(&g, 32, 32));
    }

    #[test]
    fn pairs_are_in_plink_order() {
        let g = Array2::<u8>::zeros((10, 4));
        let order: Vec<(usize, usize)> = king_robust(&bits(g)).iter().map(|p| (p.i, p.j)).collect();
        assert_eq!(order, vec![(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]);
    }
}
