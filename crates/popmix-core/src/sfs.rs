//! Site-frequency spectra.
//!
//! Two constructions, from per-population allele counts (the counted allele
//! of the genotype files; pseudo-haploid samples count one allele):
//!
//! - [`joint_sfs`], as admixtools2 `joint_sfs()`: for every combination of
//!   counts from 0 to each population's largest sample size, the number of
//!   SNPs with those counts. A SNP with missing data enters with the count
//!   among its observed alleles; a SNP where some population has no data is
//!   left out.
//! - [`projected_sfs`], as dadi `Spectrum.from_data_dict`: every SNP is
//!   projected down to the given sample sizes by sampling without
//!   replacement (hypergeometric), so SNPs with some missing data still count
//!   in full; SNPs with fewer observed alleles than the projection in some
//!   population are left out. Optionally folded, as dadi.
//!
//! [`SfsSummary`] gives dadi's one-population statistics: segregating sites,
//! pi, Watterson's theta and Tajima's D.

use std::collections::HashMap;

use statrs::function::gamma::ln_gamma;

use crate::error::{PopmixError, Result};
use crate::genotype_matrix::AlleleFrequencies;

/// A spectrum over `dims` (sample size + 1 per population), row-major with
/// the first population varying slowest, and dadi's mask.
#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    pub dims: Vec<usize>,
    pub data: Vec<f64>,
    pub mask: Vec<bool>,
}

impl Spectrum {
    fn zeros(dims: Vec<usize>) -> Self {
        let len = dims.iter().product();
        Spectrum {
            dims,
            data: vec![0.0; len],
            mask: vec![false; len],
        }
    }

    /// Counts per population at each flat index.
    fn coords(&self, mut idx: usize) -> Vec<usize> {
        let mut c = vec![0; self.dims.len()];
        for (k, &d) in self.dims.iter().enumerate().rev() {
            c[k] = idx % d;
            idx /= d;
        }
        c
    }

    fn index(&self, coords: &[usize]) -> usize {
        coords
            .iter()
            .zip(&self.dims)
            .fold(0, |i, (&c, &d)| i * d + c)
    }

    /// The index of the same entry with every population's count reversed.
    fn reversed(&self, idx: usize) -> usize {
        let c: Vec<usize> = self
            .coords(idx)
            .iter()
            .zip(&self.dims)
            .map(|(&c, &d)| d - 1 - c)
            .collect();
        self.index(&c)
    }

    /// Mask the "observed in none" and "observed in all" entries.
    fn mask_corners(&mut self) {
        let last = self.data.len() - 1;
        self.mask[0] = true;
        self.mask[last] = true;
    }

    /// dadi `Spectrum.fold`: entries with more than half the alleles derived
    /// are added onto their reverse and masked; entries with exactly half are
    /// averaged with their reverse.
    pub fn fold(&self) -> Spectrum {
        let total: usize = self.dims.iter().map(|d| d - 1).sum();
        let half = (total / 2) as f64;
        let per_entry = |i: usize| self.coords(i).iter().sum::<usize>() as f64;
        let mut out = self.clone();
        for i in 0..self.data.len() {
            let r = self.reversed(i);
            let t = per_entry(i);
            let mut v = self.data[i];
            let tr = per_entry(r);
            if tr > half {
                v += self.data[r];
            }
            if t > half {
                v = 0.0;
            }
            if t == total as f64 / 2.0 {
                v += -0.5 * self.data[i] + 0.5 * self.data[r];
            }
            out.data[i] = v;
            out.mask[i] = self.mask[i] || self.mask[r] || t > half;
        }
        out
    }
}

/// Which admixtools2 `joint_sfs()` quirks to reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct JointSfsCompat {
    /// ADMIXTOOLS2-0003: the allele count is frequency times sample size in
    /// floating point, and a SNP whose product is not exactly an integer
    /// (23/42 * 42 = 23.000000000000004) matches no entry and is dropped.
    pub drop_inexact_counts: bool,
}

impl JointSfsCompat {
    pub const CORRECTED: JointSfsCompat = JointSfsCompat {
        drop_inexact_counts: false,
    };
    pub const ADMIXTOOLS2_2_0_10: JointSfsCompat = JointSfsCompat {
        drop_inexact_counts: true,
    };
}

/// What happened to the SNPs of a joint spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JointSfsSnps {
    pub n_snps: usize,
    pub in_sfs: usize,
    /// Left out because some population has no data.
    pub missing: usize,
    /// Left out because a count is not an integer: a pseudo-haploid sample
    /// with a heterozygous call (half an allele), or, in compat mode, a count
    /// that is not exactly an integer in floating point.
    pub not_integer: usize,
}

/// The allele count, if it is (within rounding) an integer.
fn integer_count(f: f64, n: u32, exact: bool) -> Option<u32> {
    let x = f * n as f64;
    let r = x.round();
    let ok = if exact { x == r } else { (x - r).abs() < 1e-6 };
    ok.then_some(r as u32)
}

/// admixtools2 `joint_sfs()` over the populations of `af`.
pub fn joint_sfs(af: &AlleleFrequencies, compat: JointSfsCompat) -> (Spectrum, JointSfsSnps) {
    let (n_snps, n_pops) = af.frequencies.dim();
    let dims: Vec<usize> = (0..n_pops)
        .map(|p| af.counts.column(p).iter().copied().max().unwrap_or(0) as usize + 1)
        .collect();
    let mut sfs = Spectrum::zeros(dims);
    let mut snps = JointSfsSnps {
        n_snps,
        in_sfs: 0,
        missing: 0,
        not_integer: 0,
    };
    let mut coords = vec![0usize; n_pops];
    'snp: for s in 0..n_snps {
        for p in 0..n_pops {
            let (f, n) = (af.frequencies[[s, p]], af.counts[[s, p]]);
            if n == 0 || !f.is_finite() {
                snps.missing += 1;
                continue 'snp;
            }
        }
        for (p, c) in coords.iter_mut().enumerate() {
            match integer_count(
                af.frequencies[[s, p]],
                af.counts[[s, p]],
                compat.drop_inexact_counts,
            ) {
                Some(k) => *c = k as usize,
                None => {
                    snps.not_integer += 1;
                    continue 'snp;
                }
            }
        }
        let i = sfs.index(&coords);
        sfs.data[i] += 1.0;
        snps.in_sfs += 1;
    }
    (sfs, snps)
}

/// ln C(n, k), minus infinity outside 0 <= k <= n (dadi `_lncomb`).
fn ln_comb(n: i64, k: i64) -> f64 {
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

/// Probabilities of 0..=to derived alleles in `to` of `from` alleles of which
/// `hits` are derived (dadi `_cached_projection`).
fn projection(to: usize, from: usize, hits: usize) -> Vec<f64> {
    if from < to {
        return vec![0.0; to + 1];
    }
    let (to, from, hits) = (to as i64, from as i64, hits as i64);
    (0..=to)
        .map(|j| (ln_comb(to, j) + ln_comb(from - to, hits - j) - ln_comb(from, hits)).exp())
        .collect()
}

/// dadi `Spectrum.from_data_dict` (corners masked) over the populations of
/// `af`, projected to `sizes`, folded if `folded`. Returns the spectrum and
/// the number of SNPs left out for a count that is not an integer (a
/// pseudo-haploid heterozygote).
pub fn projected_sfs(
    af: &AlleleFrequencies,
    sizes: &[usize],
    folded: bool,
) -> Result<(Spectrum, usize)> {
    let (n_snps, n_pops) = af.frequencies.dim();
    if sizes.len() != n_pops {
        return Err(PopmixError::DimensionMismatch {
            context: "SFS projection".to_string(),
            expected: format!("{n_pops} sample sizes"),
            actual: format!("{}", sizes.len()),
        });
    }
    // SNPs by their (observed, derived) counts, as dadi count_data_dict
    let mut configs: HashMap<Vec<(usize, usize)>, f64> = HashMap::new();
    let mut not_integer = 0;
    'snp: for s in 0..n_snps {
        let mut key = Vec::with_capacity(n_pops);
        for p in 0..n_pops {
            let n = af.counts[[s, p]];
            let k = if n == 0 {
                0
            } else {
                match integer_count(af.frequencies[[s, p]], n, false) {
                    Some(k) => k,
                    None => {
                        not_integer += 1;
                        continue 'snp;
                    }
                }
            };
            key.push((n as usize, k as usize));
        }
        *configs.entry(key).or_insert(0.0) += 1.0;
    }
    let dims: Vec<usize> = sizes.iter().map(|s| s + 1).collect();
    let mut sfs = Spectrum::zeros(dims);
    // a fixed order, so the floating-point sums do not depend on hashing
    let mut configs: Vec<(Vec<(usize, usize)>, f64)> = configs.into_iter().collect();
    configs.sort_by(|a, b| a.0.cmp(&b.0));
    for (key, count) in configs {
        let contribs: Vec<Vec<f64>> = key
            .iter()
            .zip(sizes)
            .map(|(&(n, k), &to)| projection(to, n, k))
            .collect();
        // outer product, first population slowest
        let mut prod = vec![count];
        for c in &contribs {
            prod = prod
                .iter()
                .flat_map(|&x| c.iter().map(move |&y| x * y))
                .collect();
        }
        for (d, p) in sfs.data.iter_mut().zip(&prod) {
            *d += p;
        }
    }
    sfs.mask_corners();
    Ok((if folded { sfs.fold() } else { sfs }, not_integer))
}

/// dadi's statistics of a one-population spectrum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SfsSummary {
    /// Segregating sites (the spectrum without its corners).
    pub s: f64,
    /// Expected pairwise differences.
    pub pi: f64,
    pub watterson_theta: f64,
    pub tajima_d: f64,
}

impl SfsSummary {
    /// From an unfolded one-population spectrum.
    pub fn of(sfs: &Spectrum) -> Result<SfsSummary> {
        if sfs.dims.len() != 1 {
            return Err(PopmixError::Numerical(
                "summary statistics need a one-population spectrum".to_string(),
            ));
        }
        let n = sfs.dims[0] - 1;
        let visible = |i: usize| {
            if sfs.mask[i] || i == 0 || i == n {
                0.0
            } else {
                sfs.data[i]
            }
        };
        let s: f64 = (0..=n).map(visible).sum();
        let nf = n as f64;
        let pi = (0..=n)
            .map(|i| i as f64 * (nf - i as f64) * visible(i))
            .sum::<f64>()
            * 2.0
            / (nf * (nf - 1.0));
        let a1: f64 = (1..n).map(|i| 1.0 / i as f64).sum();
        let a2: f64 = (1..n).map(|i| 1.0 / (i as f64).powi(2)).sum();
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / (a1 * a1);
        let e1 = c1 / a1;
        let e2 = c2 / (a1 * a1 + a2);
        let watterson_theta = s / a1;
        let tajima_d = (pi - watterson_theta) / (e1 * s + e2 * s * (s - 1.0)).sqrt();
        Ok(SfsSummary {
            s,
            pi,
            watterson_theta,
            tajima_d,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{Chromosome, PopulationId, Snp};
    use ndarray::Array2;

    fn af(rows: &[&[(f64, u32)]]) -> AlleleFrequencies {
        let n_pops = rows[0].len();
        let frequencies = Array2::from_shape_fn((rows.len(), n_pops), |(s, p)| rows[s][p].0);
        let counts = Array2::from_shape_fn((rows.len(), n_pops), |(s, p)| rows[s][p].1);
        AlleleFrequencies {
            frequencies,
            counts,
            populations: (0..n_pops)
                .map(|p| PopulationId::new(format!("P{p}")))
                .collect(),
            snps: (0..rows.len())
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
                .collect(),
        }
    }

    #[test]
    fn joint_sfs_counts_each_combination() {
        let a = af(&[
            &[(0.5, 2), (0.25, 4)],
            &[(0.5, 2), (0.25, 4)],
            &[(1.0, 2), (0.0, 4)],
        ]);
        let (sfs, snps) = joint_sfs(&a, JointSfsCompat::CORRECTED);
        assert_eq!(sfs.dims, vec![3, 5]);
        assert_eq!(sfs.data[sfs.index(&[1, 1])], 2.0);
        assert_eq!(sfs.data[sfs.index(&[2, 0])], 1.0);
        assert_eq!(snps.in_sfs, 3);
    }

    #[test]
    fn missing_populations_and_half_counts_are_left_out() {
        let a = af(&[
            &[(f64::NAN, 0), (0.25, 4)],
            &[(0.25, 2), (0.25, 4)],
            &[(0.5, 2), (0.5, 4)],
        ]);
        let (_, snps) = joint_sfs(&a, JointSfsCompat::CORRECTED);
        assert_eq!((snps.missing, snps.not_integer, snps.in_sfs), (1, 1, 1));
    }

    #[test]
    fn compat_drops_counts_that_are_not_exactly_integer() {
        let f = 23.0 / 42.0;
        assert_ne!(f * 42.0, 23.0);
        let a = af(&[&[(f, 42)]]);
        assert_eq!(joint_sfs(&a, JointSfsCompat::CORRECTED).1.in_sfs, 1);
        assert_eq!(
            joint_sfs(&a, JointSfsCompat::ADMIXTOOLS2_2_0_10).1.in_sfs,
            0
        );
    }

    #[test]
    fn projection_is_hypergeometric_and_sums_to_one() {
        let p = projection(4, 10, 3);
        // in logs, as dadi, so to rounding
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-13);
        // P(1 of 4 | 3 of 10) = C(3,1) C(7,3) / C(10,4) = 3 * 35 / 210
        assert!((p[1] - 0.5).abs() < 1e-13);
        assert_eq!(projection(5, 4, 2), vec![0.0; 6]);
    }

    #[test]
    fn projecting_to_the_full_size_is_the_count() {
        let a = af(&[&[(0.5, 4)], &[(0.25, 4)], &[(0.25, 2)]]);
        let (sfs, _) = projected_sfs(&a, &[4], false).unwrap();
        // the third SNP has too few alleles
        for (got, want) in sfs.data.iter().zip([0.0, 1.0, 1.0, 0.0, 0.0]) {
            assert!((got - want).abs() < 1e-13, "{got} != {want}");
        }
        assert_eq!(sfs.mask, vec![true, false, false, false, true]);
    }

    #[test]
    fn folding_adds_the_major_allele_side() {
        let mut s = Spectrum::zeros(vec![5]);
        s.data = vec![10.0, 1.0, 2.0, 3.0, 4.0];
        s.mask_corners();
        let f = s.fold();
        assert_eq!(f.data, vec![14.0, 4.0, 2.0, 0.0, 0.0]);
        assert_eq!(f.mask, vec![true, false, false, true, true]);
    }

    #[test]
    fn summary_of_a_neutral_spectrum_has_tajima_d_near_zero() {
        // E[xi_i] = theta / i under the standard neutral model
        let n = 20;
        let mut s = Spectrum::zeros(vec![n + 1]);
        for i in 1..n {
            s.data[i] = 1000.0 / i as f64;
        }
        s.mask_corners();
        let st = SfsSummary::of(&s).unwrap();
        assert!((st.pi - st.watterson_theta).abs() / st.pi < 1e-12);
        assert!(st.tajima_d.abs() < 1e-10);
    }
}
