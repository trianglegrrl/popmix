//! READ-style kinship for pseudo-haploid (ancient) samples, as READv2
//! (Alaçamlı et al. 2024).
//!
//! For each pair of samples, P0 is the proportion of SNPs called in both at
//! which their alleles differ (a heterozygote in either sample counts as half
//! a difference). It is normalized by the P0 of unrelated pairs, by default
//! the median over all pairs, so that unrelated pairs are near 1, first-degree
//! relatives near 0.75 and duplicates near 0.5; the kinship coefficient is 1
//! minus the normalized P0. Its standard error is a jackknife over 5 Mb
//! windows (of the unweighted mean of the windows' P0, as READv2), and
//! first-degree pairs are split into parent-offspring and siblings by the
//! share of 20 Mb windows whose normalized P0 is outside the first and second
//! degree range.

use rayon::prelude::*;

use super::king::BitGenotypes;
use crate::error::{PopmixError, Result};

/// Which READv2 quirks to reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReadCompat {
    /// READV2-0001: each window after the first also counts the last SNP of
    /// the window before it.
    pub windows_overlap: bool,
    /// READV2-0002: the block jackknife leaves out each window but the last.
    pub jackknife_skips_last_window: bool,
    /// READV2-0003: the reported 95% interval is that of the normalized P0,
    /// not of the kinship coefficient.
    pub ci_of_normalized_p0: bool,
    /// Windows are split where the base-pair bin changes, ignoring
    /// chromosome boundaries (READv2's rule; on human data it splits at every
    /// chromosome boundary anyway).
    pub windows_ignore_chromosomes: bool,
}

impl ReadCompat {
    pub const CORRECTED: ReadCompat = ReadCompat {
        windows_overlap: false,
        jackknife_skips_last_window: false,
        ci_of_normalized_p0: false,
        windows_ignore_chromosomes: false,
    };
    /// READv2 2.01 exactly.
    pub const V2_01: ReadCompat = ReadCompat {
        windows_overlap: true,
        jackknife_skips_last_window: true,
        ci_of_normalized_p0: true,
        windows_ignore_chromosomes: true,
    };
}

/// How P0 is normalized.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReadNorm {
    Median,
    Mean,
    Max,
    Value(f64),
}

#[derive(Debug, Clone)]
pub struct ReadOptions {
    /// Jackknife window, base pairs (READv2 `--window_size`).
    pub window_size: u64,
    /// Windows for telling parent-offspring from siblings.
    pub first_degree_window_size: u64,
    pub norm: ReadNorm,
    /// P0 as the mean over windows, its SE their standard error (READv1,
    /// READv2 `--window_est`), instead of one genome-wide P0.
    pub window_estimate: bool,
    /// Normalized-P0 thresholds between unrelated / third, third / second,
    /// second / first degree and first degree / identical.
    pub thresholds: [f64; 4],
    pub compat: ReadCompat,
}

/// READv2's default thresholds: 1 - 1/2^(k + 1.5)... as 0.953125, 0.90625,
/// 0.8125, 0.625.
pub const READ_THRESHOLDS: [f64; 4] = [0.953125, 0.90625, 0.8125, 0.625];

/// READv2 `--2pow` thresholds: 1 - 2^-4.5, 1 - 2^-3.5, 1 - 2^-2.5, 1 - 2^-1.5.
pub fn read_thresholds_2pow() -> [f64; 4] {
    [4.5f64, 3.5, 2.5, 1.5].map(|e| 1.0 - 1.0 / 2f64.powf(e))
}

impl Default for ReadOptions {
    fn default() -> Self {
        ReadOptions {
            window_size: 5_000_000,
            first_degree_window_size: 20_000_000,
            norm: ReadNorm::Median,
            window_estimate: false,
            thresholds: READ_THRESHOLDS,
            compat: ReadCompat::CORRECTED,
        }
    }
}

/// Effective SNPs (overlapping SNPs times the normalization value) needed
/// for a third-degree call and for telling first-degree types apart.
const EFFECTIVE_SNPS_THIRD: f64 = 3000.0;
const EFFECTIVE_SNPS_FIRST: f64 = 10000.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadRelationship {
    Unrelated,
    /// In the third-degree range but with too few SNPs to call it.
    UnrelatedOrThirdDegree,
    ThirdDegree,
    SecondDegree,
    FirstDegree,
    IdenticalOrTwins,
}

impl ReadRelationship {
    /// READv2's label.
    pub fn label(self) -> &'static str {
        match self {
            ReadRelationship::Unrelated => "Unrelated",
            ReadRelationship::UnrelatedOrThirdDegree => "Unrelated/Consistent with Third Degree",
            ReadRelationship::ThirdDegree => "Third Degree",
            ReadRelationship::SecondDegree => "Second Degree",
            ReadRelationship::FirstDegree => "First Degree",
            ReadRelationship::IdenticalOrTwins => "IdenticalTwins/SameIndividual",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FirstDegreeType {
    ParentOffspring,
    Siblings,
    /// READv2 "NotApplicable": most windows outside the first/second degree
    /// range.
    NotApplicable,
}

impl FirstDegreeType {
    pub fn label(self) -> &'static str {
        match self {
            FirstDegreeType::ParentOffspring => "Parent-offspring",
            FirstDegreeType::Siblings => "Siblings",
            FirstDegreeType::NotApplicable => "NotApplicable",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReadPair {
    /// The pair, `i < j` in sample order.
    pub i: usize,
    pub j: usize,
    pub relationship: ReadRelationship,
    /// Standard errors from the boundary to less related (upper) and more
    /// related (lower) classes; None where READv2 has NA.
    pub z_upper: Option<f64>,
    pub z_lower: Option<f64>,
    pub normalized_p0: f64,
    pub normalized_p0_se: f64,
    pub p0: f64,
    pub p0_se: f64,
    pub first_degree_type: Option<FirstDegreeType>,
    /// Share of first-degree windows outside the first/second degree range.
    pub first_degree_windows_outside: f64,
    pub overlap_snps: u64,
    pub effective_snps: f64,
    pub kinship: f64,
    /// 95% interval (low, high).
    pub kinship_ci: (f64, f64),
}

#[derive(Debug, Clone)]
pub struct ReadResult {
    pub norm_value: f64,
    /// Pairs with at least one SNP called in both, `(0, 1), (0, 2), ...`.
    pub pairs: Vec<ReadPair>,
}

/// Window index ranges (inclusive) over SNPs at `positions` on `chroms`.
fn windows(
    positions: &[u64],
    chroms: &[u32],
    size: u64,
    compat: ReadCompat,
) -> Vec<(usize, usize)> {
    let n = positions.len();
    if n == 0 {
        return Vec::new();
    }
    let bin = |k: usize| positions[k] / size;
    let ends: Vec<usize> = (0..n - 1)
        .filter(|&k| {
            bin(k) != bin(k + 1)
                || (!compat.windows_ignore_chromosomes && chroms[k] != chroms[k + 1])
        })
        .collect();
    let mut out = Vec::with_capacity(ends.len() + 1);
    let mut start = 0;
    for &e in &ends {
        out.push((start, e));
        start = if compat.windows_overlap { e } else { e + 1 };
    }
    out.push((start, n - 1));
    out
}

/// Called, het and mismatch counts of a pair up to (not including) each of
/// `points` (sorted SNP indices).
fn prefix_counts(g: &BitGenotypes, i: usize, j: usize, points: &[usize]) -> Vec<[u64; 3]> {
    let (hi, ai, bi) = g.masks(i);
    let (hj, aj, bj) = g.masks(j);
    let mut out = Vec::with_capacity(points.len());
    let mut acc = [0u64; 3];
    let mut word = 0;
    let count = |w: usize, mask: u64| -> [u64; 3] {
        let ci = hi[w] | ai[w] | bi[w];
        let cj = hj[w] | aj[w] | bj[w];
        let called = ci & cj & mask;
        [
            called.count_ones() as u64,
            ((hi[w] | hj[w]) & called).count_ones() as u64,
            (((ai[w] & bj[w]) | (bi[w] & aj[w])) & mask).count_ones() as u64,
        ]
    };
    for &p in points {
        while word < p / 64 {
            let c = count(word, u64::MAX);
            (0..3).for_each(|k| acc[k] += c[k]);
            word += 1;
        }
        let rem = p % 64;
        let partial = if rem == 0 {
            [0; 3]
        } else {
            count(word, (1u64 << rem) - 1)
        };
        out.push([
            acc[0] + partial[0],
            acc[1] + partial[1],
            acc[2] + partial[2],
        ]);
    }
    out
}

/// (called, P0) over each window, None for windows with no SNP called.
fn window_p0(
    prefix: &[[u64; 3]],
    at: &dyn Fn(usize) -> usize,
    wins: &[(usize, usize)],
) -> Vec<(u64, f64)> {
    wins.iter()
        .filter_map(|&(a, b)| {
            let (lo, hi) = (prefix[at(a)], prefix[at(b + 1)]);
            let called = hi[0] - lo[0];
            if called == 0 {
                return None;
            }
            let het = (hi[1] - lo[1]) as f64 / 2.0;
            let mismatch = (hi[2] - lo[2]) as f64;
            Some((called, (mismatch + het) / called as f64))
        })
        .collect()
}

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

/// READv2 `blockJN_SE`: jackknife of the mean of `x`, leaving out each
/// element (each but the last with `skip_last`).
fn jackknife_se(x: &[f64], skip_last: bool) -> f64 {
    let n = x.len();
    let total: f64 = x.iter().sum();
    let reps: Vec<f64> = (0..if skip_last { n.saturating_sub(1) } else { n })
        .map(|k| (total - x[k]) / (n - 1) as f64)
        .collect();
    let grand = mean(&reps);
    let var = reps.iter().map(|m| (m - grand).powi(2)).sum::<f64>() * ((n - 1) as f64 / n as f64);
    var.sqrt()
}

fn median(mut x: Vec<f64>) -> f64 {
    x.sort_by(f64::total_cmp);
    let n = x.len();
    if n % 2 == 1 {
        x[n / 2]
    } else {
        (x[n / 2 - 1] + x[n / 2]) / 2.0
    }
}

/// READ-style kinship for every pair of samples.
///
/// `positions` and `chroms` give each SNP's base-pair position and
/// chromosome, in the order of `g`.
pub fn read_kinship(
    g: &BitGenotypes,
    positions: &[u64],
    chroms: &[u32],
    options: &ReadOptions,
) -> Result<ReadResult> {
    if positions.len() != g.n_snps() || chroms.len() != g.n_snps() {
        return Err(PopmixError::DimensionMismatch {
            context: "READ SNP positions".to_string(),
            expected: format!("{} SNPs", g.n_snps()),
            actual: format!(
                "{} positions, {} chromosomes",
                positions.len(),
                chroms.len()
            ),
        });
    }
    let n = g.n_samples();
    if n < 2 {
        return Err(PopmixError::Numerical(
            "READ needs at least two samples".to_string(),
        ));
    }
    let compat = options.compat;
    let jk_wins = windows(positions, chroms, options.window_size, compat);
    let fd_wins = windows(positions, chroms, options.first_degree_window_size, compat);
    let mut points: Vec<usize> = jk_wins
        .iter()
        .chain(&fd_wins)
        .flat_map(|&(a, b)| [a, b + 1])
        .chain([0, g.n_snps()])
        .collect();
    points.sort_unstable();
    points.dedup();
    let at = |p: usize| points.binary_search(&p).expect("window boundary");

    struct Raw {
        i: usize,
        j: usize,
        called: u64,
        p0: f64,
        jk: Vec<f64>,
        fd: Vec<f64>,
    }
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();
    let raw: Vec<Raw> = pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            let prefix = prefix_counts(g, i, j, &points);
            let jk = window_p0(&prefix, &at, &jk_wins);
            let fd: Vec<f64> = window_p0(&prefix, &at, &fd_wins)
                .into_iter()
                .map(|w| w.1)
                .collect();
            let (called, p0) = if options.window_estimate {
                if jk.is_empty() {
                    return None;
                }
                (
                    jk.iter().map(|w| w.0).sum(),
                    mean(&jk.iter().map(|w| w.1).collect::<Vec<_>>()),
                )
            } else {
                let all = window_p0(&prefix, &at, &[(0, g.n_snps() - 1)]);
                let &(called, p0) = all.first()?;
                (called, p0)
            };
            Some(Raw {
                i,
                j,
                called,
                p0,
                jk: jk.into_iter().map(|w| w.1).collect(),
                fd,
            })
        })
        .collect();
    if raw.is_empty() {
        return Err(PopmixError::Numerical(
            "no pair of samples has a SNP called in both".to_string(),
        ));
    }

    let p0s: Vec<f64> = raw.iter().map(|r| r.p0).collect();
    let norm_value = match options.norm {
        ReadNorm::Median => median(p0s),
        ReadNorm::Mean => mean(&p0s),
        ReadNorm::Max => p0s.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        ReadNorm::Value(v) => v,
    };
    if norm_value <= 0.0 || !norm_value.is_finite() {
        return Err(PopmixError::Numerical(format!(
            "READ normalization value is {norm_value} (are most samples duplicates?)"
        )));
    }

    let th = options.thresholds;
    let out = raw
        .into_iter()
        .map(|r| {
            let x = r.p0 / norm_value;
            let jk: Vec<f64> = r.jk.iter().map(|p| p / norm_value).collect();
            let se = if options.window_estimate {
                let m = mean(&jk);
                let sd = (jk.iter().map(|v| (v - m).powi(2)).sum::<f64>()
                    / (jk.len() as f64 - 1.0))
                    .sqrt();
                sd / (jk.len() as f64).sqrt()
            } else {
                jackknife_se(&jk, compat.jackknife_skips_last_window)
            };
            let outside = r.fd.iter().filter(|p| {
                let v = *p / norm_value;
                v > th[1] || v < th[3]
            });
            let first_degree_windows_outside = outside.count() as f64 / r.fd.len() as f64;
            let effective = r.called as f64 * norm_value;
            let z = |t: f64| ((t - x) / se).abs();
            let (relationship, z_upper, z_lower) = if x <= th[0] && x >= th[1] {
                let rel = if effective < EFFECTIVE_SNPS_THIRD {
                    ReadRelationship::UnrelatedOrThirdDegree
                } else {
                    ReadRelationship::ThirdDegree
                };
                (rel, Some(z(th[0])), Some(z(th[1])))
            } else if x <= th[1] && x >= th[2] {
                (
                    ReadRelationship::SecondDegree,
                    Some(z(th[1])),
                    Some(z(th[2])),
                )
            } else if x <= th[2] && x >= th[3] {
                (
                    ReadRelationship::FirstDegree,
                    Some(z(th[2])),
                    Some(z(th[3])),
                )
            } else if x <= th[3] {
                (ReadRelationship::IdenticalOrTwins, Some(z(th[3])), None)
            } else {
                // READv2 measures unrelated pairs against the second-degree
                // boundary
                (ReadRelationship::Unrelated, None, Some(z(th[1])))
            };
            let first_degree_type = if relationship == ReadRelationship::FirstDegree
                && effective >= EFFECTIVE_SNPS_FIRST
            {
                if first_degree_windows_outside >= 0.6 {
                    Some(FirstDegreeType::NotApplicable)
                } else if first_degree_windows_outside >= 0.35 {
                    Some(FirstDegreeType::Siblings)
                } else if first_degree_windows_outside <= 0.3 {
                    Some(FirstDegreeType::ParentOffspring)
                } else {
                    None
                }
            } else {
                None
            };
            let ci = (x - 1.96 * se, x + 1.96 * se);
            ReadPair {
                i: r.i,
                j: r.j,
                relationship,
                z_upper,
                z_lower,
                normalized_p0: x,
                normalized_p0_se: se,
                p0: r.p0,
                p0_se: se * norm_value,
                first_degree_type,
                first_degree_windows_outside,
                overlap_snps: r.called,
                effective_snps: effective,
                kinship: 1.0 - x,
                kinship_ci: if compat.ci_of_normalized_p0 {
                    ci
                } else {
                    (1.0 - ci.1, 1.0 - ci.0)
                },
            }
        })
        .collect();
    Ok(ReadResult {
        norm_value,
        pairs: out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::MatrixRows;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    fn bits(geno: &Array2<u8>) -> BitGenotypes {
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
        let m = GenotypeMatrix::new(geno.clone(), samples, snps).unwrap();
        let mut rows = MatrixRows::new(&m, (0..n_snps).collect(), (0..n).collect());
        BitGenotypes::read(&mut rows, n, n_snps).unwrap()
    }

    #[test]
    fn overlapping_windows_share_their_boundary_snp() {
        let pos = [1, 2, 6, 7, 11];
        let chr = [1; 5];
        assert_eq!(
            windows(&pos, &chr, 5, ReadCompat::CORRECTED),
            vec![(0, 1), (2, 3), (4, 4)]
        );
        assert_eq!(
            windows(&pos, &chr, 5, ReadCompat::V2_01),
            vec![(0, 1), (1, 3), (3, 4)]
        );
    }

    #[test]
    fn chromosome_boundaries_split_windows() {
        let pos = [1, 2, 3, 4];
        let chr = [1, 1, 2, 2];
        assert_eq!(
            windows(&pos, &chr, 100, ReadCompat::CORRECTED),
            vec![(0, 1), (2, 3)]
        );
        assert_eq!(windows(&pos, &chr, 100, ReadCompat::V2_01), vec![(0, 3)]);
    }

    #[test]
    fn jackknife_of_equal_windows_is_zero_and_skipping_one_changes_it() {
        assert_eq!(jackknife_se(&[0.2, 0.2, 0.2], false), 0.0);
        let x = [0.1, 0.2, 0.4, 0.3];
        // leave-one-out means 0.3, 0.2667, 0.2, 0.2333: the usual jackknife
        let full = jackknife_se(&x, false);
        let sd = {
            let m = mean(&x);
            (x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / 3.0).sqrt()
        };
        assert!((full - sd / 2.0).abs() < 1e-15);
        assert!((jackknife_se(&x, true) - full).abs() > 1e-3);
    }

    #[test]
    fn identical_samples_are_identical_and_halves_are_counted() {
        // pseudo-haploid samples: 0 and 1 identical, 2 differs at half the SNPs,
        // 3 has a heterozygote
        let n = 4000;
        let geno = Array2::from_shape_fn((n, 4), |(s, k)| match k {
            0 | 1 => 2 * ((s / 3) % 2) as u8,
            2 => 2 * ((s / 7) % 2) as u8,
            _ => {
                if s == 0 {
                    1
                } else {
                    2 * ((s / 5) % 2) as u8
                }
            }
        });
        let pos: Vec<u64> = (0..n as u64).map(|s| s * 1000).collect();
        let opts = ReadOptions {
            norm: ReadNorm::Value(0.5),
            ..Default::default()
        };
        let r = read_kinship(&bits(&geno), &pos, &vec![1; n], &opts).unwrap();
        let p = |i, j| r.pairs.iter().find(|p| p.i == i && p.j == j).unwrap();
        assert_eq!(p(0, 1).p0, 0.0);
        assert_eq!(p(0, 1).relationship, ReadRelationship::IdenticalOrTwins);
        assert_eq!(p(0, 1).kinship, 1.0);
        // the heterozygote at SNP 0 counts half a mismatch
        let mism = (1..n).filter(|&s| geno[[s, 0]] != geno[[s, 3]]).count() as f64 + 0.5;
        assert!((p(0, 3).p0 - mism / n as f64).abs() < 1e-15);
    }

    #[test]
    fn kinship_interval_is_around_the_kinship_unless_compat() {
        let n = 3000;
        let geno = Array2::from_shape_fn((n, 3), |(s, k)| 2 * (((s * (k + 3)) / 11) % 2) as u8);
        let pos: Vec<u64> = (0..n as u64).map(|s| s * 10_000).collect();
        let corrected =
            read_kinship(&bits(&geno), &pos, &vec![1; n], &ReadOptions::default()).unwrap();
        let p = &corrected.pairs[0];
        assert!(p.kinship_ci.0 < p.kinship && p.kinship < p.kinship_ci.1);
        let compat = ReadOptions {
            compat: ReadCompat::V2_01,
            ..Default::default()
        };
        let v = read_kinship(&bits(&geno), &pos, &vec![1; n], &compat).unwrap();
        let q = &v.pairs[0];
        assert!(q.kinship_ci.0 < q.normalized_p0 && q.normalized_p0 < q.kinship_ci.1);
    }
}
