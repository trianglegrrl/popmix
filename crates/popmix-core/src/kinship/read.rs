//! READ-style kinship for pseudo-haploid (ancient) genotypes.
//!
//! For every pair of samples: P0, the mean per-SNP mismatch (½ at a
//! heterozygote) over the SNPs called in both, with a block-jackknife
//! standard error over windows of `window_size` base pairs; P0 normalized
//! by one value for the run (the median P0 of the pairs by default); the
//! degree of relationship from four thresholds on the normalized P0, with
//! z scores to the class's boundaries; for first-degree pairs, parent-
//! offspring or siblings from the share of 20 Mb windows outside the
//! first/second-degree range; and the kinship θ = 1 − normalized P0 with a
//! 95 % interval.
//!
//! Provenance:
//!
//! - Monroy Kuhn, Jakobsson and Günther 2018, "Estimating genetic kin
//!   relationships in prehistoric populations", PLOS ONE 13:e0195491
//!   (READ: P0 per window, normalisation, thresholds at the mid-points of the
//!   expected normalized P0 of each degree).
//! - Alaçamlı et al. 2024, "READv2: advanced and user-friendly detection of
//!   biological relatedness in archaeogenomics", Genome Biology 25:216,
//!   Methods (genome-wide P0 with a 5 Mb block jackknife, the third-degree
//!   class with at least 3000 expected mismatches, 20 Mb windows and the
//!   0.3 / 0.35 / 0.6 proportions with at least 10 000, θ = 1 − normalized
//!   P0).
//! - The clean-room specification `docs/cleanroom/read.md`, itself from the
//!   papers, READv2's README and `--help`, and black-box runs of READv2
//!   (commit 73739bef) on hand-made datasets and the goldens under
//!   `tests/golden/readv2/`.
//! - The outputs of the specification's black-box runs of READv2 (its
//!   section 12, e1–e17) as fixtures in `read/tests/`.
//!
//! Written clean-room in C9 (TASKS.md) from the specification and the tests
//! alone, without reading READv2's source (GPL-3), an earlier version of
//! this module, or DEVIATIONS.md's READv2 entries.
//!
//! `ReadCompat::V2_01` reproduces READv2 2.01's output; the default,
//! `ReadCompat::CORRECTED`, differs only in the standard errors (windows that
//! do not overlap, respect chromosomes, and are all left out in turn) and in
//! centring the kinship interval on the kinship.

mod counts;

use rayon::prelude::*;

use crate::error::{PopmixError, Result};
use crate::kinship::king::BitGenotypes;
use counts::{segments, Count, PairBits};

/// Which of READv2 2.01's behaviours to reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReadCompat {
    /// Every window but the first also takes the last SNP of the window
    /// before it (so a SNP that ends a window is in two).
    pub windows_overlap: bool,
    /// The genome-wide jackknife leaves out every window but the last in
    /// turn (the factor still uses the number of all windows).
    pub jackknife_skips_last_window: bool,
    /// The 95 % interval is centred on the normalized P0 instead of on the
    /// kinship coefficient.
    pub ci_of_normalized_p0: bool,
    /// Windows split only where the position bin changes, not where the
    /// chromosome does.
    pub windows_ignore_chromosomes: bool,
}

impl ReadCompat {
    /// popmix's corrected form (the default).
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

/// The value P0 is normalized by.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReadNorm {
    /// The median P0 of the pairs (READv2's default).
    Median,
    /// The mean P0 of the pairs.
    Mean,
    /// The largest P0 of the pairs.
    Max,
    /// A given value.
    Value(f64),
}

/// Options of [`read_kinship`]; the defaults are READv2's. Not exhaustive:
/// build it from `ReadOptions::default()`.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ReadOptions {
    /// Base pairs per jackknife window (P0 standard error).
    pub window_size: u64,
    /// Base pairs per window of the parent-offspring / sibling call.
    pub first_degree_window_size: u64,
    /// What P0 is normalized by.
    pub norm: ReadNorm,
    /// P0 as the mean of the window P0s (READv1) instead of over all SNPs at
    /// once; its SE is then the standard error of that mean.
    pub window_estimate: bool,
    /// Normalized-P0 boundaries unrelated / third, third / second,
    /// second / first and first / identical, strictly decreasing.
    pub thresholds: [f64; 4],
    /// Which of READv2 2.01's behaviours to reproduce.
    pub compat: ReadCompat,
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

/// READ's thresholds: mid-points between the expected normalized P0 of
/// unrelated (1), third (0.9375), second (0.875), first degree (0.75) and
/// identical (0.5).
pub const READ_THRESHOLDS: [f64; 4] = [0.953125, 0.90625, 0.8125, 0.625];

/// READv2's `--2pow` thresholds, 1 − 2^−4.5 … 1 − 2^−1.5: the geometric
/// mid-points of the kinship values 1/32 … 1/2.
pub fn read_thresholds_2pow() -> [f64; 4] {
    [-4.5, -3.5, -2.5, -1.5].map(|e: f64| 1.0 - 2f64.powf(e))
}

/// Degree of relationship from the normalized P0.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadRelationship {
    Unrelated,
    /// In the third-degree range, with fewer than 3000 expected mismatches.
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

/// Type of a first-degree relationship.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FirstDegreeType {
    ParentOffspring,
    Siblings,
    /// Too many windows outside the first/second-degree range for either.
    NotApplicable,
}

impl FirstDegreeType {
    /// READv2's label.
    pub fn label(self) -> &'static str {
        match self {
            FirstDegreeType::ParentOffspring => "Parent-offspring",
            FirstDegreeType::Siblings => "Siblings",
            FirstDegreeType::NotApplicable => "NotApplicable",
        }
    }
}

/// READ's estimates for one pair of samples.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadPair {
    /// The pair, `i < j`.
    pub i: usize,
    pub j: usize,
    pub relationship: ReadRelationship,
    /// Standard errors from the normalized P0 to the class's upper (less
    /// related) boundary; `None` for unrelated pairs.
    pub z_upper: Option<f64>,
    /// Standard errors to the lower (more related) boundary, t1 for
    /// unrelated pairs; `None` for identical ones.
    pub z_lower: Option<f64>,
    /// P0 / norm value.
    pub normalized_p0: f64,
    pub normalized_p0_se: f64,
    /// Mean mismatch over the SNPs called in both (or of the window P0s
    /// with `window_estimate`).
    pub p0: f64,
    /// Jackknife standard error of the mean of the window P0s.
    pub p0_se: f64,
    /// Parent-offspring / siblings / neither, for first-degree pairs with
    /// at least 10 000 effective SNPs and a decisive window share.
    pub first_degree_type: Option<FirstDegreeType>,
    /// Share (0 to 1) of the first-degree windows whose normalized P0 is
    /// outside [t3, t1]; reported for every pair.
    pub first_degree_windows_outside: f64,
    /// SNPs called in both (in window mode, summed over the windows).
    pub overlap_snps: u64,
    /// overlap_snps × norm value: the expected mismatches.
    pub effective_snps: f64,
    /// 1 − normalized P0.
    pub kinship: f64,
    /// 95 % interval: kinship ± 1.96 SE (normalized P0 ± 1.96 SE with
    /// `ci_of_normalized_p0`).
    pub kinship_ci: (f64, f64),
}

/// READ's result: the norm value and every pair with a SNP called in both.
#[derive(Debug, Clone)]
pub struct ReadResult {
    pub norm_value: f64,
    pub pairs: Vec<ReadPair>,
}

/// Expected mismatches needed for the third-degree class.
const THIRD_DEGREE_MIN_EFFECTIVE: f64 = 3000.0;
/// Expected mismatches needed for the first-degree type.
const FIRST_TYPE_MIN_EFFECTIVE: f64 = 10_000.0;

/// The windows of `size` base pairs over SNPs in input order, as inclusive
/// index pairs: a new window at each change of ⌊position / size⌋ and (unless
/// `windows_ignore_chromosomes`) of chromosome; with `windows_overlap`
/// every window but the first also takes the last SNP of the one before.
fn windows(
    positions: &[u64],
    chroms: &[u32],
    size: u64,
    compat: ReadCompat,
) -> Vec<(usize, usize)> {
    let len = positions.len();
    if len == 0 {
        return Vec::new();
    }
    let mut starts = vec![0];
    for k in 1..len {
        let bin_changes = positions[k] / size != positions[k - 1] / size;
        let chrom_changes = !compat.windows_ignore_chromosomes && chroms[k] != chroms[k - 1];
        if bin_changes || chrom_changes {
            starts.push(k);
        }
    }
    (0..starts.len())
        .map(|r| {
            let end = starts.get(r + 1).map_or(len - 1, |s| s - 1);
            let start = if compat.windows_overlap && r > 0 {
                starts[r] - 1
            } else {
                starts[r]
            };
            (start, end)
        })
        .collect()
}

/// Delete-one jackknife standard error of the mean of `values`; with
/// `skip_last` the last value is never left out (the factor still uses
/// all of them).
fn jackknife_se(values: &[f64], skip_last: bool) -> f64 {
    let g = values.len();
    if g == 0 {
        return f64::NAN;
    }
    let sum: f64 = values.iter().sum();
    let gm1 = (g - 1) as f64;
    let kept = if skip_last { g - 1 } else { g };
    let loo: Vec<f64> = values[..kept].iter().map(|x| (sum - x) / gm1).collect();
    let q = if loo.is_empty() {
        0.0
    } else {
        let m = mean(&loo);
        loo.iter().map(|l| (l - m) * (l - m)).sum::<f64>()
    };
    (gm1 / g as f64 * q).sqrt()
}

/// The arithmetic mean, summed in order.
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// The median (the mean of the two middle values for an even count).
fn median(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        f64::midpoint(v[n / 2 - 1], v[n / 2])
    }
}

/// A pair's P0 before normalisation.
struct RawPair {
    i: usize,
    j: usize,
    p0: f64,
    p0_se: f64,
    overlap_snps: u64,
    /// P0 of each first-degree window with a SNP called in both.
    first_windows_p0: Vec<f64>,
}

fn raw_pair(
    bits: &PairBits,
    (i, j): (usize, usize),
    jack: &[counts::Segments],
    first: &[counts::Segments],
    options: &ReadOptions,
) -> Option<RawPair> {
    let total = bits.total();
    if total.n == 0 {
        return None;
    }
    let window_counts = bits.windows(jack);
    let window_p0: Vec<f64> = window_counts.iter().filter_map(|c| c.p0()).collect();
    let (p0, p0_se, overlap_snps) = if options.window_estimate {
        let overlap = window_counts.iter().map(|c: &Count| c.n).sum();
        (mean(&window_p0), jackknife_se(&window_p0, false), overlap)
    } else {
        let se = jackknife_se(&window_p0, options.compat.jackknife_skips_last_window);
        (total.p0()?, se, total.n)
    };
    let first_windows_p0 = bits.windows(first).iter().filter_map(|c| c.p0()).collect();
    Some(RawPair {
        i,
        j,
        p0,
        p0_se,
        overlap_snps,
        first_windows_p0,
    })
}

fn classify(x: f64, effective: f64, t: &[f64; 4]) -> ReadRelationship {
    if x > t[0] {
        ReadRelationship::Unrelated
    } else if x >= t[1] {
        if effective >= THIRD_DEGREE_MIN_EFFECTIVE {
            ReadRelationship::ThirdDegree
        } else {
            ReadRelationship::UnrelatedOrThirdDegree
        }
    } else if x >= t[2] {
        ReadRelationship::SecondDegree
    } else if x >= t[3] {
        ReadRelationship::FirstDegree
    } else {
        ReadRelationship::IdenticalOrTwins
    }
}

/// The class's (upper, lower) boundaries; unrelated pairs measure the
/// lower one to t1.
fn boundaries(rel: ReadRelationship, t: &[f64; 4]) -> (Option<f64>, Option<f64>) {
    match rel {
        ReadRelationship::Unrelated => (None, Some(t[1])),
        ReadRelationship::UnrelatedOrThirdDegree | ReadRelationship::ThirdDegree => {
            (Some(t[0]), Some(t[1]))
        }
        ReadRelationship::SecondDegree => (Some(t[1]), Some(t[2])),
        ReadRelationship::FirstDegree => (Some(t[2]), Some(t[3])),
        ReadRelationship::IdenticalOrTwins => (Some(t[3]), None),
    }
}

fn first_degree_type(rel: ReadRelationship, effective: f64, q: f64) -> Option<FirstDegreeType> {
    if rel != ReadRelationship::FirstDegree || effective < FIRST_TYPE_MIN_EFFECTIVE {
        return None;
    }
    if q <= 0.3 {
        Some(FirstDegreeType::ParentOffspring)
    } else if (0.35..0.6).contains(&q) {
        Some(FirstDegreeType::Siblings)
    } else if q >= 0.6 {
        Some(FirstDegreeType::NotApplicable)
    } else {
        None
    }
}

fn finish(raw: RawPair, norm: f64, options: &ReadOptions) -> ReadPair {
    let t = &options.thresholds;
    let x = raw.p0 / norm;
    let se_n = raw.p0_se / norm;
    let effective = raw.overlap_snps as f64 * norm;
    let relationship = classify(x, effective, t);
    let (upper, lower) = boundaries(relationship, t);
    let outside = raw
        .first_windows_p0
        .iter()
        .filter(|&&p| {
            let y = p / norm;
            !(t[3] <= y && y <= t[1])
        })
        .count();
    let q = outside as f64 / raw.first_windows_p0.len() as f64;
    let kinship = 1.0 - x;
    let centre = if options.compat.ci_of_normalized_p0 {
        x
    } else {
        kinship
    };
    ReadPair {
        i: raw.i,
        j: raw.j,
        relationship,
        z_upper: upper.map(|u| (u - x) / se_n),
        z_lower: lower.map(|b| (x - b) / se_n),
        normalized_p0: x,
        normalized_p0_se: se_n,
        p0: raw.p0,
        p0_se: raw.p0_se,
        first_degree_type: first_degree_type(relationship, effective, q),
        first_degree_windows_outside: q,
        overlap_snps: raw.overlap_snps,
        effective_snps: effective,
        kinship,
        kinship_ci: (centre - 1.96 * se_n, centre + 1.96 * se_n),
    }
}

fn validate(
    g: &BitGenotypes,
    positions: &[u64],
    chroms: &[u32],
    options: &ReadOptions,
) -> Result<()> {
    let m = g.n_snps();
    for (what, len) in [("positions", positions.len()), ("chroms", chroms.len())] {
        if len != m {
            return Err(PopmixError::DimensionMismatch {
                context: format!("READ kinship {what}"),
                expected: m.to_string(),
                actual: len.to_string(),
            });
        }
    }
    for (name, size) in [
        ("window_size", options.window_size),
        ("first_degree_window_size", options.first_degree_window_size),
    ] {
        if size == 0 {
            return Err(PopmixError::parameter(name, "must be at least 1 base pair"));
        }
    }
    if let ReadNorm::Value(v) = options.norm {
        if !v.is_finite() || v <= 0.0 {
            return Err(PopmixError::parameter(
                "norm_value",
                format!("must be a positive number, got {v}"),
            ));
        }
    }
    let t = &options.thresholds;
    if !t.iter().all(|v| v.is_finite()) || !t.windows(2).all(|w| w[0] > w[1]) {
        return Err(PopmixError::parameter(
            "thresholds",
            format!("must be finite and strictly decreasing, got {t:?}"),
        ));
    }
    if g.n_samples() < 2 {
        return Err(PopmixError::too_few(
            "READ kinship",
            "samples",
            2,
            g.n_samples(),
        ));
    }
    Ok(())
}

/// READ kinship between every pair of samples of `g`, the SNPs at
/// `positions` (base pairs) on `chroms` (only a change from one SNP to the
/// next matters), in input order.
///
/// Pairs with no SNP called in both are left out. The pairs come in the
/// order (0, 1), (0, 2), …, (n−2, n−1), whatever the thread count.
pub fn read_kinship(
    g: &BitGenotypes,
    positions: &[u64],
    chroms: &[u32],
    options: &ReadOptions,
) -> Result<ReadResult> {
    validate(g, positions, chroms, options)?;
    let jack = segments(&windows(
        positions,
        chroms,
        options.window_size,
        options.compat,
    ));
    let first = segments(&windows(
        positions,
        chroms,
        options.first_degree_window_size,
        options.compat,
    ));
    let n = g.n_samples();
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();
    let raw: Vec<RawPair> = pairs
        .par_iter()
        .filter_map(|&(i, j)| raw_pair(&PairBits::new(g, i, j), (i, j), &jack, &first, options))
        .collect();
    let norm = match options.norm {
        ReadNorm::Value(v) => v,
        _ if raw.is_empty() => return Err(PopmixError::data(
            "READ: no pair of samples has a SNP called in both, so there is no P0 to normalize by",
        )),
        ReadNorm::Median => median(&raw.iter().map(|p| p.p0).collect::<Vec<_>>()),
        ReadNorm::Mean => mean(&raw.iter().map(|p| p.p0).collect::<Vec<_>>()),
        ReadNorm::Max => raw.iter().map(|p| p.p0).fold(f64::NEG_INFINITY, f64::max),
    };
    if norm == 0.0 || !norm.is_finite() {
        return Err(PopmixError::data(format!(
            "READ: the normalization value is {norm}: are most samples duplicates?"
        )));
    }
    let pairs = raw.into_iter().map(|p| finish(p, norm, options)).collect();
    Ok(ReadResult {
        norm_value: norm,
        pairs,
    })
}

#[cfg(test)]
mod tests;
