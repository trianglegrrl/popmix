//! ALDER (Loh et al. 2013): weighted LD curves for dating admixture and
//! testing for it.
//!
//! With two references A and B, the test population C's LD is weighted by
//! the allele frequency difference between A and B; with one reference A, by
//! A's frequencies through ALDER's unbiased single-reference estimator. The
//! curve is fit from the distance at which C's LD stops being correlated
//! with the references' (at least 0.5 cM) as `A exp(-n d) + C`, the constant
//! taken from SNP pairs on different chromosomes, with standard errors from a
//! jackknife over chromosomes. With two references the admixture test checks
//! that both one-reference curves and the two-reference curve are
//! significant.
//!
//! Genotypes are read as called (a pseudo-haploid sample counts as
//! homozygous), as ALDER reads them.

pub mod fit;
pub mod ldcorr;
pub mod wld;

use rayon::prelude::*;

use crate::error::{PopmixError, Result};
use crate::genotypes::{Chromosome, Snp};
use crate::pca::PackedGenotypes;
pub use fit::{fit_decay, CurveFit, DecayFitParams, Jackknifed};

/// Which ALDER quirks to reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AlderCompat {
    /// ALDER-0001: the inter-chromosome affine term divides by the pairs of
    /// complete SNPs, counting SNPs without a weight (a reference with no
    /// data there) that its numerator leaves out.
    pub affine_counts_unweighted_snps: bool,
}

impl AlderCompat {
    pub const CORRECTED: AlderCompat = AlderCompat {
        affine_counts_unweighted_snps: false,
    };
    pub const V1_03: AlderCompat = AlderCompat {
        affine_counts_unweighted_snps: true,
    };
}

#[derive(Debug, Clone)]
pub struct AlderOptions {
    /// Distance bin (Morgans).
    pub binsize: f64,
    /// Largest distance fit (Morgans).
    pub maxdis: f64,
    /// Test samples that must be called at a SNP for the two-reference curve.
    pub mincount: usize,
    /// Fit from this distance instead of finding where LD stops being
    /// correlated with the references'.
    pub mindis: Option<f64>,
    pub compat: AlderCompat,
}

impl Default for AlderOptions {
    fn default() -> Self {
        AlderOptions {
            binsize: 0.0005,
            maxdis: 0.5,
            mincount: 4,
            mindis: None,
            compat: AlderCompat::CORRECTED,
        }
    }
}

/// The smallest distance ALDER fits from.
pub const DEFAULT_FIT_START: f64 = 0.005;

/// A weighted LD curve: distance (Morgans; the last row, at infinity, the
/// inter-chromosome affine term if computed), weighted LD and SNP pairs.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct WeightedLdCurve {
    pub d: Vec<f64>,
    pub wld: Vec<f64>,
    pub count: Vec<f64>,
}

/// One weighted LD computation and its fits.
#[derive(Debug, Clone)]
pub struct AlderRun {
    /// The references whose weights were used (indices into the input).
    pub refs: Vec<usize>,
    pub fit_start: f64,
    pub curve: WeightedLdCurve,
    /// The curve with each chromosome left out (labels in `AlderResult`).
    pub jackknife_curves: Vec<WeightedLdCurve>,
    /// Fits from fit_start - 0.2 cM to fit_start + 0.2 cM in 0.1 cM steps.
    pub fits: Vec<CurveFit>,
    /// The fit from fit_start.
    pub fit_index: usize,
}

impl AlderRun {
    pub fn fit(&self) -> &CurveFit {
        &self.fits[self.fit_index]
    }
}

/// A decay-rate comparison between two curves.
#[derive(Debug, Clone, PartialEq)]
pub struct DecayDiff {
    pub z: f64,
    pub percent: f64,
}

/// ALDER's test for admixture with two references.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmixtureTest {
    /// Both one-reference curves significant.
    pub pretest: bool,
    pub tworef_significant: bool,
    /// One-reference A and B minus two-reference, and B minus A.
    pub decay_diffs: [DecayDiff; 3],
    /// All within 25%.
    pub consistent: bool,
    pub succeeds: bool,
    pub z: f64,
    pub p: f64,
}

#[derive(Debug, Clone)]
pub struct AlderResult {
    /// Chromosome labels, in order (the jackknife replicates).
    pub chromosomes: Vec<String>,
    /// Per reference, the distance to which LD is correlated (Morgans;
    /// infinite if past 2 cM), before the 0.5 cM floor; None with `mindis`.
    pub ld_corr_stops: Vec<Option<f64>>,
    pub run: AlderRun,
    /// With two references: the one-reference run with each.
    pub oneref_runs: Vec<AlderRun>,
    pub test: Option<AdmixtureTest>,
    /// With one reference: the lower bound on the mixture fraction (value,
    /// jackknife SE).
    pub mixture_fraction: Option<(f64, f64)>,
}

/// Per-SNP data of the test population and the chromosome layout.
struct Prepared<'a> {
    test: &'a PackedGenotypes,
    n: usize,
    chroms: Vec<Vec<usize>>,
    labels: Vec<String>,
    pos: Vec<f64>,
    missing: Vec<usize>,
    sum: Vec<f64>,
    sum2: Vec<f64>,
}

fn chrom_label(c: &Chromosome) -> String {
    c.to_string()
}

fn prepare<'a>(test: &'a PackedGenotypes, snps: &[Snp]) -> Result<Prepared<'a>> {
    if snps.len() != test.n_snps() {
        return Err(PopmixError::DimensionMismatch {
            context: "ALDER SNPs".to_string(),
            expected: format!("{} SNPs", test.n_snps()),
            actual: format!("{}", snps.len()),
        });
    }
    let mut chroms: Vec<Vec<usize>> = Vec::new();
    let mut labels = Vec::new();
    for (s, snp) in snps.iter().enumerate() {
        let new = s == 0 || snp.chromosome != snps[s - 1].chromosome;
        if !new && snp.position_morgans < snps[s - 1].position_morgans {
            return Err(PopmixError::Numerical(format!(
                "SNPs must be sorted by genetic position (SNP {} before {})",
                snps[s - 1].id,
                snp.id
            )));
        }
        if new {
            let label = chrom_label(&snp.chromosome);
            if labels.contains(&label) {
                return Err(PopmixError::Numerical(format!(
                    "SNPs of chromosome {label} are not contiguous"
                )));
            }
            labels.push(label);
            chroms.push(Vec::new());
        }
        chroms.last_mut().expect("chromosome").push(s);
    }
    let n = test.n_samples();
    let mut missing = vec![0; snps.len()];
    let mut sum = vec![0.0; snps.len()];
    let mut sum2 = vec![0.0; snps.len()];
    for s in 0..snps.len() {
        for i in 0..n {
            let g = test.get(s, i);
            if g > 2 {
                missing[s] += 1;
            } else {
                sum[s] += g as f64;
                sum2[s] += (g * g) as f64;
            }
        }
    }
    Ok(Prepared {
        test,
        n,
        chroms,
        labels,
        pos: snps.iter().map(|s| s.position_morgans).collect(),
        missing,
        sum,
        sum2,
    })
}

/// Allele frequency of the counted allele per SNP (NaN without data).
fn frequencies(g: &PackedGenotypes) -> Vec<f64> {
    (0..g.n_snps())
        .map(|s| {
            let (mut t, mut c) = (0u32, 0u32);
            for i in 0..g.n_samples() {
                let x = g.get(s, i);
                if x <= 2 {
                    t += x as u32;
                    c += 1;
                }
            }
            0.5 * t as f64 / c as f64
        })
        .collect()
}

/// Compute a weighted LD curve and fit it.
fn run_curve(
    p: &Prepared,
    num_refs: usize,
    refs: Vec<usize>,
    weights: &[f64],
    fit_start: f64,
    options: &AlderOptions,
) -> Result<(AlderRun, Vec<bool>)> {
    if num_refs == 1 && p.n < 4 {
        return Err(PopmixError::Numerical(
            "the one-reference curve needs at least 4 test samples".to_string(),
        ));
    }
    let numbins = (options.maxdis / options.binsize) as usize;
    let maxmissing = if num_refs == 2 {
        p.n.saturating_sub(options.mincount)
    } else {
        0
    };
    let ignore: Vec<bool> = (0..p.pos.len())
        .map(|s| p.missing[s] > maxmissing || weights[s].is_nan())
        .collect();
    let per_chrom: Vec<wld::ChromCurve> = p
        .chroms
        .par_iter()
        .map(|snps| {
            let first = p.pos[snps[0]];
            let bins: Vec<usize> = snps
                .iter()
                .map(|&s| ((p.pos[s] - first) / options.binsize) as usize)
                .collect();
            let take = |v: &[f64]| snps.iter().map(|&s| v[s]).collect::<Vec<f64>>();
            let (w, sum, sum2) = (take(weights), take(&p.sum), take(&p.sum2));
            let missing: Vec<usize> = snps.iter().map(|&s| p.missing[s]).collect();
            let ign: Vec<bool> = snps.iter().map(|&s| ignore[s]).collect();
            let data = wld::ChromData {
                geno: p.test,
                snps,
                bins: &bins,
                weights: &w,
                missing: &missing,
                sum: &sum,
                sum2: &sum2,
                ignore: &ign,
                count_ignored_in_affine: options.compat.affine_counts_unweighted_snps,
            };
            wld::run_chrom(&data, p.n, num_refs, numbins)
        })
        .collect();

    let nc = p.chroms.len();
    let jackknife = nc > 1;
    let inter_chrom = if jackknife { nc > 2 } else { nc > 1 };
    let curve_without = |skip: Option<usize>| -> WeightedLdCurve {
        let mut c = WeightedLdCurve {
            d: (0..numbins).map(|b| b as f64 * options.binsize).collect(),
            wld: vec![0.0; numbins],
            count: vec![0.0; numbins],
        };
        for (k, ch) in per_chrom.iter().enumerate() {
            if Some(k) == skip {
                continue;
            }
            for b in 1..numbins {
                c.wld[b] += ch.wld[b];
                c.count[b] += ch.count[b];
            }
        }
        c.d[0] = 0.0;
        for b in 1..numbins {
            c.wld[b] /= c.count[b];
        }
        if inter_chrom {
            let affs: Vec<&wld::AffineData> = per_chrom
                .iter()
                .enumerate()
                .filter(|(k, _)| Some(*k) != skip)
                .map(|(_, ch)| &ch.affine)
                .collect();
            let (aff, pairs) = wld::inter_chrom_affine(&affs, p.n, num_refs);
            c.d.push(f64::INFINITY);
            c.wld.push(aff);
            c.count.push(pairs);
        }
        c
    };
    let curve = curve_without(None);
    let jackknife_curves: Vec<WeightedLdCurve> = (0..nc).map(|k| curve_without(Some(k))).collect();

    let mut fits = Vec::new();
    let mut fit_index = 0;
    if fit_start.is_finite() {
        // as ALDER: starts accumulated in steps of 0.001
        let mut mindis = fit_start - 0.002;
        while mindis <= fit_start + 0.0021 {
            if (mindis - fit_start).abs() < 1e-9 {
                fit_index = fits.len();
            }
            fits.push(fit_all(
                &curve,
                &jackknife_curves,
                jackknife,
                mindis,
                options.maxdis,
            ));
            mindis += 0.001;
        }
    }
    Ok((
        AlderRun {
            refs,
            fit_start,
            curve,
            jackknife_curves,
            fits,
            fit_index,
        },
        ignore,
    ))
}

/// Fit the full curve and each jackknife curve from `mindis`.
fn fit_all(
    full: &WeightedLdCurve,
    jacks: &[WeightedLdCurve],
    jackknife: bool,
    mindis: f64,
    maxdis: f64,
) -> CurveFit {
    let one = |c: &WeightedLdCurve| -> (DecayFitParams, bool) {
        let bmin =
            c.d.iter()
                .position(|&d| d >= mindis - 1e-9)
                .unwrap_or(c.d.len());
        let inter = c.d.last().is_some_and(|d| d.is_infinite());
        let end = if inter { c.d.len() - 1 } else { c.d.len() };
        let affine = inter.then(|| *c.wld.last().expect("affine row"));
        let fit = if bmin < end {
            fit_decay(&c.d[bmin..end], &c.wld[bmin..end], affine)
        } else {
            DecayFitParams {
                decay: f64::NAN,
                amp_exp: f64::NAN,
                amp_aff: f64::NAN,
                success: false,
            }
        };
        (fit, inter)
    };
    let (full_fit, inter) = one(full);
    let mut failures = 0;
    let reps: Vec<DecayFitParams> = if jackknife {
        jacks
            .par_iter()
            .map(|c| one(c).0)
            .collect::<Vec<_>>()
            .into_iter()
            .map(|f| {
                if f.success {
                    f
                } else {
                    failures += 1;
                    DecayFitParams {
                        decay: f64::INFINITY,
                        amp_exp: f64::INFINITY,
                        amp_aff: f64::INFINITY,
                        success: false,
                    }
                }
            })
            .collect()
    } else {
        Vec::new()
    };
    let jk = |f: &dyn Fn(&DecayFitParams) -> f64| Jackknifed {
        value: f(&full_fit),
        reps: reps.iter().map(f).collect(),
    };
    let tot = |p: &DecayFitParams| {
        if p.success || !p.amp_exp.is_infinite() {
            p.amp_exp + p.amp_aff / 2.0
        } else {
            f64::INFINITY
        }
    };
    CurveFit {
        mindis,
        maxdis,
        inter_chrom_affine: inter,
        decay: jk(&|p| p.decay),
        amp_tot: jk(&tot),
        amp_exp: jk(&|p| p.amp_exp),
        amp_aff: jk(&|p| p.amp_aff),
        failures,
    }
}

fn decay_diff(a: &CurveFit, b: &CurveFit) -> DecayDiff {
    let d = a.decay.minus(&b.decay);
    let mean = (a.decay.value + b.decay.value) / 2.0;
    DecayDiff {
        z: d.value / d.se(),
        percent: 100.0 * d.value / mean,
    }
}

/// The mixture-fraction lower bound from a one-reference fit (ALDER
/// `mix_frac_bound`): amp / f2(C, A)^2 = 2a / (1 - a), solved for each
/// jackknife replicate.
fn mixture_fraction(
    p: &Prepared,
    reference: &PackedGenotypes,
    ignore: &[bool],
    fit: &CurveFit,
) -> (f64, f64) {
    let nc = p.chroms.len();
    let mut f2 = vec![0.0; nc];
    let mut count = vec![0.0; nc];
    let alleles = |g: &PackedGenotypes, s: usize| {
        let (mut a, mut b) = (0.0, 0.0);
        for i in 0..g.n_samples() {
            let x = g.get(s, i);
            if x <= 2 {
                a += x as f64;
                b += 2.0 - x as f64;
            }
        }
        (a, b)
    };
    for (c, snps) in p.chroms.iter().enumerate() {
        for &s in snps {
            if ignore[s] {
                continue;
            }
            let (a1, b1) = alleles(p.test, s);
            let (a2, b2) = alleles(reference, s);
            if a1 + b1 <= 1.0 || a2 + b2 <= 1.0 {
                continue;
            }
            let (p1, p2) = (a1 / (a1 + b1), a2 / (a2 + b2));
            let bias1 = a1 * b1 / ((a1 + b1) * (a1 + b1) * (a1 + b1 - 1.0));
            let bias2 = a2 * b2 / ((a2 + b2) * (a2 + b2) * (a2 + b2 - 1.0));
            f2[c] += (p1 - p2) * (p1 - p2) - bias1 - bias2;
            count[c] += 1.0;
        }
    }
    let f2_without = |skip: Option<usize>| {
        let (mut s, mut n) = (0.0, 0.0);
        for c in 0..nc {
            if Some(c) != skip {
                s += f2[c];
                n += count[c];
            }
        }
        s / n
    };
    let alpha = |amp: f64, f2: f64| {
        let y = amp / f2.powi(2);
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        while hi - lo > 1e-6 {
            let a = (lo + hi) / 2.0;
            if 2.0 * a / (1.0 - a) < y {
                lo = a;
            } else {
                hi = a;
            }
        }
        (lo + hi) / 2.0
    };
    let j = Jackknifed {
        value: alpha(fit.amp_tot.value, f2_without(None)),
        reps: (0..nc)
            .map(|c| alpha(fit.amp_tot.reps[c], f2_without(Some(c))))
            .collect(),
    };
    (j.value, j.se())
}

/// ALDER with one or two references. `test` and each of `refs` hold the
/// samples of those populations at the same SNPs, `snps` (autosomes, sorted
/// by chromosome and genetic position).
pub fn alder(
    test: &PackedGenotypes,
    refs: &[&PackedGenotypes],
    snps: &[Snp],
    options: &AlderOptions,
) -> Result<AlderResult> {
    if refs.is_empty() || refs.len() > 2 {
        return Err(PopmixError::Numerical(format!(
            "ALDER takes one or two reference populations, got {}",
            refs.len()
        )));
    }
    if let Some(r) = refs.iter().find(|r| r.n_snps() != test.n_snps()) {
        return Err(PopmixError::DimensionMismatch {
            context: "ALDER reference SNPs".to_string(),
            expected: format!("{}", test.n_snps()),
            actual: format!("{}", r.n_snps()),
        });
    }
    if options.mincount < 2 || options.mincount > test.n_samples() {
        return Err(PopmixError::Numerical(format!(
            "mincount must be between 2 and the {} test samples",
            test.n_samples()
        )));
    }
    let p = prepare(test, snps)?;
    let freqs: Vec<Vec<f64>> = refs.iter().map(|r| frequencies(r)).collect();

    let ld_corr_stops: Vec<Option<f64>> = match options.mindis {
        Some(_) => vec![None; refs.len()],
        None => refs
            .iter()
            .map(|r| {
                Some(ldcorr::ld_corr_stop(
                    test,
                    r,
                    &p.chroms,
                    &p.pos,
                    options.binsize,
                ))
            })
            .collect(),
    };
    let fit_starts: Vec<f64> = ld_corr_stops
        .iter()
        .map(|s| match (options.mindis, s) {
            (Some(m), _) => m,
            (None, Some(d)) if *d < DEFAULT_FIT_START => DEFAULT_FIT_START,
            (None, Some(d)) => *d,
            (None, None) => unreachable!(),
        })
        .collect();
    let fit_start = fit_starts.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    if refs.len() == 1 {
        let (run, ignore) = run_curve(&p, 1, vec![0], &freqs[0], fit_start, options)?;
        let mixture = (p.chroms.len() > 1 && !run.fits.is_empty())
            .then(|| mixture_fraction(&p, refs[0], &ignore, run.fit()));
        return Ok(AlderResult {
            chromosomes: p.labels.clone(),
            ld_corr_stops,
            run,
            oneref_runs: Vec::new(),
            test: None,
            mixture_fraction: mixture,
        });
    }
    let weights: Vec<f64> = freqs[0].iter().zip(&freqs[1]).map(|(a, b)| a - b).collect();
    let (run, _) = run_curve(&p, 2, vec![0, 1], &weights, fit_start, options)?;
    let mut oneref_runs = Vec::new();
    for (r, f) in freqs.iter().enumerate() {
        oneref_runs.push(run_curve(&p, 1, vec![r], f, fit_starts[r], options)?.0);
    }
    let test_result = (p.chroms.len() > 1
        && !run.fits.is_empty()
        && oneref_runs.iter().all(|r| !r.fits.is_empty()))
    .then(|| {
        let (two, a, b) = (run.fit(), oneref_runs[0].fit(), oneref_runs[1].fit());
        let pretest = a.is_significant(1.0) && b.is_significant(1.0);
        let tworef_significant = two.is_significant(1.0);
        let diffs = [decay_diff(a, two), decay_diff(b, two), decay_diff(b, a)];
        let consistent = diffs.iter().all(|d| d.percent.abs() < 25.0);
        AdmixtureTest {
            pretest,
            tworef_significant,
            consistent,
            succeeds: pretest && tworef_significant,
            z: two.decay.z().min(two.amp_exp.z()),
            p: two.p_value(1.0),
            decay_diffs: diffs,
        }
    });
    Ok(AlderResult {
        chromosomes: p.labels,
        ld_corr_stops,
        run,
        oneref_runs,
        test: test_result,
        mixture_fraction: None,
    })
}
