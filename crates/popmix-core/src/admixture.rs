//! Maximum-likelihood ancestry proportions, as ADMIXTURE (Alexander,
//! Novembre and Lange 2009).
//!
//! Each sample i has ancestry proportions q_i over K populations and each SNP
//! j an allele frequency p_kj in each; with h_ij copies of the allele whose
//! frequency is p (out of 2), the log-likelihood is
//!
//! ```text
//! L = sum_ij h_ij ln(sum_k q_ik p_kj) + (2 - h_ij) ln(sum_k q_ik (1 - p_kj))
//! ```
//!
//! over the called genotypes, exactly ADMIXTURE's (no binomial constant).
//! As ADMIXTURE, p is the frequency of the allele the genotype codes do not
//! count, and Q and P are kept within [1e-5, 1 - 1e-5].
//!
//! The maximum is found by block relaxation, as ADMIXTURE: a Newton step on
//! each sample's Q with P fixed (a small quadratic program on the simplex),
//! then on each SNP's P with Q fixed (a box-constrained one), each halved
//! until its block's log-likelihood does not decrease; the iteration is
//! accelerated with SQUAREM (Varadhan and Roland 2008), where ADMIXTURE uses
//! quasi-Newton steps, so the path differs and so can the local maximum
//! reached. With P fixed
//! (projection) the problem is convex and the Q found is the unique maximum.

mod newton;
mod qp;

use newton::*;
use qp::*;

use rayon::prelude::*;

use nalgebra::{DMatrix, DVector};

use crate::error::{PopmixError, Result};
use crate::pca::PackedGenotypes;

/// ADMIXTURE's bounds on Q and P.
pub const BOUND: f64 = 1e-5;

/// Options for fitting.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct AdmixtureOptions {
    /// Stop when an iteration improves the log-likelihood by less than this
    /// (ADMIXTURE `-C`, default 1e-4), or by less than 1e-13 of it.
    pub tolerance: f64,
    pub max_iterations: usize,
    /// Seed for the random starting values.
    pub seed: u64,
}

impl Default for AdmixtureOptions {
    fn default() -> Self {
        AdmixtureOptions {
            tolerance: 1e-4,
            max_iterations: 10_000,
            seed: 1,
        }
    }
}

/// A fit: Q (samples x K, row-major), P (SNPs x K, row-major).
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct AdmixtureFit {
    pub k: usize,
    pub q: Vec<f64>,
    pub p: Vec<f64>,
    pub loglik: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// What is estimated: everything, Q with P fixed, or P and the Q of
/// unlabeled samples with labeled samples fixed to their population.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Free,
    FixedP,
}

/// The genotype data a fit uses: genotypes, and optionally the entries held
/// out (for cross-validation).
struct Data<'a> {
    geno: &'a PackedGenotypes,
    held_out: Option<&'a HeldOut>,
}

impl Data<'_> {
    /// The genotype code g (0, 1, 2; h = 2 - g copies of the allele P refers
    /// to); None if missing or held out.
    #[inline]
    fn g(&self, snp: usize, sample: usize) -> Option<u8> {
        let g = self.geno.get(snp, sample);
        if g > 2 || self.held_out.is_some_and(|m| m.contains(snp, sample)) {
            None
        } else {
            Some(g)
        }
    }
}

/// A genotype's log-likelihood term h ln(a) + (2 - h) ln(b), h = 2 - g, with
/// one logarithm: a homozygote's other term is zero, a heterozygote's two
/// are ln(ab). The logarithms are most of the fit's time.
#[inline]
fn term(g: u8, a: f64, b: f64) -> f64 {
    // a^h b^(2-h), chosen without a branch: genotypes vary from sample to
    // sample, and a mispredicted branch costs more than the logarithm
    [a * a, a * b, b * b][g as usize].ln()
}

/// Genotype entries held out of a fit.
pub struct HeldOut {
    n_samples: usize,
    bits: Vec<u64>,
}

impl HeldOut {
    fn contains(&self, snp: usize, sample: usize) -> bool {
        let i = snp * self.n_samples + sample;
        self.bits[i / 64] >> (i % 64) & 1 == 1
    }
}

fn clamp(x: f64) -> f64 {
    x.clamp(BOUND, 1.0 - BOUND)
}

/// Q rows onto the simplex with every entry at least `BOUND`.
fn normalize_q(q: &mut [f64], k: usize) {
    for row in q.chunks_mut(k) {
        for x in row.iter_mut() {
            *x = x.max(BOUND);
        }
        let s: f64 = row.iter().sum();
        row.iter_mut().for_each(|x| *x /= s);
    }
}

/// Relative rounding error of a log-likelihood summed over many genotypes:
/// a change smaller than this fraction of it is noise.
const ROUNDING: f64 = 1e-13;

/// SNPs per parallel chunk. Fixed, and partial sums are added in order, so
/// results do not depend on the number of threads.
const CHUNK: usize = 256;

/// Log-likelihood at (q, p).
fn loglik(d: &Data, q: &[f64], p: &[f64], k: usize) -> f64 {
    let n = d.geno.n_samples();
    let m = d.geno.n_snps();
    let parts: Vec<f64> = (0..m.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| {
            let mut l = 0.0;
            for j in c * CHUNK..((c + 1) * CHUNK).min(m) {
                let pj = &p[j * k..(j + 1) * k];
                for i in 0..n {
                    if let Some(g) = d.g(j, i) {
                        let qi = &q[i * k..(i + 1) * k];
                        let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                        let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                        l += term(g, a, b);
                    }
                }
            }
            l
        })
        .collect();
    parts.iter().sum()
}

/// One block-relaxation iteration (ADMIXTURE's method): a Newton step on Q
/// with P fixed, then on P with the new Q fixed; with the log-likelihood
/// there.
fn block_step(
    d: &Data,
    q: &[f64],
    p: &[f64],
    k: usize,
    mode: Mode,
    fixed_q: &[bool],
) -> (Vec<f64>, Vec<f64>, f64) {
    let q1 = q_newton(d, q, p, k, fixed_q);
    if mode == Mode::FixedP {
        let ll = loglik(d, &q1, p, k);
        (q1, p.to_vec(), ll)
    } else {
        let (p1, ll) = p_newton(d, &q1, p, k);
        (q1, p1, ll)
    }
}

/// Block relaxation from (q, p), accelerated with SQUAREM, until the
/// log-likelihood improves by less than the tolerance.
fn fit(
    d: &Data,
    mut q: Vec<f64>,
    mut p: Vec<f64>,
    k: usize,
    mode: Mode,
    fixed_q: &[bool],
    options: &AdmixtureOptions,
) -> AdmixtureFit {
    let mut ll = loglik(d, &q, &p, k);
    let mut iterations = 0;
    let mut converged = false;
    let fix = |q: &mut Vec<f64>, p: &mut Vec<f64>, q0: &[f64]| {
        for (i, f) in fixed_q.iter().enumerate() {
            if *f {
                q[i * k..(i + 1) * k].copy_from_slice(&q0[i * k..(i + 1) * k]);
            }
        }
        normalize_q(q, k);
        p.iter_mut().for_each(|x| *x = clamp(*x));
    };
    while iterations < options.max_iterations {
        iterations += 1;
        let (q1, p1, _) = block_step(d, &q, &p, k, mode, fixed_q);
        let (q2, p2, l2) = block_step(d, &q1, &p1, k, mode, fixed_q);
        // r = theta1 - theta0, v = theta2 - 2 theta1 + theta0
        let diff = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x - y).collect::<Vec<f64>>();
        let (rq, rp) = (diff(&q1, &q), diff(&p1, &p));
        let (vq, vp) = (
            q2.iter()
                .zip(&q1)
                .zip(&rq)
                .map(|((a, b), r)| a - b - r)
                .collect::<Vec<f64>>(),
            p2.iter()
                .zip(&p1)
                .zip(&rp)
                .map(|((a, b), r)| a - b - r)
                .collect::<Vec<f64>>(),
        );
        let norm = |a: &[f64], b: &[f64]| (a.iter().chain(b).map(|x| x * x).sum::<f64>()).sqrt();
        let (nr, nv) = (norm(&rq, &rp), norm(&vq, &vp));
        let mut candidate = None;
        if nv > 0.0 {
            let alpha = -(nr / nv).max(1.0);
            let step = |x: &[f64], r: &[f64], v: &[f64]| -> Vec<f64> {
                x.iter()
                    .zip(r)
                    .zip(v)
                    .map(|((x, r), v)| x - 2.0 * alpha * r + alpha * alpha * v)
                    .collect()
            };
            let (mut qa, mut pa) = (step(&q, &rq, &vq), step(&p, &rp, &vp));
            fix(&mut qa, &mut pa, &q);
            // a stabilizing step from the extrapolated point
            let (qs, ps, ls) = block_step(d, &qa, &pa, k, mode, fixed_q);
            candidate = Some((qs, ps, ls));
        }
        let (qn, pn, ln) = match candidate {
            Some((qs, ps, ls)) if ls.is_finite() && ls >= l2 => (qs, ps, ls),
            _ => (q2, p2, l2),
        };
        let gain = ln - ll;
        q = qn;
        p = pn;
        ll = ln;
        // below 1e-13 of the log-likelihood, changes are rounding noise in
        // its sum over millions of genotypes
        if gain.abs() < options.tolerance.max(ROUNDING * ll.abs()) {
            converged = true;
            break;
        }
    }
    AdmixtureFit {
        k,
        q,
        p,
        loglik: ll,
        iterations,
        converged,
    }
}

/// A small deterministic generator (SplitMix64) for starting values and
/// cross-validation folds.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_start(n: usize, m: usize, k: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = Rng(seed);
    let mut q: Vec<f64> = (0..n * k).map(|_| rng.uniform() + 0.1).collect();
    normalize_q(&mut q, k);
    let p: Vec<f64> = (0..m * k).map(|_| 0.05 + 0.9 * rng.uniform()).collect();
    (q, p)
}

fn check(geno: &PackedGenotypes, k: usize) -> Result<()> {
    if k < 1 {
        return Err(PopmixError::parameter("K", "must be at least 1"));
    }
    if geno.n_samples() < 2 || geno.n_snps() < 1 {
        return Err(PopmixError::data(format!(
            "{} samples and {} SNPs: too few to fit",
            geno.n_samples(),
            geno.n_snps()
        )));
    }
    Ok(())
}

/// Unsupervised fit with K populations.
pub fn admixture(
    geno: &PackedGenotypes,
    k: usize,
    options: &AdmixtureOptions,
) -> Result<AdmixtureFit> {
    check(geno, k)?;
    let (q, p) = random_start(geno.n_samples(), geno.n_snps(), k, options.seed);
    let d = Data {
        geno,
        held_out: None,
    };
    Ok(fit(
        &d,
        q,
        p,
        k,
        Mode::Free,
        &vec![false; geno.n_samples()],
        options,
    ))
}

/// Supervised fit: `labels[i]` is sample i's population (0..K) or None to
/// estimate; labeled samples keep Q at their population (within the
/// bounds), as ADMIXTURE `--supervised`.
pub fn admixture_supervised(
    geno: &PackedGenotypes,
    k: usize,
    labels: &[Option<usize>],
    options: &AdmixtureOptions,
) -> Result<AdmixtureFit> {
    check(geno, k)?;
    if labels.len() != geno.n_samples() {
        return Err(PopmixError::DimensionMismatch {
            context: "supervised labels".to_string(),
            expected: format!("{} samples", geno.n_samples()),
            actual: format!("{}", labels.len()),
        });
    }
    if let Some(bad) = labels.iter().flatten().find(|&&l| l >= k) {
        return Err(PopmixError::out_of_range("label", *bad, k));
    }
    let (mut q, p) = random_start(geno.n_samples(), geno.n_snps(), k, options.seed);
    for (i, l) in labels.iter().enumerate() {
        if let Some(l) = l {
            for t in 0..k {
                q[i * k + t] = if t == *l { 1.0 } else { 0.0 };
            }
        }
    }
    normalize_q(&mut q, k);
    let fixed: Vec<bool> = labels.iter().map(|l| l.is_some()).collect();
    let d = Data {
        geno,
        held_out: None,
    };
    Ok(fit(&d, q, p, k, Mode::Free, &fixed, options))
}

/// Projection: Q of every sample with P fixed (SNPs x K, row-major), as
/// ADMIXTURE `-P`.
pub fn admixture_projection(
    geno: &PackedGenotypes,
    p: &[f64],
    k: usize,
    options: &AdmixtureOptions,
) -> Result<AdmixtureFit> {
    check(geno, k)?;
    if p.len() != geno.n_snps() * k {
        return Err(PopmixError::DimensionMismatch {
            context: "projection P".to_string(),
            expected: format!("{} x {k}", geno.n_snps()),
            actual: format!("{} values", p.len()),
        });
    }
    let n = geno.n_samples();
    let q = vec![1.0 / k as f64; n * k];
    let p: Vec<f64> = p.iter().map(|&x| clamp(x)).collect();
    let d = Data {
        geno,
        held_out: None,
    };
    Ok(fit(&d, q, p, k, Mode::FixedP, &vec![false; n], options))
}

/// Log-likelihood of given Q (samples x K) and P (SNPs x K).
pub fn admixture_loglik(geno: &PackedGenotypes, q: &[f64], p: &[f64], k: usize) -> Result<f64> {
    if q.len() != geno.n_samples() * k || p.len() != geno.n_snps() * k {
        return Err(PopmixError::DimensionMismatch {
            context: "admixture log-likelihood".to_string(),
            expected: format!("Q {} x {k}, P {} x {k}", geno.n_samples(), geno.n_snps()),
            actual: format!("{} and {} values", q.len(), p.len()),
        });
    }
    Ok(loglik(
        &Data {
            geno,
            held_out: None,
        },
        q,
        p,
        k,
    ))
}

/// Cross-validation error: the called genotypes are split into `folds` at
/// random; for each fold, Q and P are refit from `full` without it, and the
/// held-out genotypes predicted. The error is half the mean deviance of the
/// predictions, g ln(g / mu) + (2 - g) ln((2 - g) / (2 - mu)) with
/// mu = 2 sum_k q p (the expected count), as ADMIXTURE reports it.
pub fn admixture_cv(
    geno: &PackedGenotypes,
    full: &AdmixtureFit,
    folds: usize,
    options: &AdmixtureOptions,
) -> Result<f64> {
    // each genotype's fold is kept in a byte
    if !(2..=256).contains(&folds) {
        return Err(PopmixError::parameter(
            "folds",
            format!("cross-validation takes 2 to 256 folds, not {folds}"),
        ));
    }
    let (n, m, k) = (geno.n_samples(), geno.n_snps(), full.k);
    let mut rng = Rng(options.seed ^ 0xC0FF_EE00_D15E_A5E5);
    let fold_of: Vec<u8> = (0..n * m)
        .map(|_| (rng.next() % folds as u64) as u8)
        .collect();
    let (mut dev, mut count) = (0.0, 0usize);
    for f in 0..folds {
        let mut bits = vec![0u64; (n * m).div_ceil(64)];
        for (idx, &fo) in fold_of.iter().enumerate() {
            if fo as usize == f {
                bits[idx / 64] |= 1 << (idx % 64);
            }
        }
        let held = HeldOut { n_samples: n, bits };
        let d = Data {
            geno,
            held_out: Some(&held),
        };
        let refit = fit(
            &d,
            full.q.clone(),
            full.p.clone(),
            k,
            Mode::Free,
            &vec![false; n],
            options,
        );
        for j in 0..m {
            for i in 0..n {
                let g = geno.get(j, i);
                if g > 2 || !held.contains(j, i) {
                    continue;
                }
                let h = 2.0 - g as f64;
                let a: f64 = (0..k)
                    .map(|t| refit.q[i * k + t] * refit.p[j * k + t])
                    .sum();
                let mu = 2.0 * a;
                let term = |x: f64, y: f64| if x > 0.0 { x * (x / y).ln() } else { 0.0 };
                dev += term(h, mu) + term(2.0 - h, 2.0 - mu);
                count += 1;
            }
        }
    }
    Ok(dev / count as f64)
}

#[cfg(test)]
mod tests;
