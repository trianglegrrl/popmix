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

use rayon::prelude::*;

use nalgebra::{DMatrix, DVector};

use crate::error::{PopmixError, Result};
use crate::pca::PackedGenotypes;

/// ADMIXTURE's bounds on Q and P.
pub const BOUND: f64 = 1e-5;

/// Options for fitting.
#[derive(Debug, Clone)]
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
    /// h = 2 - g, the count of the allele P refers to; None if not used.
    #[inline]
    fn h(&self, snp: usize, sample: usize) -> Option<f64> {
        let g = self.geno.get(snp, sample);
        if g > 2 || self.held_out.is_some_and(|m| m.contains(snp, sample)) {
            None
        } else {
            Some(2.0 - g as f64)
        }
    }
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
                    if let Some(h) = d.h(j, i) {
                        let qi = &q[i * k..(i + 1) * k];
                        let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                        let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                        l += h * a.ln() + (2.0 - h) * b.ln();
                    }
                }
            }
            l
        })
        .collect();
    parts.iter().sum()
}

/// Minimize 1/2 x'Ax - g'x over lo <= x <= hi (lo <= 0 <= hi), A positive
/// definite (K x K, row-major): an active-set method over the bounds.
fn box_qp(a: &[f64], g: &[f64], lo: &[f64], hi: &[f64]) -> Vec<f64> {
    let k = g.len();
    let mut x = vec![0.0; k];
    // -1 at the lower bound, 1 at the upper, 0 free
    let mut active = vec![0i8; k];
    for _ in 0..4 * k + 4 {
        let free: Vec<usize> = (0..k).filter(|&t| active[t] == 0).collect();
        let nf = free.len();
        let mut target = x.clone();
        for t in 0..k {
            match active[t] {
                -1 => target[t] = lo[t],
                1 => target[t] = hi[t],
                _ => {}
            }
        }
        if nf > 0 {
            let m = DMatrix::from_fn(nf, nf, |r, c| a[free[r] * k + free[c]]);
            let rhs = DVector::from_fn(nf, |r, _| {
                let t = free[r];
                g[t] - (0..k)
                    .filter(|&u| active[u] != 0)
                    .map(|u| a[t * k + u] * target[u])
                    .sum::<f64>()
            });
            let Some(sol) = m.cholesky().map(|c| c.solve(&rhs)) else {
                break;
            };
            for (r, &t) in free.iter().enumerate() {
                target[t] = sol[r];
            }
        }
        // step towards the target as far as the bounds allow
        let mut alpha = 1.0f64;
        let mut blocking = None;
        for &t in &free {
            let (bound, side) = if target[t] < lo[t] {
                (lo[t], -1)
            } else if target[t] > hi[t] {
                (hi[t], 1)
            } else {
                continue;
            };
            let step = (bound - x[t]) / (target[t] - x[t]);
            if step < alpha {
                alpha = step.max(0.0);
                blocking = Some((t, side));
            }
        }
        for t in 0..k {
            x[t] += alpha * (target[t] - x[t]);
        }
        if let Some((t, side)) = blocking {
            x[t] = if side < 0 { lo[t] } else { hi[t] };
            active[t] = side;
            continue;
        }
        // release a bound the minimum lies away from
        let release = (0..k).filter(|&t| active[t] != 0).find(|&t| {
            let grad: f64 = (0..k).map(|u| a[t * k + u] * x[u]).sum::<f64>() - g[t];
            (active[t] < 0 && grad < -1e-12) || (active[t] > 0 && grad > 1e-12)
        });
        match release {
            Some(t) => active[t] = 0,
            None => break,
        }
    }
    x
}

/// Minimize 1/2 x'Ax - g'x subject to sum(x) = 0 and x >= lo (lo <= 0), A
/// positive definite: an active-set method over the lower bounds.
fn simplex_qp(a: &[f64], g: &[f64], lo: &[f64]) -> Vec<f64> {
    let k = g.len();
    let mut x = vec![0.0; k];
    let mut active = vec![false; k];
    for _ in 0..4 * k + 4 {
        // equality-constrained minimizer over the free coordinates, the
        // active ones held at their bound
        let free: Vec<usize> = (0..k).filter(|&t| !active[t]).collect();
        let nf = free.len();
        if nf == 0 {
            break;
        }
        let fixed_sum: f64 = (0..k).filter(|&t| active[t]).map(|t| lo[t]).sum();
        let mut m = DMatrix::<f64>::zeros(nf + 1, nf + 1);
        let mut rhs = DVector::<f64>::zeros(nf + 1);
        for (r, &t) in free.iter().enumerate() {
            for (c, &u) in free.iter().enumerate() {
                m[(r, c)] = a[t * k + u];
            }
            m[(r, nf)] = 1.0;
            m[(nf, r)] = 1.0;
            let bound: f64 = (0..k)
                .filter(|&u| active[u])
                .map(|u| a[t * k + u] * lo[u])
                .sum();
            rhs[r] = g[t] - bound;
        }
        rhs[nf] = -fixed_sum;
        let Some(sol) = m.lu().solve(&rhs) else { break };
        let mut target = x.clone();
        for t in 0..k {
            if active[t] {
                target[t] = lo[t];
            }
        }
        for (r, &t) in free.iter().enumerate() {
            target[t] = sol[r];
        }
        // step towards the target as far as the bounds allow
        let mut alpha = 1.0f64;
        let mut blocking = None;
        for &t in &free {
            if target[t] < lo[t] {
                let step = (lo[t] - x[t]) / (target[t] - x[t]);
                if step < alpha {
                    alpha = step.max(0.0);
                    blocking = Some(t);
                }
            }
        }
        for t in 0..k {
            x[t] += alpha * (target[t] - x[t]);
        }
        if let Some(t) = blocking {
            x[t] = lo[t];
            active[t] = true;
            continue;
        }
        // release an active bound whose multiplier says the minimum is inside
        let lambda = sol[nf];
        let release = (0..k).filter(|&t| active[t]).find(|&t| {
            let grad: f64 = (0..k).map(|u| a[t * k + u] * x[u]).sum::<f64>() - g[t];
            grad + lambda < -1e-12
        });
        match release {
            Some(t) => active[t] = false,
            None => break,
        }
    }
    x
}

/// Per-sample log-likelihood contributions, and (with `derivatives`) each
/// sample's gradient and negated Hessian of its log-likelihood in Q.
fn q_terms(
    d: &Data,
    q: &[f64],
    p: &[f64],
    k: usize,
    derivatives: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = d.geno.n_samples();
    let m = d.geno.n_snps();
    // larger chunks: each holds a K x K matrix per sample
    let chunk = 4 * CHUNK;
    let zero = || {
        (
            vec![0.0; n],
            vec![0.0; if derivatives { n * k } else { 0 }],
            vec![0.0; if derivatives { n * k * k } else { 0 }],
        )
    };
    let parts: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..m.div_ceil(chunk))
        .into_par_iter()
        .map(|c| {
            let (mut ll, mut grad, mut hess) = zero();
            for j in c * chunk..((c + 1) * chunk).min(m) {
                let pj = &p[j * k..(j + 1) * k];
                for i in 0..n {
                    let Some(h) = d.h(j, i) else { continue };
                    let qi = &q[i * k..(i + 1) * k];
                    let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                    let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                    ll[i] += h * a.ln() + (2.0 - h) * b.ln();
                    if derivatives {
                        let (wa, wb) = (h / a, (2.0 - h) / b);
                        let (wa2, wb2) = (wa / a, wb / b);
                        for t in 0..k {
                            grad[i * k + t] += wa * pj[t] + wb * (1.0 - pj[t]);
                            for u in t..k {
                                hess[(i * k + t) * k + u] +=
                                    wa2 * pj[t] * pj[u] + wb2 * (1.0 - pj[t]) * (1.0 - pj[u]);
                            }
                        }
                    }
                }
            }
            (ll, grad, hess)
        })
        .collect();
    let mut out = zero();
    for y in parts {
        out.0.iter_mut().zip(&y.0).for_each(|(a, b)| *a += b);
        out.1.iter_mut().zip(&y.1).for_each(|(a, b)| *a += b);
        out.2.iter_mut().zip(&y.2).for_each(|(a, b)| *a += b);
    }
    out
}

/// One Newton step on every sample's Q with P fixed, each within the
/// simplex and bounds, halved until that sample's log-likelihood does not
/// decrease.
fn q_newton(d: &Data, q: &[f64], p: &[f64], k: usize, fixed_q: &[bool]) -> Vec<f64> {
    let n = d.geno.n_samples();
    let (ll0, grad, mut hess) = q_terms(d, q, p, k, true);
    for i in 0..n {
        for t in 0..k {
            for u in 0..t {
                hess[(i * k + t) * k + u] = hess[(i * k + u) * k + t];
            }
        }
    }
    let mut step: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            if fixed_q[i] {
                return vec![0.0; k];
            }
            let qi = &q[i * k..(i + 1) * k];
            let a = &hess[i * k * k..(i + 1) * k * k];
            let g = &grad[i * k..(i + 1) * k];
            let lo: Vec<f64> = qi.iter().map(|x| (BOUND - x).min(0.0)).collect();
            simplex_qp(a, g, &lo)
        })
        .collect();
    let mut out = q.to_vec();
    let mut pending: Vec<usize> = (0..n).filter(|&i| !fixed_q[i]).collect();
    for _ in 0..30 {
        if pending.is_empty() {
            break;
        }
        for &i in &pending {
            for t in 0..k {
                out[i * k + t] = q[i * k + t] + step[i * k + t];
            }
        }
        normalize_q(&mut out, k);
        let (ll1, _, _) = q_terms(d, &out, p, k, false);
        pending.retain(|&i| ll1[i] < ll0[i]);
        for &i in &pending {
            step[i * k..(i + 1) * k].iter_mut().for_each(|x| *x *= 0.5);
        }
    }
    for &i in &pending {
        out[i * k..(i + 1) * k].copy_from_slice(&q[i * k..(i + 1) * k]);
    }
    out
}

/// One Newton step on every SNP's P with Q fixed, within the bounds, halved
/// until that SNP's log-likelihood does not decrease.
fn p_newton(d: &Data, q: &[f64], p: &[f64], k: usize) -> Vec<f64> {
    let n = d.geno.n_samples();
    let m = d.geno.n_snps();
    let snp_ll = |j: usize, pj: &[f64]| -> f64 {
        (0..n)
            .filter_map(|i| {
                let h = d.h(j, i)?;
                let qi = &q[i * k..(i + 1) * k];
                let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                Some(h * a.ln() + (2.0 - h) * b.ln())
            })
            .sum()
    };
    (0..m)
        .into_par_iter()
        .flat_map_iter(|j| {
            let pj = &p[j * k..(j + 1) * k];
            let mut g = vec![0.0; k];
            let mut a = vec![0.0; k * k];
            let mut ll0 = 0.0;
            let mut any = false;
            for i in 0..n {
                let Some(h) = d.h(j, i) else { continue };
                any = true;
                let qi = &q[i * k..(i + 1) * k];
                let x: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                let y: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                ll0 += h * x.ln() + (2.0 - h) * y.ln();
                let w = h / x - (2.0 - h) / y;
                let w2 = h / (x * x) + (2.0 - h) / (y * y);
                for t in 0..k {
                    g[t] += qi[t] * w;
                    for u in 0..k {
                        a[t * k + u] += qi[t] * qi[u] * w2;
                    }
                }
            }
            if !any {
                return pj.to_vec();
            }
            // a tiny ridge keeps the system positive definite where some
            // population has no weight at this SNP
            for t in 0..k {
                a[t * k + t] += 1e-12 * (1.0 + a[t * k + t]);
            }
            let lo: Vec<f64> = pj.iter().map(|x| BOUND - x).collect();
            let hi: Vec<f64> = pj.iter().map(|x| 1.0 - BOUND - x).collect();
            let mut step = box_qp(&a, &g, &lo, &hi);
            for _ in 0..30 {
                let cand: Vec<f64> = pj.iter().zip(&step).map(|(x, s)| clamp(x + s)).collect();
                if snp_ll(j, &cand) >= ll0 {
                    return cand;
                }
                step.iter_mut().for_each(|x| *x *= 0.5);
            }
            pj.to_vec()
        })
        .collect()
}

/// One block-relaxation iteration (ADMIXTURE's method): a Newton step on Q
/// with P fixed, then on P with the new Q fixed.
fn block_step(
    d: &Data,
    q: &[f64],
    p: &[f64],
    k: usize,
    mode: Mode,
    fixed_q: &[bool],
) -> (Vec<f64>, Vec<f64>) {
    let q1 = q_newton(d, q, p, k, fixed_q);
    let p1 = if mode == Mode::FixedP {
        p.to_vec()
    } else {
        p_newton(d, &q1, p, k)
    };
    (q1, p1)
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
        let (q1, p1) = block_step(d, &q, &p, k, mode, fixed_q);
        let (q2, p2) = block_step(d, &q1, &p1, k, mode, fixed_q);
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
            let (qs, ps) = block_step(d, &qa, &pa, k, mode, fixed_q);
            let ls = loglik(d, &qs, &ps, k);
            candidate = Some((qs, ps, ls));
        }
        let l2 = loglik(d, &q2, &p2, k);
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
        if gain.abs() < options.tolerance.max(1e-13 * ll.abs()) {
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
        return Err(PopmixError::Numerical("K must be at least 1".to_string()));
    }
    if geno.n_samples() < 2 || geno.n_snps() < 1 {
        return Err(PopmixError::Numerical(format!(
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
        return Err(PopmixError::Numerical(format!(
            "label {bad} out of range for K = {k}"
        )));
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
        return Err(PopmixError::Validation(format!(
            "cross-validation takes 2 to 256 folds, not {folds}"
        )));
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
mod tests {
    use super::*;
    use crate::fstats::MatrixRows;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    /// Two source populations with different frequencies, and samples that
    /// are mixtures of them with known proportions.
    fn simulated(qs: &[f64], m: usize) -> (PackedGenotypes, Vec<f64>) {
        let mut rng = Rng(7);
        let p: Vec<[f64; 2]> = (0..m)
            .map(|_| [0.05 + 0.9 * rng.uniform(), 0.05 + 0.9 * rng.uniform()])
            .collect();
        let n = qs.len();
        let data = Array2::from_shape_fn((m, n), |(j, i)| {
            let f = qs[i] * p[j][0] + (1.0 - qs[i]) * p[j][1];
            // count of the non-P allele
            (0..2).filter(|_| rng.uniform() >= f).count() as u8
        });
        let samples: Vec<Sample> = (0..n).map(|i| Sample::new(format!("s{i}"), "P")).collect();
        let snps: Vec<Snp> = (0..m)
            .map(|j| {
                Snp::new(
                    format!("rs{j}"),
                    Chromosome::Autosome(1),
                    j as u32,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let mat = GenotypeMatrix::new(data, samples, snps).unwrap();
        let mut rows = MatrixRows::new(&mat, (0..m).collect(), (0..n).collect());
        let flat: Vec<f64> = p.iter().flat_map(|x| x.iter().copied()).collect();
        (PackedGenotypes::read(&mut rows, n, m).unwrap(), flat)
    }

    #[test]
    fn unsupervised_recovers_the_proportions() {
        let qs = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.25, 0.75, 1.0, 0.0, 0.5];
        let (g, _) = simulated(&qs, 4000);
        let opts = AdmixtureOptions {
            tolerance: 1e-8,
            ..Default::default()
        };
        let fit = admixture(&g, 2, &opts).unwrap();
        assert!(fit.converged);
        // cluster labels are arbitrary
        let col = if fit.q[0] > 0.5 { 0 } else { 1 };
        for (i, &want) in qs.iter().enumerate() {
            assert!(
                (fit.q[i * 2 + col] - want).abs() < 0.05,
                "sample {i}: {} vs {want}",
                fit.q[i * 2 + col]
            );
        }
    }

    #[test]
    fn the_loglikelihood_never_decreases_and_matches_its_definition() {
        let qs = [1.0, 0.0, 0.5, 0.3, 0.9, 0.1];
        let (g, _) = simulated(&qs, 500);
        let fit = admixture(&g, 2, &AdmixtureOptions::default()).unwrap();
        let again = admixture_loglik(&g, &fit.q, &fit.p, 2).unwrap();
        assert!((again - fit.loglik).abs() < 1e-9 * fit.loglik.abs());
        let (q0, p0) = random_start(6, 500, 2, 1);
        assert!(admixture_loglik(&g, &q0, &p0, 2).unwrap() < fit.loglik);
    }

    #[test]
    fn projection_with_the_true_frequencies_recovers_the_proportions() {
        let qs = [0.2, 0.8, 0.5, 1.0];
        let (g, p) = simulated(&qs, 20000);
        let opts = AdmixtureOptions {
            tolerance: 1e-10,
            ..Default::default()
        };
        let fit = admixture_projection(&g, &p, 2, &opts).unwrap();
        for (i, &want) in qs.iter().enumerate() {
            assert!(
                (fit.q[i * 2] - want).abs() < 0.02,
                "sample {i}: {} vs {want}",
                fit.q[i * 2]
            );
        }
        assert_eq!(fit.p.len(), p.len());
    }

    #[test]
    fn supervised_samples_keep_their_population() {
        let qs = [1.0, 1.0, 0.0, 0.0, 0.4];
        let (g, _) = simulated(&qs, 2000);
        let labels = [Some(0), Some(0), Some(1), Some(1), None];
        let opts = AdmixtureOptions {
            tolerance: 1e-8,
            ..Default::default()
        };
        let fit = admixture_supervised(&g, 2, &labels, &opts).unwrap();
        // sample 0 is in population 0, sample 2 in population 1
        assert!(fit.q[0] > 0.9999 && fit.q[4] < 1e-4 && fit.q[5] > 0.9999);
        assert!((fit.q[8] - 0.4).abs() < 0.05, "{}", fit.q[8]);
    }

    #[test]
    fn cross_validation_prefers_the_true_k() {
        let qs: Vec<f64> = (0..20).map(|i| if i < 10 { 1.0 } else { 0.0 }).collect();
        let (g, _) = simulated(&qs, 1500);
        let opts = AdmixtureOptions {
            tolerance: 1e-6,
            ..Default::default()
        };
        let cv = |k| admixture_cv(&g, &admixture(&g, k, &opts).unwrap(), 5, &opts).unwrap();
        assert!(cv(2) < cv(1));
    }

    #[test]
    fn cross_validation_folds_are_checked() {
        let qs: Vec<f64> = (0..4).map(|i| if i < 2 { 1.0 } else { 0.0 }).collect();
        let (g, _) = simulated(&qs, 50);
        let opts = AdmixtureOptions::default();
        let fit = admixture(&g, 1, &opts).unwrap();
        for folds in [1, 257, 1000] {
            let err = admixture_cv(&g, &fit, folds, &opts)
                .unwrap_err()
                .to_string();
            assert!(err.contains("folds"), "{folds}: {err}");
        }
    }

    #[test]
    fn box_qp_solves_inside_and_on_the_bounds() {
        let a = [2.0, 0.5, 0.5, 1.0];
        // unconstrained minimum A^-1 g
        let x = box_qp(&a, &[1.0, 1.0], &[-10.0, -10.0], &[10.0, 10.0]);
        assert!(
            (2.0 * x[0] + 0.5 * x[1] - 1.0).abs() < 1e-12
                && (0.5 * x[0] + x[1] - 1.0).abs() < 1e-12
        );
        let x = box_qp(&a, &[1.0, 1.0], &[-10.0, -10.0], &[10.0, 0.1]);
        assert_eq!(x[1], 0.1);
        assert!((x[0] - (1.0 - 0.05) / 2.0).abs() < 1e-12);
    }

    #[test]
    fn simplex_qp_keeps_the_sum_and_the_bounds() {
        let a = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        // pull hard towards the first coordinate; the third may only drop by 0.1
        let x = simplex_qp(&a, &[3.0, 0.0, -3.0], &[-0.5, -0.5, -0.1]);
        assert!(x.iter().sum::<f64>().abs() < 1e-12);
        // x3 = -0.1 and then x2 = -0.5 hit their bounds; x1 takes the rest
        for (got, want) in x.iter().zip([0.6, -0.5, -0.1]) {
            assert!((got - want).abs() < 1e-12, "{x:?}");
        }
        let free = simplex_qp(&a, &[1.0, 0.0, -1.0], &[-5.0, -5.0, -5.0]);
        assert!((free[0] - 1.0).abs() < 1e-12 && free[1].abs() < 1e-12);
    }

    #[test]
    fn bad_inputs_are_errors() {
        let (g, p) = simulated(&[0.5, 0.5], 10);
        assert!(admixture(&g, 0, &AdmixtureOptions::default()).is_err());
        assert!(admixture_projection(&g, &p[..4], 2, &AdmixtureOptions::default()).is_err());
        assert!(
            admixture_supervised(&g, 2, &[Some(2), None], &AdmixtureOptions::default()).is_err()
        );
    }
}
