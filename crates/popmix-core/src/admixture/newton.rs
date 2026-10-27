//! The Newton steps of the block relaxation: each sample's Q with P fixed,
//! then each SNP's P with Q fixed, their sums specialized on K (2-8).

use super::*;

/// Call `$f::<K>(args)` with K the compile-time value of `$k` for the
/// usual numbers of populations (so the loops over them unroll), K = 0 (the
/// value read at run time) otherwise.
macro_rules! with_k {
    ($k:expr, $f:ident($($arg:expr),*)) => {
        match $k {
            2 => $f::<2>($($arg),*),
            3 => $f::<3>($($arg),*),
            4 => $f::<4>($($arg),*),
            5 => $f::<5>($($arg),*),
            6 => $f::<6>($($arg),*),
            7 => $f::<7>($($arg),*),
            8 => $f::<8>($($arg),*),
            _ => $f::<0>($($arg),*),
        }
    };
}

/// [`q_terms`] over the SNPs `snps`, added to `out`; `K` is k, or 0.
pub(super) fn q_terms_chunk<const K: usize>(
    d: &Data,
    q: &[f64],
    p: &[f64],
    k: usize,
    snps: std::ops::Range<usize>,
    derivatives: bool,
    out: &mut (Vec<f64>, Vec<f64>, Vec<f64>),
) {
    let k = if K == 0 { k } else { K };
    let n = d.geno.n_samples();
    let (ll, grad, hess) = out;
    // each SNP's p, 1 - p and their products, shared by every sample
    let (mut np, mut pp, mut nn) = (vec![0.0; k], vec![0.0; k * k], vec![0.0; k * k]);
    for j in snps {
        let pj = &p[j * k..(j + 1) * k];
        for t in 0..k {
            np[t] = 1.0 - pj[t];
        }
        if derivatives {
            for t in 0..k {
                for u in t..k {
                    pp[t * k + u] = pj[t] * pj[u];
                    nn[t * k + u] = np[t] * np[u];
                }
            }
        }
        for i in 0..n {
            let Some(g) = d.g(j, i) else { continue };
            let qi = &q[i * k..(i + 1) * k];
            let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
            let b: f64 = qi.iter().zip(&np).map(|(x, y)| x * y).sum();
            ll[i] += term(g, a, b);
            if !derivatives {
                continue;
            }
            let gi = &mut grad[i * k..(i + 1) * k];
            let hi = &mut hess[i * k * k..(i + 1) * k * k];
            let h = 2.0 - f64::from(g);
            let (ra, rb) = (1.0 / a, 1.0 / b);
            let (wa, wb) = (h * ra, (2.0 - h) * rb);
            let (wa2, wb2) = (wa * ra, wb * rb);
            for t in 0..k {
                gi[t] += wa * pj[t] + wb * np[t];
                for u in t..k {
                    hi[t * k + u] += wa2 * pp[t * k + u] + wb2 * nn[t * k + u];
                }
            }
        }
    }
}

/// Per-sample log-likelihood contributions, and (with `derivatives`) each
/// sample's gradient and negated Hessian of its log-likelihood in Q.
pub(super) fn q_terms(
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
            let mut out = zero();
            let snps = c * chunk..((c + 1) * chunk).min(m);
            with_k!(k, q_terms_chunk(d, q, p, k, snps, derivatives, &mut out));
            out
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

/// The log-likelihood contributions of `samples` at (q, p), summed as
/// [`q_terms`] sums them, so each equals its entry there.
pub(super) fn sample_ll(d: &Data, q: &[f64], p: &[f64], k: usize, samples: &[usize]) -> Vec<f64> {
    let m = d.geno.n_snps();
    let chunk = 4 * CHUNK;
    let parts: Vec<Vec<f64>> = (0..m.div_ceil(chunk))
        .into_par_iter()
        .map(|c| {
            let mut ll = vec![0.0; samples.len()];
            for j in c * chunk..((c + 1) * chunk).min(m) {
                let pj = &p[j * k..(j + 1) * k];
                for (l, &i) in ll.iter_mut().zip(samples) {
                    let Some(g) = d.g(j, i) else { continue };
                    let qi = &q[i * k..(i + 1) * k];
                    let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                    let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                    *l += term(g, a, b);
                }
            }
            ll
        })
        .collect();
    let mut out = vec![0.0; samples.len()];
    for part in parts {
        out.iter_mut().zip(&part).for_each(|(a, b)| *a += b);
    }
    out
}

/// One Newton step on every sample's Q with P fixed, each within the
/// simplex and bounds, halved until that sample's log-likelihood does not
/// decrease.
pub(super) fn q_newton(d: &Data, q: &[f64], p: &[f64], k: usize, fixed_q: &[bool]) -> Vec<f64> {
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
    // a step predicted to gain less than the rounding error of the sample's
    // log-likelihood cannot be told from no step: keep its Q
    let (mut pending, mut kept): (Vec<usize>, Vec<usize>) =
        (0..n).filter(|&i| !fixed_q[i]).partition(|&i| {
            let r = i * k..(i + 1) * k;
            let gain = predicted_gain(
                &hess[i * k * k..(i + 1) * k * k],
                &grad[r.clone()],
                &step[r],
            );
            gain > ROUNDING * ll0[i].abs()
        });
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
        // only the samples still pending: the others are settled
        let ll1 = sample_ll(d, &out, p, k, &pending);
        let mut r = 0;
        pending.retain(|&i| {
            r += 1;
            ll1[r - 1] < ll0[i]
        });
        for &i in &pending {
            step[i * k..(i + 1) * k].iter_mut().for_each(|x| *x *= 0.5);
        }
    }
    kept.append(&mut pending);
    for &i in &kept {
        out[i * k..(i + 1) * k].copy_from_slice(&q[i * k..(i + 1) * k]);
    }
    out
}

/// SNP `j`'s log-likelihood at `pj`, its gradient in P added to `g` and its
/// negated Hessian to `a` (K x K, filled whole); None if no sample is
/// called. `qq` holds each sample's q q' (upper triangle); `K` is k, or 0.
#[allow(clippy::too_many_arguments)]
pub(super) fn p_terms<const K: usize>(
    d: &Data,
    q: &[f64],
    qq: &[f64],
    pj: &[f64],
    k: usize,
    j: usize,
    g: &mut [f64],
    a: &mut [f64],
) -> Option<f64> {
    let k = if K == 0 { k } else { K };
    let n = d.geno.n_samples();
    let np: Vec<f64> = pj.iter().map(|x| 1.0 - x).collect();
    let mut ll0 = 0.0;
    let mut any = false;
    for i in 0..n {
        let Some(gc) = d.g(j, i) else { continue };
        any = true;
        let h = 2.0 - f64::from(gc);
        let qi = &q[i * k..(i + 1) * k];
        let x: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
        let y: f64 = qi.iter().zip(&np).map(|(x, y)| x * y).sum();
        ll0 += term(gc, x, y);
        let (rx, ry) = (1.0 / x, 1.0 / y);
        let (wx, wy) = (h * rx, (2.0 - h) * ry);
        let w = wx - wy;
        let w2 = wx * rx + wy * ry;
        let qqi = &qq[i * k * k..(i + 1) * k * k];
        for t in 0..k {
            g[t] += qi[t] * w;
            for u in t..k {
                a[t * k + u] += qqi[t * k + u] * w2;
            }
        }
    }
    for t in 0..k {
        for u in 0..t {
            a[t * k + u] = a[u * k + t];
        }
    }
    any.then_some(ll0)
}

/// One Newton step on every SNP's P with Q fixed, within the bounds, halved
/// until that SNP's log-likelihood does not decrease; with the
/// log-likelihood at the new P, summed SNP by SNP in order.
pub(super) fn p_newton(d: &Data, q: &[f64], p: &[f64], k: usize) -> (Vec<f64>, f64) {
    let n = d.geno.n_samples();
    let m = d.geno.n_snps();
    let snp_ll = |j: usize, pj: &[f64]| -> f64 {
        (0..n)
            .filter_map(|i| {
                let g = d.g(j, i)?;
                let qi = &q[i * k..(i + 1) * k];
                let a: f64 = qi.iter().zip(pj).map(|(x, y)| x * y).sum();
                let b: f64 = qi.iter().zip(pj).map(|(x, y)| x * (1.0 - y)).sum();
                Some(term(g, a, b))
            })
            .sum()
    };
    // each sample's q q' (upper triangle), shared by every SNP
    let qq: Vec<f64> = (0..n)
        .flat_map(|i| {
            let qi = &q[i * k..(i + 1) * k];
            (0..k * k).map(move |tu| {
                let (t, u) = (tu / k, tu % k);
                if u >= t {
                    qi[t] * qi[u]
                } else {
                    0.0
                }
            })
        })
        .collect();
    let per_snp: Vec<(Vec<f64>, f64)> = (0..m)
        .into_par_iter()
        .map(|j| {
            let pj = &p[j * k..(j + 1) * k];
            let mut g = vec![0.0; k];
            let mut a = vec![0.0; k * k];
            let Some(ll0) = with_k!(k, p_terms(d, q, &qq, pj, k, j, &mut g, &mut a)) else {
                return (pj.to_vec(), 0.0);
            };
            // a tiny ridge keeps the system positive definite where some
            // population has no weight at this SNP
            for t in 0..k {
                a[t * k + t] += 1e-12 * (1.0 + a[t * k + t]);
            }
            let lo: Vec<f64> = pj.iter().map(|x| BOUND - x).collect();
            let hi: Vec<f64> = pj.iter().map(|x| 1.0 - BOUND - x).collect();
            let mut step = box_qp(&a, &g, &lo, &hi);
            if predicted_gain(&a, &g, &step) <= ROUNDING * ll0.abs() {
                return (pj.to_vec(), ll0);
            }
            for _ in 0..30 {
                let cand: Vec<f64> = pj.iter().zip(&step).map(|(x, s)| clamp(x + s)).collect();
                let ll = snp_ll(j, &cand);
                if ll >= ll0 {
                    return (cand, ll);
                }
                for x in &mut step {
                    *x *= 0.5;
                }
            }
            (pj.to_vec(), ll0)
        })
        .collect();
    let ll = per_snp.iter().map(|x| x.1).sum();
    (per_snp.into_iter().flat_map(|x| x.0).collect(), ll)
}
