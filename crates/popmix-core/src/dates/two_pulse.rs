//! Is one pulse enough? A check of the single-pulse model (TASKS.md H8).
//!
//! A single admixture pulse n generations ago gives an LD decay
//! A·e^(−n·d) + K; admixture over a span of time, or two pulses, give a sum
//! of such terms (Moorjani et al. 2011, Loh et al. 2013, ALDER's test for
//! multiple admixture; Chintalapati, Patterson & Moorjani 2022 for DATES).
//! popmix fits A₁·e^(−n₁·d) + A₂·e^(−n₂·d) + K to the same bins as the
//! single fit and takes the leave-one-chromosome-out jackknife of n₁, n₂,
//! A₁ and A₂.
//! Two pulses are supported only if both amplitudes are more than
//! [`Z_SUPPORT`] SEs above zero and the two dates are more than
//! [`Z_SUPPORT`] SEs of their difference apart; otherwise the curve gives no
//! evidence against one pulse (which is not evidence for one: a short span
//! of admixture fits one exponential too).
//!
//! The fit: for fixed decays the amplitudes and K are a linear
//! least-squares solution, so the error is minimised over the two decays per
//! bin, on a log-spaced grid of pairs and then by alternating golden-section
//! searches around the best pair.

use super::exponential_fit::{weighted_jackknife, DecayFitOptions};
use super::ld_decay::LdBin;
use crate::error::Result;

/// SEs a quantity must clear for two pulses to be supported.
pub const Z_SUPPORT: f64 = 2.0;
const GRID: usize = 60;
/// The dates searched, in generations: a term decaying more slowly than 1
/// generation is nearly constant over the fitted range, so the intercept
/// takes it.
const N_MIN: f64 = 1.0;
const N_MAX: f64 = 2000.0;
const ROUNDS: usize = 40;

/// A two-pulse fit: the more recent pulse first.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwoPulseFit {
    /// n₁ < n₂, in generations.
    pub generations: [f64; 2],
    /// A₁, A₂: each term's value at the first fitted bin.
    pub amplitudes: [f64; 2],
    pub intercept: f64,
    /// Mean squared error.
    pub rss: f64,
}

/// The check: the two-pulse fit with jackknife SEs, and whether it supports
/// two pulses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PulseCheck {
    pub fit: TwoPulseFit,
    pub generations_se: [f64; 2],
    pub amplitudes_se: [f64; 2],
    /// SE of n₂ − n₁.
    pub difference_se: f64,
    /// Both amplitudes and the difference of the dates more than
    /// [`Z_SUPPORT`] SEs from zero.
    pub multiple_pulses: bool,
}

/// The least-squares (A₁, A₂, K) and mean squared error for decays per bin
/// u₁, u₂, from the normal equations by Cramer's rule; +∞ error if they
/// are singular.
fn solve(y: &[f64], u: [f64; 2]) -> ([f64; 3], f64) {
    let r = [(-u[0]).exp(), (-u[1]).exp()];
    let (mut s11, mut s12, mut s22, mut s1, mut s2) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut b1, mut b2, mut b3) = (0.0, 0.0, 0.0);
    let mut x = r;
    for &yk in y {
        s11 += x[0] * x[0];
        s12 += x[0] * x[1];
        s22 += x[1] * x[1];
        s1 += x[0];
        s2 += x[1];
        b1 += x[0] * yk;
        b2 += x[1] * yk;
        b3 += yk;
        x = [x[0] * r[0], x[1] * r[1]];
    }
    let m = y.len() as f64;
    let det3 = |a: [[f64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let n = [[s11, s12, s1], [s12, s22, s2], [s1, s2, m]];
    let det = det3(n);
    if det.abs() <= 1e-14 * (s11 * s22 * m).abs() || det == 0.0 {
        return ([f64::NAN; 3], f64::INFINITY);
    }
    let rhs = [b1, b2, b3];
    let mut beta = [0.0; 3];
    for (c, out) in beta.iter_mut().enumerate() {
        let mut a = n;
        for (row, &v) in a.iter_mut().zip(&rhs) {
            row[c] = v;
        }
        *out = det3(a) / det;
    }
    let mut x = r;
    let mut sse = 0.0;
    for &yk in y {
        let e = yk - (beta[0] * x[0] + beta[1] * x[1] + beta[2]);
        sse += e * e;
        x = [x[0] * r[0], x[1] * r[1]];
    }
    let mse = sse / m;
    (beta, if mse.is_finite() { mse } else { f64::INFINITY })
}

/// Golden-section minimum of f on [lo, hi].
fn golden(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let g = (5f64.sqrt() - 1.0) / 2.0;
    let (mut c, mut d) = (hi - g * (hi - lo), lo + g * (hi - lo));
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..100 {
        if hi - lo <= 1e-12 * hi {
            break;
        }
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - g * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + g * (hi - lo);
            fd = f(d);
        }
    }
    if fc <= fd {
        c
    } else {
        d
    }
}

/// The two-pulse fit of `y`, the correlations of evenly spaced bins `step`
/// Morgans apart (decays per bin converted to generations), or None with
/// fewer than 6 bins or a singular fit.
pub fn fit_two_pulse_values(y: &[f64], step: f64) -> Option<TwoPulseFit> {
    if y.len() < 6 {
        return None;
    }
    let (u_min, u_max) = (N_MIN * step, N_MAX * step);
    let ratio = (u_max / u_min).powf(1.0 / (GRID - 1) as f64);
    let grid: Vec<f64> = (0..GRID).map(|t| u_min * ratio.powi(t as i32)).collect();
    let mut best = ([f64::NAN; 2], f64::INFINITY);
    for i in 0..GRID {
        for j in i + 2..GRID {
            let e = solve(y, [grid[i], grid[j]]).1;
            if e < best.1 {
                best = ([grid[i], grid[j]], e);
            }
        }
    }
    let mut u = best.0;
    if !u[0].is_finite() {
        return None;
    }
    for _ in 0..ROUNDS {
        let before = u;
        // u[0] < u[1] throughout, kept a grid step apart so that the two
        // terms never merge into one with cancelling amplitudes
        for k in 0..2 {
            let other = u[1 - k];
            let (lo, hi) = if k == 0 {
                ((u[0] / ratio).max(u_min), (u[0] * ratio).min(other / ratio))
            } else {
                ((u[1] / ratio).max(other * ratio), (u[1] * ratio).min(u_max))
            };
            if lo < hi {
                u[k] = golden(
                    |v| {
                        let mut w = [other; 2];
                        w[k] = v;
                        solve(y, w).1
                    },
                    lo,
                    hi,
                );
            }
        }
        if (u[0] - before[0]).abs() <= 1e-12 * u[0] && (u[1] - before[1]).abs() <= 1e-12 * u[1] {
            break;
        }
    }
    let (beta, mse) = solve(y, u);
    if !mse.is_finite() {
        return None;
    }
    let r = [(-u[0]).exp(), (-u[1]).exp()];
    Some(TwoPulseFit {
        generations: [u[0] / step, u[1] / step],
        amplitudes: [beta[0] * r[0], beta[1] * r[1]],
        intercept: beta[2],
        rss: mse,
    })
}

/// The two-pulse fit of the bins in `options`' fitting range, or None if
/// there are too few or the fit is singular.
pub fn fit_two_pulse(ld_bins: &[LdBin], options: &DecayFitOptions) -> Option<TwoPulseFit> {
    let y: Vec<f64> = ld_bins
        .iter()
        .skip_while(|b| b.distance_cm < options.min_distance_cm)
        .take_while(|b| b.distance_cm <= options.max_distance_cm)
        .map(|b| b.correlation)
        .collect();
    if y.iter().any(|v| !v.is_finite()) {
        return None;
    }
    fit_two_pulse_values(&y, options.bin_size_cm / 100.0)
}

/// The check from the full fit and the leave-one-block-out fits with their
/// block weights.
pub fn check_pulses(full: &TwoPulseFit, left_out: &[(f64, TwoPulseFit)]) -> Result<PulseCheck> {
    let jk = |f: &dyn Fn(&TwoPulseFit) -> f64| -> Result<f64> {
        let estimates: Vec<(f64, f64)> = left_out.iter().map(|(w, t)| (*w, f(t))).collect();
        Ok(weighted_jackknife(f(full), &estimates)?.1)
    };
    let generations_se = [jk(&|t| t.generations[0])?, jk(&|t| t.generations[1])?];
    let amplitudes_se = [jk(&|t| t.amplitudes[0])?, jk(&|t| t.amplitudes[1])?];
    let difference_se = jk(&|t| t.generations[1] - t.generations[0])?;
    let above = |x: f64, se: f64| x > Z_SUPPORT * se;
    let multiple_pulses = above(full.amplitudes[0], amplitudes_se[0])
        && above(full.amplitudes[1], amplitudes_se[1])
        && above(full.generations[1] - full.generations[0], difference_se);
    Ok(PulseCheck {
        fit: *full,
        generations_se,
        amplitudes_se,
        difference_se,
        multiple_pulses,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 0.2 cM bins, few, so that the tests run fast unoptimised.
    const STEP: f64 = 0.002;

    /// Correlations of 160 bins from 0.5 cM: Σ aᵢ e^(−nᵢ d) + k, with a
    /// deterministic wiggle of size `noise` and phase `seed`.
    fn curve(terms: &[(f64, f64)], k: f64, noise: f64, seed: f64) -> Vec<f64> {
        (0..160)
            .map(|i: i32| {
                let i = f64::from(i);
                let d = 0.005 + i * STEP;
                let wiggle = noise * ((i * 1.7 + seed).sin() + (i * 0.31 + 2.0 * seed).cos());
                terms.iter().map(|(a, n)| a * (-n * d).exp()).sum::<f64>() + k + wiggle
            })
            .collect()
    }

    #[test]
    fn two_pulses_are_recovered() {
        let y = curve(&[(0.01, 8.0), (0.006, 80.0)], 1e-4, 0.0, 0.0);
        let f = fit_two_pulse_values(&y, STEP).unwrap();
        assert!((f.generations[0] - 8.0).abs() < 1e-3, "{f:?}");
        assert!((f.generations[1] - 80.0).abs() < 1e-2, "{f:?}");
    }

    fn check(terms: &[(f64, f64)]) -> PulseCheck {
        let full = fit_two_pulse_values(&curve(terms, 1e-4, 2e-5, 0.0), STEP).unwrap();
        let left_out: Vec<(f64, TwoPulseFit)> = (1..=8)
            .map(|b| {
                (
                    1.0,
                    fit_two_pulse_values(&curve(terms, 1e-4, 2e-5, f64::from(b)), STEP).unwrap(),
                )
            })
            .collect();
        check_pulses(&full, &left_out).unwrap()
    }

    #[test]
    fn two_distinct_pulses_are_supported() {
        assert!(check(&[(0.01, 8.0), (0.006, 80.0)]).multiple_pulses);
    }

    #[test]
    fn one_pulse_is_not_taken_for_two() {
        assert!(!check(&[(0.012, 20.0)]).multiple_pulses);
    }
}
