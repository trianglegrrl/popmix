//! Polishing the best admixture weights past a zero-bound kink
//! (TASKS.md G14, DEVIATIONS.md ADMIXTOOLS2-0009).
//!
//! Where a fitted drift edge reaches its bound of 0 the score has a kink.
//! L-BFGS-B with finite-difference gradients, as R's `optim` (and so
//! admixtools2), straddles it and stops short from every start. Starting
//! from the best start's weights, each weight in turn is moved by a
//! golden-section search over a window around it (clipped to its bounds),
//! keeping a move only if it lowers the score, until a round improves the
//! score by less than a relative 1e-12. Near a kink the score is unimodal
//! in each weight, which golden-section search needs.
//!
//! Provenance: golden-section search (Kiefer 1953); written here.

use crate::Result;

/// Half-width of the window searched around each weight.
const WINDOW: f64 = 0.01;
/// Golden-section stops when the bracket is this narrow.
const TOLERANCE: f64 = 1e-12;
/// Rounds over the weights at most.
const MAX_ROUNDS: usize = 50;

/// The weights `w` (within `lo`..`hi`) moved, one at a time, to lower
/// `score`, and the score there.
pub(crate) fn polish(
    score: &mut dyn FnMut(&[f64]) -> Result<f64>,
    w: &[f64],
    lo: &[f64],
    hi: &[f64],
) -> Result<(Vec<f64>, f64)> {
    let mut w = w.to_vec();
    let mut best = score(&w)?;
    for _ in 0..MAX_ROUNDS {
        let before = best;
        for j in 0..w.len() {
            let (a, b) = ((w[j] - WINDOW).max(lo[j]), (w[j] + WINDOW).min(hi[j]));
            let mut at = |x: f64| -> Result<f64> {
                let mut v = w.clone();
                v[j] = x;
                score(&v)
            };
            let (x, s) = golden_section(&mut at, a, b)?;
            if s < best {
                w[j] = x;
                best = s;
            }
        }
        if before - best <= TOLERANCE * before.abs() {
            break;
        }
    }
    Ok((w, best))
}

/// The minimum of a unimodal `f` on [a, b] by golden-section search.
fn golden_section(
    f: &mut dyn FnMut(f64) -> Result<f64>,
    mut a: f64,
    mut b: f64,
) -> Result<(f64, f64)> {
    let r = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = b - r * (b - a);
    let mut d = a + r * (b - a);
    let (mut fc, mut fd) = (f(c)?, f(d)?);
    while (b - a).abs() > TOLERANCE {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - r * (b - a);
            fc = f(c)?;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + r * (b - a);
            fd = f(d)?;
        }
    }
    let x = f64::midpoint(a, b);
    Ok((x, f(x)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polish_reaches_the_minimum_at_a_kink() {
        // a kink at 0.45087, approached from the stall point 0.4502
        let mut f = |w: &[f64]| -> Result<f64> {
            Ok(0.677059 + (w[0] - 0.45087).abs() + 0.1 * (w[1] - 0.3).powi(2))
        };
        let (w, s) = polish(&mut f, &[0.4502, 0.31], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
        assert!((w[0] - 0.45087).abs() < 1e-9, "{w:?}");
        assert!((w[1] - 0.3).abs() < 1e-5, "{w:?}");
        assert!((s - 0.677059).abs() < 1e-9);
    }

    #[test]
    fn polish_stays_within_the_bounds() {
        let mut f = |w: &[f64]| -> Result<f64> { Ok(-w[0]) };
        let (w, _) = polish(&mut f, &[0.995], &[0.0], &[1.0]).unwrap();
        assert!(w[0] <= 1.0);
        assert!(w[0] > 0.999_999);
    }

    #[test]
    fn polish_never_makes_the_score_worse() {
        let mut f = |w: &[f64]| -> Result<f64> { Ok((w[0] - 0.2).powi(2)) };
        let (w, s) = polish(&mut f, &[0.2], &[0.0], &[1.0]).unwrap();
        assert_eq!(w, vec![0.2]);
        assert_eq!(s, 0.0);
    }
}
