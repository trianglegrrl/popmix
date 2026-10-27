//! The two small quadratic programs of each Newton step: over a box (a
//! SNP's frequencies) and over the simplex (a sample's proportions).

use super::*;

/// The increase of the quadratic model g's - 1/2 s'As (A the negated
/// Hessian, K x K row-major) for the step s.
pub(super) fn predicted_gain(a: &[f64], g: &[f64], s: &[f64]) -> f64 {
    let k = g.len();
    let mut quad = 0.0;
    for t in 0..k {
        for u in 0..k {
            quad += s[t] * a[t * k + u] * s[u];
        }
    }
    g.iter().zip(s).map(|(x, y)| x * y).sum::<f64>() - 0.5 * quad
}

/// Minimize 1/2 x'Ax - g'x over lo <= x <= hi (lo <= 0 <= hi), A positive
/// definite (K x K, row-major): an active-set method over the bounds.
pub(super) fn box_qp(a: &[f64], g: &[f64], lo: &[f64], hi: &[f64]) -> Vec<f64> {
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
pub(super) fn simplex_qp(a: &[f64], g: &[f64], lo: &[f64]) -> Vec<f64> {
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
