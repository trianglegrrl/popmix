//! The low-rank fit behind qpWave, qpAdm, `qpadm_p`, rotation and pairs.
//!
//! Written from the clean-room specification `docs/cleanroom/qpadm-fit.md`
//! and the published method: Patterson et al. 2012 (Genetics, appendix) for
//! the rank test of an f4 matrix, and Haak et al. 2015 (Nature, SI 10) for
//! the qpAdm weights.
//!
//! The m x n f4 matrix X (rows: left populations L1..Lm against L0; columns:
//! right populations R1..Rn against R0) is approximated by a product A B of
//! rank r, in the metric of the inverse jackknife covariance W of vec(X).
//! Every vectorisation here is row-major: entry (i, j) of an m x n matrix is
//! element i·n + j of its vector.
//!
//! The fit alternates between two weighted least-squares problems (A with B
//! fixed, then B with A fixed), each a small linear system with a ridge of
//! `fudge` times its trace. The misfit chi² = (x − vec(AB))ᵀ W (x − vec(AB))
//! is compared with a chi-squared distribution on (m − r)(n − r) degrees of
//! freedom. For qpAdm (target L0, sources L1..Lm, rank m − 1) the weights are
//! the combination of the sources' rows of A that vanishes, normalised to sum
//! to one.

use std::cmp::Ordering;

use nalgebra::{DMatrix, DVector};
use ndarray::Array2;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use super::compat::Admixtools2Compat;
use crate::{PopmixError, Result};

/// Number of alternating rounds (A-step then B-step) in [`fit_rank`].
const FIT_ROUNDS: usize = 20;

/// Convert an ndarray matrix into a nalgebra one (same shape and entries).
pub(crate) fn to_dmatrix(a: &Array2<f64>) -> DMatrix<f64> {
    let (nr, nc) = a.dim();
    DMatrix::from_fn(nr, nc, |i, j| a[[i, j]])
}

/// Whether an error means "this model cannot be fitted" (as opposed to bad
/// input): a singular system or an infeasible weight problem. Callers that
/// scan many models skip these.
pub(crate) fn is_model_failure(e: &PopmixError) -> bool {
    matches!(
        e,
        PopmixError::SingularMatrix { .. } | PopmixError::Infeasible { .. }
    )
}

fn singular(context: &str) -> PopmixError {
    PopmixError::SingularMatrix {
        context: context.to_string(),
    }
}

/// Add `fudge · tr(M)` to every diagonal element of a square matrix.
fn add_trace_ridge(m: &mut DMatrix<f64>, fudge: f64) {
    let ridge = fudge * m.trace();
    for i in 0..m.nrows() {
        m[(i, i)] += ridge;
    }
}

/// The weight matrix W = (Q + f·tr(Q)·I)⁻¹ of the f4 entries (spec section 1).
pub(crate) fn fudged_inverse(cov: &Array2<f64>, fudge: f64) -> Result<DMatrix<f64>> {
    let mut q = to_dmatrix(cov);
    add_trace_ridge(&mut q, fudge);
    q.try_inverse()
        .filter(|inv| inv.iter().all(|v| v.is_finite()))
        .ok_or_else(|| singular("f4 covariance matrix"))
}

/// A rank-r approximation X ≈ A B: `a` is m x r and `b` is r x n.
#[derive(Debug, Clone)]
pub(crate) struct RankFit {
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
}

/// Row-major vectorisation of a matrix.
fn vec_row_major(x: &DMatrix<f64>) -> DVector<f64> {
    let (nr, nc) = x.shape();
    DVector::from_fn(nr * nc, |k, _| x[(k / nc, k % nc)])
}

/// Solve the ridged normal equations (DᵀWD + f·tr·I) θ = DᵀWx for the
/// weighted least-squares problem min (x − Dθ)ᵀ W (x − Dθ).
fn ridged_wls(
    design: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    x: &DVector<f64>,
    fudge: f64,
) -> Result<DVector<f64>> {
    let dt_w = design.transpose() * qinv;
    let mut normal = &dt_w * design;
    add_trace_ridge(&mut normal, fudge);
    let rhs = dt_w * x;
    normal
        .lu()
        .solve(&rhs)
        .filter(|s| s.iter().all(|v| v.is_finite()))
        .ok_or_else(|| singular("qpWave/qpAdm fit"))
}

/// Fit X ≈ A B at rank `rnk` in the metric `qinv` (spec section 2).
///
/// Starts from the leading `rnk` right singular vectors of X, then runs
/// [`FIT_ROUNDS`] rounds of alternating ridged weighted least squares.
pub(crate) fn fit_rank(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
) -> Result<RankFit> {
    let (m, n) = xmat.shape();
    if rnk == 0 {
        return Ok(RankFit {
            a: DMatrix::zeros(m, 1),
            b: DMatrix::zeros(1, n),
        });
    }
    let k = m.min(n);
    if rnk > k {
        return Err(PopmixError::Validation(format!(
            "rank {rnk} exceeds the f4 matrix's {k} singular values"
        )));
    }

    // Start: the rnk right singular vectors with the largest singular values.
    let svd = xmat.clone().svd(false, true);
    let v_t = svd
        .v_t
        .ok_or_else(|| PopmixError::Numerical("SVD of the f4 matrix failed".to_string()))?;
    let mut order: Vec<usize> = (0..svd.singular_values.len()).collect();
    order.sort_by(|&i, &j| {
        svd.singular_values[j]
            .partial_cmp(&svd.singular_values[i])
            .unwrap_or(Ordering::Equal)
    });
    let mut b = DMatrix::from_fn(rnk, n, |r, j| v_t[(order[r], j)]);
    let mut a = xmat * b.transpose();

    let x = vec_row_major(xmat);
    for _ in 0..FIT_ROUNDS {
        // A-step: vec(AB) = D_A vec(A), D_A[(i·n + j), (i·r + c)] = B[c, j].
        let mut d_a = DMatrix::zeros(m * n, m * rnk);
        for i in 0..m {
            for j in 0..n {
                for c in 0..rnk {
                    d_a[(i * n + j, i * rnk + c)] = b[(c, j)];
                }
            }
        }
        let theta = ridged_wls(&d_a, qinv, &x, fudge)?;
        a = DMatrix::from_fn(m, rnk, |i, c| theta[i * rnk + c]);

        // B-step: vec(AB) = D_B vec(B), D_B[(i·n + j), (c·n + j)] = A[i, c].
        let mut d_b = DMatrix::zeros(m * n, rnk * n);
        for i in 0..m {
            for j in 0..n {
                for c in 0..rnk {
                    d_b[(i * n + j, c * n + j)] = a[(i, c)];
                }
            }
        }
        let theta = ridged_wls(&d_b, qinv, &x, fudge)?;
        b = DMatrix::from_fn(rnk, n, |c, j| theta[c * n + j]);
    }
    Ok(RankFit { a, b })
}

/// The misfit chi² of a rank fit and its degrees of freedom (m − r)(n − r)
/// (spec section 3).
pub(crate) fn evaluate_fit(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    fit: &RankFit,
    rnk: usize,
) -> (f64, usize) {
    let (m, n) = xmat.shape();
    let resid = vec_row_major(&(xmat - &fit.a * &fit.b));
    let chisq = resid.dot(&(qinv * &resid));
    let dof = m.saturating_sub(rnk) * n.saturating_sub(rnk);
    (chisq, dof)
}

/// Upper tail of the chi-squared distribution, with R `pchisq(x, df,
/// lower.tail = FALSE)`'s conventions at the edges: with df 0 (a point mass
/// at zero) it is 1 for x < 0 and 0 otherwise; for x <= 0 it is 1.
fn chisq_upper_tail(x: f64, df: usize) -> f64 {
    if df == 0 {
        return if x < 0.0 { 1.0 } else { 0.0 };
    }
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 1.0;
    }
    match ChiSquared::new(df as f64) {
        Ok(dist) => dist.sf(x),
        Err(_) => f64::NAN,
    }
}

/// p-value of a fit: the chi-squared upper tail, NaN for a saturated model
/// (0 degrees of freedom), which is not tested.
pub(crate) fn chisq_p_value(chisq: f64, dof: usize) -> f64 {
    if dof == 0 {
        f64::NAN
    } else {
        chisq_upper_tail(chisq, dof)
    }
}

/// Order p-values largest first, with NaN (untested) after every number.
pub(crate) fn by_p_descending(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => b.partial_cmp(&a).unwrap_or(Ordering::Equal),
    }
}

/// qpAdm weights from the fit at rank `rnk` (spec section 4).
pub(crate) fn qpadm_weights(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
    constrained: bool,
) -> Result<Vec<f64>> {
    let fit = fit_rank(xmat, qinv, rnk, fudge)?;
    weights_from_fit(&fit, rnk, constrained)
}

/// The weights w with wᵀA = 0 and Σw = 1, in the least-squares sense
/// (spec section 4): with C = [Aᵀ; 1ᵀ] and e = (0, ..., 0, 1), solve
/// (CᵀC) w = Cᵀe, or (constrained) the same quadratic over the simplex
/// 0 <= w <= 1, Σw = 1; then normalise to sum 1.
pub(crate) fn weights_from_fit(fit: &RankFit, rnk: usize, constrained: bool) -> Result<Vec<f64>> {
    let m = fit.a.nrows();
    if rnk == 0 && m == 1 {
        return Ok(vec![1.0]);
    }
    let r = fit.a.ncols();
    let mut c = DMatrix::from_element(r + 1, m, 1.0);
    for k in 0..r {
        for i in 0..m {
            c[(k, i)] = fit.a[(i, k)];
        }
    }
    let mut e = DVector::zeros(r + 1);
    e[r] = 1.0;
    let ctc = c.transpose() * &c;
    let cte = c.transpose() * e;

    let w = if constrained {
        simplex_qp(&ctc, &cte).map_err(|err| match err {
            PopmixError::Infeasible { .. } => PopmixError::Infeasible {
                context: "constrained qpAdm weights".to_string(),
            },
            other => other,
        })?
    } else {
        ctc.lu()
            .solve(&cte)
            .filter(|s| s.iter().all(|v| v.is_finite()))
            .ok_or_else(|| singular("qpAdm weights"))?
    };
    let total = w.sum();
    Ok(w.iter().map(|v| v / total).collect())
}

/// Minimise ½ wᵀDw − dᵀw subject to Σw = 1 and 0 <= w <= 1.
///
/// A primal active-set method: bounds in the working set hold their
/// variables exactly at 0 or 1; the remaining (free) variables solve the
/// equality-constrained problem through its KKT system. A step that would
/// leave the box stops at the first bound it meets and adds it; at a
/// stationary point a bound whose multiplier has the wrong sign is released.
/// D is assumed positive semi-definite (it is CᵀC for the qpAdm weights).
pub(crate) fn simplex_qp(d_mat: &DMatrix<f64>, d_vec: &DVector<f64>) -> Result<DVector<f64>> {
    let m = d_vec.len();
    let infeasible = || PopmixError::Infeasible {
        context: "simplex quadratic program".to_string(),
    };
    if m == 0 || d_mat.shape() != (m, m) {
        return Err(infeasible());
    }
    if m == 1 {
        return Ok(DVector::from_element(1, 1.0));
    }

    #[derive(Clone, Copy, PartialEq)]
    enum State {
        Free,
        Lower,
        Upper,
    }
    let mut state = vec![State::Free; m];
    let mut w = DVector::from_element(m, 1.0 / m as f64);
    let tol = 1e-12;
    let max_iter = 50 * (m + 1);

    for _ in 0..max_iter {
        let free: Vec<usize> = (0..m).filter(|&i| state[i] == State::Free).collect();
        let nf = free.len();
        let fixed_sum = state.iter().filter(|&&s| s == State::Upper).count() as f64;

        // Equality-constrained optimum over the free variables:
        // D_FF w_F − λ1 = d_F − D_FX w_X,   1ᵀ w_F = 1 − Σ w_X.
        let (target_free, lambda) = {
            let mut kkt = DMatrix::zeros(nf + 1, nf + 1);
            let mut rhs = DVector::zeros(nf + 1);
            for (a, &i) in free.iter().enumerate() {
                for (b, &j) in free.iter().enumerate() {
                    kkt[(a, b)] = d_mat[(i, j)];
                }
                kkt[(a, nf)] = -1.0;
                kkt[(nf, a)] = 1.0;
                let mut v = d_vec[i];
                for j in 0..m {
                    if state[j] != State::Free {
                        v -= d_mat[(i, j)] * w[j];
                    }
                }
                rhs[a] = v;
            }
            rhs[nf] = 1.0 - fixed_sum;
            let sol = kkt
                .lu()
                .solve(&rhs)
                .filter(|s| s.iter().all(|v| v.is_finite()))
                .ok_or_else(|| singular("constrained qpAdm weights"))?;
            let mut target: Vec<f64> = (0..nf).map(|a| sol[a]).collect();
            if nf == 1 {
                // A single free variable is pinned by the sum constraint.
                target[0] = 1.0 - fixed_sum;
            }
            (target, sol[nf])
        };

        // Step from w towards the free optimum, stopping at the first bound.
        let mut alpha = 1.0;
        let mut blocking: Option<(usize, State)> = None;
        for (a, &i) in free.iter().enumerate() {
            let p = target_free[a] - w[i];
            if p < -tol && target_free[a] < 0.0 {
                let t = w[i] / -p;
                if t < alpha {
                    alpha = t;
                    blocking = Some((i, State::Lower));
                }
            } else if p > tol && target_free[a] > 1.0 {
                let t = (1.0 - w[i]) / p;
                if t < alpha {
                    alpha = t;
                    blocking = Some((i, State::Upper));
                }
            }
        }
        for (a, &i) in free.iter().enumerate() {
            w[i] += alpha * (target_free[a] - w[i]);
        }
        if let Some((i, s)) = blocking {
            state[i] = s;
            w[i] = if s == State::Upper { 1.0 } else { 0.0 };
            if s == State::Upper {
                // Σw = 1 with one weight at 1 forces every other to 0.
                for j in 0..m {
                    if j != i {
                        state[j] = State::Lower;
                        w[j] = 0.0;
                    }
                }
                // Keep one variable free so the multiplier is defined.
                state[i] = State::Free;
            }
            continue;
        }
        for (a, &i) in free.iter().enumerate() {
            w[i] = target_free[a].clamp(0.0, 1.0);
        }

        // Stationary on the working set: check the bound multipliers.
        // Gradient of the Lagrangian g = Dw − d − λ1; at a lower bound it
        // must be >= 0, at an upper bound <= 0.
        let g = d_mat * &w - d_vec;
        let mut worst: Option<(usize, f64)> = None;
        for i in 0..m {
            let gi = g[i] - lambda;
            let violation = match state[i] {
                State::Free => 0.0,
                State::Lower => -gi,
                State::Upper => gi,
            };
            if violation > tol * (1.0 + g[i].abs().max(lambda.abs()))
                && worst.is_none_or(|(_, v)| violation > v)
            {
                worst = Some((i, violation));
            }
        }
        match worst {
            None => {
                if w.iter().any(|v| !v.is_finite()) {
                    return Err(infeasible());
                }
                return Ok(w);
            }
            Some((i, _)) => state[i] = State::Free,
        }
    }
    Err(PopmixError::ConvergenceFailure {
        iterations: max_iter,
        context: "constrained qpAdm weights".to_string(),
    })
}

/// Jackknife standard errors of the qpAdm weights (spec section 5): the
/// weights of every leave-one-block-out f4 matrix with the same W, then
/// sqrt(Σ_k (w_k − w̄)²) per source. NaN with fewer than two replicates.
pub(crate) fn qpadm_weight_std_errors(
    leave_one_out: &[Array2<f64>],
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
    constrained: bool,
) -> Result<Vec<f64>> {
    if leave_one_out.len() < 2 {
        let m = leave_one_out.first().map_or(rnk + 1, |x| x.nrows());
        return Ok(vec![f64::NAN; m]);
    }
    let reps = leave_one_out
        .iter()
        .map(|x| qpadm_weights(&to_dmatrix(x), qinv, rnk, fudge, constrained))
        .collect::<Result<Vec<_>>>()?;
    let m = reps[0].len();
    let nrep = reps.len() as f64;
    Ok((0..m)
        .map(|s| {
            let mean = reps.iter().map(|w| w[s]).sum::<f64>() / nrep;
            reps.iter()
                .map(|w| (w[s] - mean).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect())
}

/// One model of qpAdm's `popdrop`: the fit to a subset of the sources.
#[derive(Debug, Clone)]
pub struct PopDropRow {
    /// One character per source: '0' kept, '1' dropped.
    pub pattern: String,
    /// (sources − 1) − f4rank.
    pub wt: usize,
    /// Degrees of freedom of the fit.
    pub dof: usize,
    /// Chi-squared of the fit.
    pub chi_sq: f64,
    /// p-value of the fit (NaN when dof is 0).
    pub p_value: f64,
    /// Rank of the fit (kept sources − 1).
    pub f4rank: usize,
    /// Per source: its weight, or None if dropped.
    pub weights: Vec<Option<f64>>,
    /// Whether every weight lies in [0, 1].
    pub feasible: bool,
    /// Whether the model is on the nested path.
    pub best: bool,
    /// dof difference from the model it is compared with.
    pub dofdiff: Option<usize>,
    /// Chi-squared difference from the model it is compared with.
    pub chisqdiff: Option<f64>,
    /// p-value of that nested comparison.
    pub p_nested: Option<f64>,
}

/// All k-subsets of 0..n in lexicographic order (as R `combn`).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if k > n {
        return out;
    }
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        out.push(idx.clone());
        // Rightmost position that can still advance.
        let Some(pos) = (0..k).rev().find(|&p| idx[p] < n - k + p) else {
            return out;
        };
        idx[pos] += 1;
        for p in pos + 1..k {
            idx[p] = idx[p - 1] + 1;
        }
    }
}

/// Fit every non-empty subset of the sources (spec section 6).
pub(crate) fn drop_pops(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    fudge: f64,
    constrained: bool,
    compat: Admixtools2Compat,
) -> Result<Vec<PopDropRow>> {
    let (m, n) = xmat.shape();
    if m == 0 {
        return Ok(Vec::new());
    }
    let full_rank = m - 1;
    let mut rows = Vec::with_capacity((1usize << m.min(63)) - 1);
    for size in 1..=m {
        for subset in combinations(m, size) {
            let xs = DMatrix::from_fn(size, n, |i, j| xmat[(subset[i], j)]);
            let entries: Vec<usize> = subset
                .iter()
                .flat_map(|&i| (0..n).map(move |j| i * n + j))
                .collect();
            let ws = DMatrix::from_fn(entries.len(), entries.len(), |a, b| {
                qinv[(entries[a], entries[b])]
            });
            let rnk = size - 1;
            let fit = fit_rank(&xs, &ws, rnk, fudge)?;
            let (chi_sq, dof) = evaluate_fit(&xs, &ws, &fit, rnk);
            let w = weights_from_fit(&fit, rnk, constrained)?;
            let feasible = w.iter().all(|v| (0.0..=1.0).contains(v));
            let mut weights = vec![None; m];
            let mut pattern = vec!['1'; m];
            for (k, &i) in subset.iter().enumerate() {
                weights[i] = Some(w[k]);
                pattern[i] = '0';
            }
            rows.push(PopDropRow {
                pattern: pattern.into_iter().collect(),
                wt: full_rank - rnk,
                dof,
                chi_sq,
                p_value: chisq_p_value(chi_sq, dof),
                f4rank: rnk,
                weights,
                feasible,
                best: false,
                dofdiff: None,
                chisqdiff: None,
                p_nested: None,
            });
        }
    }
    mark_nested_path(&mut rows, full_rank, compat);
    rows.sort_by(|a, b| a.dof.cmp(&b.dof).then_with(|| a.pattern.cmp(&b.pattern)));
    Ok(rows)
}

/// Record on `rows[target]` its comparison with model `other`.
fn record_comparison(rows: &mut [PopDropRow], target: usize, other: usize) {
    let dofdiff = rows[target].dof.abs_diff(rows[other].dof);
    let chisqdiff = rows[target].chi_sq - rows[other].chi_sq;
    let t = &mut rows[target];
    t.dofdiff = Some(dofdiff);
    t.chisqdiff = Some(chisqdiff);
    t.p_nested = Some(chisq_upper_tail(chisqdiff, dofdiff));
}

/// Mark the `best` models and their nested comparisons (spec section 6).
/// `rows` are in subset order; `full_rank` is (sources − 1).
fn mark_nested_path(rows: &mut [PopDropRow], full_rank: usize, compat: Admixtools2Compat) {
    if compat.popdrop_sibling_path {
        // admixtools2 2.0.10: every one-drop model, each compared with the
        // next in order of decreasing dof.
        if full_rank == 0 {
            return;
        }
        let mut sib: Vec<usize> = (0..rows.len())
            .filter(|&i| rows[i].f4rank == full_rank - 1)
            .collect();
        sib.sort_by(|&a, &b| rows[b].dof.cmp(&rows[a].dof));
        for &i in &sib {
            rows[i].best = true;
        }
        for pair in sib.windows(2) {
            record_comparison(rows, pair[0], pair[1]);
        }
        return;
    }

    // Corrected: a greedy nested path from the full model downwards.
    let Some(mut current) = rows
        .iter()
        .position(|r| r.f4rank == full_rank && r.pattern.bytes().all(|c| c == b'0'))
    else {
        return;
    };
    rows[current].best = true;
    while rows[current].f4rank > 0 {
        let cur_rank = rows[current].f4rank;
        let cur_pat = rows[current].pattern.as_bytes().to_vec();
        let mut chosen: Option<usize> = None;
        for (i, r) in rows.iter().enumerate() {
            let nested = r.f4rank + 1 == cur_rank
                && r.feasible
                && r.pattern.len() == cur_pat.len()
                && r.pattern
                    .bytes()
                    .zip(&cur_pat)
                    .all(|(c, &p)| p != b'1' || c == b'1');
            if nested && chosen.is_none_or(|c| r.chi_sq < rows[c].chi_sq) {
                chosen = Some(i);
            }
        }
        let Some(child) = chosen else {
            break;
        };
        rows[child].best = true;
        record_comparison(rows, child, current);
        current = child;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn untested_models_rank_after_tested_ones() {
        let mut p = [0.2, f64::NAN, 0.9, 0.01, f64::NAN, 0.5];
        p.sort_by(|a, b| by_p_descending(*a, *b));
        assert_eq!(&p[..4], &[0.9, 0.5, 0.2, 0.01]);
        assert!(p[4..].iter().all(|x| x.is_nan()));
    }

    #[test]
    fn a_saturated_model_has_no_p_value() {
        assert!(chisq_p_value(0.3, 0).is_nan());
        assert!(chisq_p_value(0.0, 0).is_nan());
        // R pchisq for nested differences keeps its point mass at 0
        assert_eq!(chisq_upper_tail(0.3, 0), 0.0);
    }

    #[test]
    fn test_rank_one_matrix_fits_rank_one_exactly() {
        // x = u v' has rank 1, so the rank-1 residual is (near) zero.
        let x = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0]);
        let qinv = DMatrix::<f64>::identity(6, 6);
        let fit = fit_rank(&x, &qinv, 1, 0.0).unwrap();
        let (chisq, dof) = evaluate_fit(&x, &qinv, &fit, 1);
        assert!(chisq < 1e-20, "chisq {chisq}");
        assert_eq!(dof, 2);
    }

    #[test]
    fn test_rank_zero_chisq_is_the_full_quadratic_form() {
        let x = DMatrix::from_row_slice(1, 2, &[0.5, -1.0]);
        let qinv = DMatrix::from_row_slice(2, 2, &[2.0, 0.5, 0.5, 1.0]);
        let fit = fit_rank(&x, &qinv, 0, 1e-4).unwrap();
        let (chisq, dof) = evaluate_fit(&x, &qinv, &fit, 0);
        // [0.5, -1] Q [0.5, -1]' = 0.5 - 0.5 + 1
        assert_relative_eq!(chisq, 1.0, max_relative = 1e-14);
        assert_eq!(dof, 2);
    }

    #[test]
    fn test_fudge_scales_with_trace() {
        let cov = ndarray::arr2(&[[2.0, 0.0], [0.0, 6.0]]);
        let q = fudged_inverse(&cov, 0.5).unwrap();
        // diag + 0.5 * 8 = [6, 10]
        assert_relative_eq!(q[(0, 0)], 1.0 / 6.0, max_relative = 1e-14);
        assert_relative_eq!(q[(1, 1)], 1.0 / 10.0, max_relative = 1e-14);
    }

    #[test]
    fn test_simplex_qp_matches_unconstrained_inside_and_clamps_outside() {
        let d = DMatrix::from_row_slice(2, 2, &[2.0, 0.5, 0.5, 1.0]);
        // Unconstrained optimum with sum 1 inside the box: same as the KKT solution
        let w = simplex_qp(&d, &DVector::from_vec(vec![1.0, 1.0])).unwrap();
        assert_relative_eq!(w.sum(), 1.0, max_relative = 1e-14);
        assert!(w.iter().all(|&v| (0.0..=1.0).contains(&v)));
        // A strong pull towards the first weight pins it at 1
        let w = simplex_qp(&d, &DVector::from_vec(vec![10.0, -10.0])).unwrap();
        assert_eq!(w.as_slice(), &[1.0, 0.0]);
    }

    #[test]
    fn test_combinations_follow_r_combn_order() {
        assert_eq!(combinations(3, 2), vec![vec![0, 1], vec![0, 2], vec![1, 2]]);
        assert_eq!(combinations(3, 3), vec![vec![0, 1, 2]]);
        assert_eq!(combinations(2, 1), vec![vec![0], vec![1]]);
    }

    #[test]
    fn test_chisq_upper_tail_with_zero_df() {
        // R: pchisq(-3, 0, lower.tail = FALSE) = 1; pchisq(3, 0, ...) = 0
        assert_eq!(chisq_upper_tail(-3.0, 0), 1.0);
        assert_eq!(chisq_upper_tail(3.0, 0), 0.0);
    }

    fn row(pattern: &str, f4rank: usize, dof: usize, chi_sq: f64, feasible: bool) -> PopDropRow {
        PopDropRow {
            pattern: pattern.to_string(),
            wt: 2 - f4rank,
            dof,
            chi_sq,
            p_value: chisq_p_value(chi_sq, dof),
            f4rank,
            weights: Vec::new(),
            feasible,
            best: false,
            dofdiff: None,
            chisqdiff: None,
            p_nested: None,
        }
    }

    /// Three sources, in admixtools2 power_set order. Full model 000; of the
    /// one-drop models 010 fits best but is infeasible, so the path takes
    /// 001; of 001's children (011, 101) 101 fits best. 110 fits better
    /// still but is not nested in 001.
    fn three_source_rows() -> Vec<PopDropRow> {
        vec![
            row("011", 0, 5, 100.0, true),
            row("101", 0, 5, 50.0, true),
            row("110", 0, 5, 30.0, true),
            row("001", 1, 4, 12.0, true),
            row("010", 1, 4, 10.0, false),
            row("100", 1, 4, 15.0, true),
            row("000", 2, 3, 8.0, false),
        ]
    }

    /// (pattern, dofdiff, chisqdiff, p_nested) of each model on the path.
    type Nested<'a> = (&'a str, Option<usize>, Option<f64>, Option<f64>);

    fn nested_columns(rows: &[PopDropRow]) -> Vec<Nested<'_>> {
        rows.iter()
            .filter(|r| r.best)
            .map(|r| (r.pattern.as_str(), r.dofdiff, r.chisqdiff, r.p_nested))
            .collect()
    }

    #[test]
    fn test_corrected_popdrop_path_is_nested_and_greedy() {
        // ADMIXTOOLS2-0001 oracle: each step drops one more source, taking
        // the feasible child with the lowest chi-squared, and is compared
        // with the model it came from (1 degree of freedom).
        let mut rows = three_source_rows();
        mark_nested_path(&mut rows, 2, Admixtools2Compat::CORRECTED);
        let mut got = nested_columns(&rows);
        got.sort_by(|a, b| a.0.cmp(b.0));
        assert_eq!(
            got,
            vec![
                ("000", None, None, None),
                ("001", Some(1), Some(4.0), Some(chisq_upper_tail(4.0, 1))),
                ("101", Some(1), Some(38.0), Some(chisq_upper_tail(38.0, 1))),
            ]
        );
    }

    #[test]
    fn test_corrected_popdrop_path_stops_without_a_feasible_child() {
        let mut rows = three_source_rows();
        for r in rows.iter_mut().filter(|r| r.f4rank == 1) {
            r.feasible = false;
        }
        mark_nested_path(&mut rows, 2, Admixtools2Compat::CORRECTED);
        assert_eq!(nested_columns(&rows), vec![("000", None, None, None)]);
    }

    #[test]
    fn test_v2_0_10_popdrop_compares_the_one_drop_siblings() {
        // ADMIXTOOLS2-0001 compat: admixtools2 2.0.10 marks every model at
        // full rank - 1 and compares each with the next (same dof, so
        // dofdiff 0 and pchisq(x, 0) is 0 for x >= 0, 1 for x < 0).
        let mut rows = three_source_rows();
        mark_nested_path(&mut rows, 2, Admixtools2Compat::V2_0_10);
        assert_eq!(
            nested_columns(&rows),
            vec![
                ("001", Some(0), Some(2.0), Some(0.0)),
                ("010", Some(0), Some(-5.0), Some(1.0)),
                ("100", None, None, None),
            ]
        );
    }

    #[test]
    fn test_popdrop_compat_flag_only_changes_the_nested_columns() {
        // ADMIXTOOLS2-0001 isolation.
        let x = DMatrix::from_fn(3, 5, |i, j| {
            ((i * 5 + j) as f64 * 1.7).sin() + 0.3 * i as f64
        });
        let qinv = DMatrix::from_fn(15, 15, |i, j| if i == j { 2.0 } else { 0.1 });
        let fixed = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::CORRECTED).unwrap();
        let compat = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::V2_0_10).unwrap();
        assert_eq!(fixed.len(), 7);
        for (a, b) in fixed.iter().zip(&compat) {
            assert_eq!(a.pattern, b.pattern);
            assert_eq!(
                (a.wt, a.dof, a.f4rank, a.feasible),
                (b.wt, b.dof, b.f4rank, b.feasible)
            );
            assert_eq!(a.chi_sq.to_bits(), b.chi_sq.to_bits());
            assert_eq!(a.p_value.to_bits(), b.p_value.to_bits());
            assert_eq!(a.weights, b.weights);
        }
        assert_ne!(nested_columns(&fixed), nested_columns(&compat));
    }

    #[test]
    fn test_constrained_popdrop_with_three_sources() {
        // ADMIXTOOLS2-0002: admixtools2 qpadm(constrained = TRUE) stops here.
        // Every subset model is fitted at full rank, where the constrained
        // weights are well defined, and the chi-squareds do not depend on
        // the weights at all.
        let x = DMatrix::from_fn(3, 5, |i, j| {
            ((i * 5 + j) as f64 * 1.7).sin() + 0.3 * i as f64
        });
        let qinv = DMatrix::from_fn(15, 15, |i, j| if i == j { 2.0 } else { 0.1 });
        let free = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::CORRECTED).unwrap();
        let con = drop_pops(&x, &qinv, 1e-4, true, Admixtools2Compat::CORRECTED).unwrap();
        for (f, c) in free.iter().zip(&con) {
            assert_eq!(f.chi_sq.to_bits(), c.chi_sq.to_bits());
            assert!(
                c.feasible,
                "pat {}: constrained weights {:?}",
                c.pattern, c.weights
            );
            let sum: f64 = c.weights.iter().flatten().sum();
            assert!((sum - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_p_value_tail_is_exact() {
        // pchisq(129.6103737533115, 4, lower.tail = FALSE) in R
        assert_relative_eq!(
            chisq_p_value(129.6103737533115, 4),
            4.717642937981588e-27,
            max_relative = 1e-10
        );
    }
}
