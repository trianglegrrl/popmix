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
    pub(crate) a: DMatrix<f64>,
    pub(crate) b: DMatrix<f64>,
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
        return Err(PopmixError::parameter(
            "rank",
            format!("{rnk} exceeds the f4 matrix's {k} singular values"),
        ));
    }

    // Start: the rnk right singular vectors with the largest singular values.
    let svd = xmat.clone().svd(false, true);
    let v_t = svd
        .v_t
        .ok_or_else(|| PopmixError::numerical("qpAdm", "SVD of the f4 matrix failed"))?;

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
mod tests;
