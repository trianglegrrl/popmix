//! qpWave/qpAdm rank fitting, as in admixtools2.
//!
//! Reproduces admixtools2 (`qpadm`, `qpadm_fit`, `qpadm_evaluate_fit`, and
//! the C++ `cpp_qpadm_weights`, `cpp_opt_A`, `cpp_opt_B`):
//!
//! 1. Q = covariance of the f4 matrix with `fudge * trace(Q)` added to the
//!    diagonal, then inverted.
//! 2. A rank-r model f4 ≈ A B (A: rows × r, B: r × cols) starts from the top r
//!    right singular vectors of the f4 matrix and alternates covariance-
//!    weighted least squares for A and for B, 20 times. Each solve also adds
//!    `fudge * trace` to its coefficient diagonal.
//! 3. chi² = res' Q⁻¹ res for the residual res = f4 - A B, vectorised by rows,
//!    with (rows - r)(cols - r) degrees of freedom.
//! 4. qpAdm weights w solve [A 1]' w = (0, ..., 0, 1) in least squares, then
//!    are scaled to sum to 1. Constrained weights (0 <= w <= 1, sum 1) solve
//!    the same quadratic program admixtools2 passes to quadprog.
//! 5. qpAdm weight standard errors refit the weights on every leave-one-block-
//!    out f4 matrix with the same Q⁻¹ (`cpp_get_weights_covariance`).

use nalgebra::{DMatrix, DVector};
use ndarray::Array2;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::{PopmixError, Result};

/// admixtools2 qpadm_weights `iterations`.
const ITERATIONS: usize = 20;

pub(crate) fn to_dmatrix(a: &Array2<f64>) -> DMatrix<f64> {
    DMatrix::from_fn(a.nrows(), a.ncols(), |i, j| a[[i, j]])
}

/// Row-major vector of a matrix (armadillo `vectorise(x.t())`, R `c(t(x))`).
fn rows_vec(x: &DMatrix<f64>) -> DVector<f64> {
    DVector::from_iterator(x.len(), x.transpose().iter().copied())
}

fn add_fudge(m: &mut DMatrix<f64>, fudge: f64) {
    let add = fudge * m.trace();
    for i in 0..m.nrows() {
        m[(i, i)] += add;
    }
}

fn solve(a: DMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>> {
    a.lu()
        .solve(b)
        .ok_or_else(|| PopmixError::Numerical("singular system in qpWave/qpAdm fit".to_string()))
}

/// Inverse of the f4 covariance after adding `fudge * trace` to its diagonal.
pub(crate) fn fudged_inverse(cov: &Array2<f64>, fudge: f64) -> Result<DMatrix<f64>> {
    let mut q = to_dmatrix(cov);
    add_fudge(&mut q, fudge);
    q.try_inverse()
        .ok_or_else(|| PopmixError::Numerical("f4 covariance matrix is singular".to_string()))
}

/// cpp_opt_A: A minimising the Q-weighted error for fixed B.
fn opt_a(
    b: &DMatrix<f64>,
    xvec: &DVector<f64>,
    qinv: &DMatrix<f64>,
    nr: usize,
    fudge: f64,
) -> Result<DMatrix<f64>> {
    let b2 = DMatrix::<f64>::identity(nr, nr).kronecker(b);
    let mut coeffs = &b2 * qinv * b2.transpose();
    let rhs = &b2 * qinv * xvec;
    add_fudge(&mut coeffs, fudge);
    let a2 = solve(coeffs, &rhs)?;
    // mat(&A2[0], B.n_rows, nr).t(): column-major (rnk x nr), transposed
    Ok(DMatrix::from_column_slice(b.nrows(), nr, a2.as_slice()).transpose())
}

/// cpp_opt_B: B minimising the Q-weighted error for fixed A.
fn opt_b(
    a: &DMatrix<f64>,
    xvec: &DVector<f64>,
    qinv: &DMatrix<f64>,
    nc: usize,
    fudge: f64,
) -> Result<DMatrix<f64>> {
    let a2 = a.kronecker(&DMatrix::<f64>::identity(nc, nc));
    let mut coeffs = a2.transpose() * qinv * &a2;
    let rhs = a2.transpose() * qinv * xvec;
    add_fudge(&mut coeffs, fudge);
    let b2 = solve(coeffs, &rhs)?;
    // mat(&B2[0], nc, A.n_cols).t()
    Ok(DMatrix::from_column_slice(nc, a.ncols(), b2.as_slice()).transpose())
}

/// A rank-r fit of the f4 matrix: f4 ≈ A B.
pub(crate) struct RankFit {
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
}

/// cpp_qpadm_weights without the weights: the alternating fit for rank `rnk`.
/// Rank 0 is A = 0, B = 0.
pub(crate) fn fit_rank(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
) -> Result<RankFit> {
    let (nr, nc) = xmat.shape();
    if rnk == 0 {
        return Ok(RankFit {
            a: DMatrix::zeros(nr, 1),
            b: DMatrix::zeros(1, nc),
        });
    }
    let svd = xmat.clone().svd(false, true);
    let v_t = svd
        .v_t
        .ok_or_else(|| PopmixError::Numerical("SVD of f4 matrix failed".to_string()))?;
    // Singular values in descending order, as armadillo returns them
    let mut order: Vec<usize> = (0..svd.singular_values.len()).collect();
    order.sort_by(|&i, &j| svd.singular_values[j].total_cmp(&svd.singular_values[i]));
    if rnk > order.len() {
        return Err(PopmixError::Validation(format!(
            "rank {rnk} exceeds the f4 matrix's {} singular values",
            order.len()
        )));
    }
    let mut b = DMatrix::from_fn(rnk, nc, |i, j| v_t[(order[i], j)]);
    let mut a = xmat * b.transpose();
    let xvec = rows_vec(xmat);
    for _ in 0..ITERATIONS {
        a = opt_a(&b, &xvec, qinv, nr, fudge)?;
        b = opt_b(&a, &xvec, qinv, nc, fudge)?;
    }
    Ok(RankFit { a, b })
}

/// qpadm_evaluate_fit: chi-squared and degrees of freedom of a rank-r fit.
pub(crate) fn evaluate_fit(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    fit: &RankFit,
    rnk: usize,
) -> (f64, usize) {
    let res = rows_vec(&(xmat - &fit.a * &fit.b));
    let chisq = (res.transpose() * qinv * &res)[(0, 0)];
    let (nr, nc) = xmat.shape();
    let dof = nr.saturating_sub(rnk) * nc.saturating_sub(rnk);
    (chisq, dof)
}

/// cpp_qpadm_weights: admixture weights of the rank `rnk` fit, summing to 1.
/// Rank 0 (one source) has weight 1.
pub(crate) fn qpadm_weights(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
    constrained: bool,
) -> Result<Vec<f64>> {
    if rnk == 0 {
        return Ok(vec![1.0]);
    }
    weights_from_fit(&fit_rank(xmat, qinv, rnk, fudge)?, rnk, constrained)
}

/// The weights step of cpp_qpadm_weights for an existing rank `rnk` fit.
fn weights_from_fit(fit: &RankFit, rnk: usize, constrained: bool) -> Result<Vec<f64>> {
    if rnk == 0 {
        return Ok(vec![1.0]);
    }
    let nr = fit.a.nrows();
    // x = [A 1]' ((rnk + 1) x nr), y = (0, ..., 0, 1)
    let x = DMatrix::from_fn(
        rnk + 1,
        nr,
        |i, j| if i < rnk { fit.a[(j, i)] } else { 1.0 },
    );
    let mut y = DVector::zeros(rnk + 1);
    y[rnk] = 1.0;
    let rhs = x.transpose() * &x;
    let lhs = x.transpose() * &y;
    let w = if constrained {
        simplex_qp(&rhs, &lhs)?
    } else {
        solve(rhs, &lhs)?
    };
    let total = w.sum();
    Ok(w.iter().map(|v| v / total).collect())
}

/// One row of admixtools2 qpadm `popdrop`: the model fitted to a subset of
/// the sources.
#[derive(Debug, Clone)]
pub struct PopDropRow {
    /// One character per source: '0' if used, '1' if dropped.
    pub pattern: String,
    /// Sources dropped relative to the full model: full rank - `f4rank`.
    pub wt: usize,
    /// Degrees of freedom.
    pub dof: usize,
    /// Chi-squared of the subset model.
    pub chi_sq: f64,
    /// P-value (1 when `dof` is 0, where admixtools2 reports NA).
    pub p_value: f64,
    /// Rank of the subset model: number of sources used - 1.
    pub f4rank: usize,
    /// Weight of each source, `None` for dropped sources.
    pub weights: Vec<Option<f64>>,
    /// Whether every weight lies in [0, 1].
    pub feasible: bool,
    /// Whether admixtools2 marks this model as on its nested path.
    pub best: bool,
    /// Nested-path comparison with the next model on the path.
    pub dofdiff: Option<usize>,
    /// Chi-squared difference with the next model on the path.
    pub chisqdiff: Option<f64>,
    /// P-value of `chisqdiff` with `dofdiff` degrees of freedom.
    pub p_nested: Option<f64>,
}

/// admixtools2 `drop_pops`: fit every non-empty subset of the sources (rows
/// of `xmat`), using the matching block of the full model's Q⁻¹.
///
/// The nested path reproduces admixtools2 2.0.10 as it behaves, not as it
/// is presumably meant to: its `child_patterns` only finds children when one
/// model is on the path, but the path starts with every model at rank
/// (full rank - 1), so for two or more sources nothing below that rank is
/// ever added. `best`, `dofdiff`, `chisqdiff` and `p_nested` therefore cover
/// only the models that drop exactly one source, compared with each other in
/// order. See DEVIATIONS.md.
pub(crate) fn drop_pops(
    xmat: &DMatrix<f64>,
    qinv: &DMatrix<f64>,
    fudge: f64,
    constrained: bool,
) -> Result<Vec<PopDropRow>> {
    let (nr, nc) = xmat.shape();
    let full_rank = nr - 1;
    let mut rows = Vec::new();
    // power_set: subsets by size, each size in combn (lexicographic) order
    for size in 1..=nr {
        for subset in combinations(nr, size) {
            let sub_x = DMatrix::from_fn(size, nc, |i, j| xmat[(subset[i], j)]);
            let idx: Vec<usize> = subset
                .iter()
                .flat_map(|&r| (0..nc).map(move |c| r * nc + c))
                .collect();
            let sub_q = DMatrix::from_fn(idx.len(), idx.len(), |i, j| qinv[(idx[i], idx[j])]);
            let rnk = size - 1;
            let fit = fit_rank(&sub_x, &sub_q, rnk, fudge)?;
            let (chi_sq, dof) = evaluate_fit(&sub_x, &sub_q, &fit, rnk);
            let w = weights_from_fit(&fit, rnk, constrained)?;
            let mut weights = vec![None; nr];
            for (k, &r) in subset.iter().enumerate() {
                weights[r] = Some(w[k]);
            }
            rows.push(PopDropRow {
                pattern: (0..nr)
                    .map(|r| if subset.contains(&r) { '0' } else { '1' })
                    .collect(),
                wt: full_rank - rnk,
                dof,
                chi_sq,
                p_value: chisq_p_value(chi_sq, dof),
                f4rank: rnk,
                feasible: w.iter().all(|v| (0.0..=1.0).contains(v)),
                weights,
                best: false,
                dofdiff: None,
                chisqdiff: None,
                p_nested: None,
            });
        }
    }

    // The nested path: models at rank full_rank - 1 (see above), arranged by
    // decreasing dof (stable), each compared with the next.
    if full_rank >= 1 {
        let mut nested: Vec<usize> = (0..rows.len())
            .filter(|&i| rows[i].f4rank == full_rank - 1)
            .collect();
        nested.sort_by_key(|&i| std::cmp::Reverse(rows[i].dof));
        for (k, &i) in nested.iter().enumerate() {
            rows[i].best = true;
            if let Some(&next) = nested.get(k + 1) {
                let dofdiff = rows[i].dof - rows[next].dof;
                let chisqdiff = rows[i].chi_sq - rows[next].chi_sq;
                rows[i].dofdiff = Some(dofdiff);
                rows[i].chisqdiff = Some(chisqdiff);
                rows[i].p_nested = Some(chisq_upper_tail(chisqdiff, dofdiff));
            }
        }
    }

    // arrange(dof, pat)
    rows.sort_by(|a, b| a.dof.cmp(&b.dof).then_with(|| a.pattern.cmp(&b.pattern)));
    Ok(rows)
}

/// All `k`-element subsets of 0..n in lexicographic order (R `combn`).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        out.push(idx.clone());
        let Some(i) = (0..k).rev().find(|&i| idx[i] != i + n - k) else {
            return out;
        };
        idx[i] += 1;
        for j in i + 1..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

/// R `pchisq(x, df, lower.tail = FALSE)`, including df = 0 (a point mass
/// at 0) and negative x.
fn chisq_upper_tail(x: f64, df: usize) -> f64 {
    if df == 0 {
        return if x < 0.0 { 1.0 } else { 0.0 };
    }
    if x <= 0.0 {
        return 1.0;
    }
    chisq_p_value(x, df)
}

/// Minimise ½ w'Dw - d'w subject to sum(w) = 1 and 0 <= w <= 1: the problem
/// admixtools2 gives quadprog for constrained qpAdm.
///
/// Solved exactly by trying every assignment of each weight to 0, 1 or free
/// and solving the equality-constrained problem on the free weights; the
/// feasible candidate with the lowest objective is the optimum (D is positive
/// semi-definite, so the problem is convex). 3^n candidates: fine for the few
/// sources qpAdm models have.
fn simplex_qp(d_mat: &DMatrix<f64>, d_vec: &DVector<f64>) -> Result<DVector<f64>> {
    const TOL: f64 = 1e-12;
    let n = d_vec.len();
    let objective = |w: &DVector<f64>| 0.5 * (w.transpose() * d_mat * w)[(0, 0)] - d_vec.dot(w);
    let mut best: Option<(f64, DVector<f64>)> = None;
    for code in 0..3usize.pow(n as u32) {
        // state per weight: 0 = free, 1 = fixed at 0, 2 = fixed at 1
        let state: Vec<usize> = (0..n).map(|i| code / 3usize.pow(i as u32) % 3).collect();
        let free: Vec<usize> = (0..n).filter(|&i| state[i] == 0).collect();
        let mut w: DVector<f64> = DVector::from_fn(n, |i, _| if state[i] == 2 { 1.0 } else { 0.0 });
        let remaining = 1.0 - w.sum();
        if free.is_empty() {
            if remaining.abs() > TOL {
                continue;
            }
        } else {
            // KKT: [D_ff 1; 1' 0] [w_f; mu] = [d_f - D_fx w_x; remaining]
            let m = free.len();
            let fixed_term = d_mat * &w;
            let kkt = DMatrix::from_fn(m + 1, m + 1, |i, j| match (i < m, j < m) {
                (true, true) => d_mat[(free[i], free[j])],
                (true, false) | (false, true) => 1.0,
                (false, false) => 0.0,
            });
            let rhs = DVector::from_fn(m + 1, |i, _| {
                if i < m {
                    d_vec[free[i]] - fixed_term[free[i]]
                } else {
                    remaining
                }
            });
            let Some(sol) = kkt.lu().solve(&rhs) else {
                continue;
            };
            for (k, &i) in free.iter().enumerate() {
                w[i] = sol[k];
            }
            if free.iter().any(|&i| w[i] < -TOL || w[i] > 1.0 + TOL) {
                continue;
            }
        }
        let obj = objective(&w);
        if best.as_ref().is_none_or(|(b, _)| obj < *b) {
            best = Some((obj, w));
        }
    }
    best.map(|(_, w)| w)
        .ok_or_else(|| PopmixError::Numerical("no feasible constrained qpAdm weights".to_string()))
}

/// cpp_get_weights_covariance (jackknife): weight standard errors from the
/// weights refitted on each leave-one-block-out f4 matrix.
pub(crate) fn qpadm_weight_std_errors(
    leave_one_out: &[Array2<f64>],
    qinv: &DMatrix<f64>,
    rnk: usize,
    fudge: f64,
    constrained: bool,
) -> Result<Vec<f64>> {
    let numreps = leave_one_out.len();
    let wmat = leave_one_out
        .iter()
        .map(|lo| qpadm_weights(&to_dmatrix(lo), qinv, rnk, fudge, constrained))
        .collect::<Result<Vec<_>>>()?;
    if numreps < 2 {
        return Ok(vec![f64::NAN; wmat.first().map_or(0, |w| w.len())]);
    }
    // wmat *= sqrt(numreps - 1); var = cov(wmat), the sample covariance
    let scale = ((numreps - 1) as f64).sqrt();
    let n_weights = wmat[0].len();
    Ok((0..n_weights)
        .map(|j| {
            let col: Vec<f64> = wmat.iter().map(|w| w[j] * scale).collect();
            let mean = col.iter().sum::<f64>() / numreps as f64;
            let var = col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (numreps - 1) as f64;
            var.sqrt()
        })
        .collect())
}

/// Upper tail probability of a chi-squared distribution. 1 when dof is 0.
pub(crate) fn chisq_p_value(chisq: f64, dof: usize) -> f64 {
    if dof == 0 {
        return 1.0;
    }
    ChiSquared::new(dof as f64)
        .map(|d| d.sf(chisq))
        .unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
