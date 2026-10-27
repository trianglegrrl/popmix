//! admixtools2 `qpadm_p`: the fit of one qpAdm model at one rank, without
//! standard errors or `popdrop`, as `qpadm_rotate`, `qpwave_pairs` and
//! `qpadm_multi(full_results = FALSE)` evaluate many models.

use super::f4_matrix::F4Matrix;
use super::fit;
use crate::fstats::F2Blocks;
use crate::{PopmixError, Result};

/// Result of [`qpadm_p`].
#[derive(Debug, Clone)]
pub struct QpAdmP {
    /// Rank of the fit.
    pub f4rank: usize,
    pub dof: usize,
    pub chi_sq: f64,
    /// NaN when `dof` is 0 (admixtools2 NA).
    pub p_value: f64,
    /// Weights of the left populations after the first (with a target, the
    /// sources), summing to 1, when the fit determines them: at rank 0 with
    /// one such population, and at rank n - 1 or more. admixtools2 reports
    /// them in these cases only too (between, its solve is singular).
    pub weights: Option<Vec<f64>>,
    /// Every weight in [0, 1] (true without weights, as admixtools2).
    pub feasible: bool,
}

/// One qpAdm model's fit (admixtools2 `qpadm_p`): the f4 matrix of `left`
/// (with `target`, if any, put first) against `right`, fitted at rank `rnk`
/// (default: one less than the number of left populations other than the
/// target).
pub fn qpadm_p(
    f2_blocks: &F2Blocks,
    left: &[&str],
    right: &[&str],
    target: Option<&str>,
    rnk: Option<usize>,
    fudge: f64,
    constrained: bool,
) -> Result<QpAdmP> {
    let others: Vec<&str> = left
        .iter()
        .copied()
        .filter(|p| Some(*p) != target)
        .collect();
    let mut full: Vec<&str> = target.into_iter().collect();
    full.extend(&others);
    let rnk = match rnk {
        Some(r) => r,
        None => others.len().checked_sub(1).ok_or_else(|| {
            PopmixError::too_few("qpadm_p", "left populations besides the target", 1, 0)
        })?,
    };
    let f4 = F4Matrix::from_f2_blocks(f2_blocks, &full, right)?;
    let xmat = fit::to_dmatrix(&f4.values);
    let qinv = fit::fudged_inverse(&f4.covariance, fudge)?;
    let rank_fit = fit::fit_rank(&xmat, &qinv, rnk, fudge)?;
    let (chi_sq, dof) = fit::evaluate_fit(&xmat, &qinv, &rank_fit, rnk);
    let p_value = fit::chisq_p_value(chi_sq, dof);
    let nr = xmat.nrows();
    let weights = if rnk == 0 {
        (nr == 1).then(|| vec![1.0])
    } else if rnk + 1 >= nr {
        Some(fit::weights_from_fit(&rank_fit, rnk, constrained)?)
    } else {
        None
    };
    let feasible = weights
        .as_ref()
        .is_none_or(|w| w.iter().all(|v| (0.0..=1.0).contains(v)));
    Ok(QpAdmP {
        f4rank: rnk,
        dof,
        chi_sq,
        p_value,
        weights,
        feasible,
    })
}
