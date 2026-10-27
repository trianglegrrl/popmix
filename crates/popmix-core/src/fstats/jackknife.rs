//! The weighted block jackknife (delete-m jackknife for unequal m).
//!
//! Blocks k = 1..g have lengths m_k (SNPs), n = Σ m_k and h_k = n / m_k. For
//! a statistic that is a weighted mean of per-block values x_k the full-data
//! value is θ = Σ m_k x_k / n and the leave-one-block-out value is
//! θ_(k) = (θ − x_k m_k / n) / (1 − m_k / n). The jackknife estimate is
//! θ_J = g θ − Σ_k (1 − m_k / n) θ_(k) and its variance, through the
//! pseudo-values τ_k = h_k θ − (h_k − 1) θ_(k), is
//! var = (1/g) Σ_k (τ_k − θ_J)² / (h_k − 1).
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/jackknife.md`, after Busing, Meijer and van der Leeden
//! 1999 (Statistics and Computing 9:3, "Delete-m jackknife for unequal m",
//! eqs. 5 and 7) as used for f-statistics by Patterson et al. 2012, and
//! checked against admixtools2's black-box outputs.

use ndarray::Array2;

/// θ = Σ m_k x_k / n over the given blocks, and n.
fn weighted_mean(values: &[f64], lengths: &[f64]) -> (f64, f64) {
    let (s, n) = values
        .iter()
        .zip(lengths)
        .fold((0.0, 0.0), |(s, n), (x, m)| (s + x * m, n + m));
    (s / n, n)
}

/// The blocks whose value is finite: their values and lengths.
fn finite_blocks(values: &[f64], lengths: &[f64]) -> (Vec<f64>, Vec<f64>) {
    values
        .iter()
        .zip(lengths)
        .filter(|(x, _)| x.is_finite())
        .map(|(&x, &m)| (x, m))
        .unzip()
}

/// The scaled deviations u_k = (τ_k − θ_J) / sqrt(h_k − 1) and θ_J, from
/// leave-one-out values, the full-data value and the lengths of the blocks
/// in use (all finite). The variance is (1/g) Σ u_k².
fn deviations(loo: &[f64], lengths: &[f64], tot: f64) -> (Vec<f64>, f64) {
    let g = loo.len() as f64;
    let n: f64 = lengths.iter().sum();
    let est = g * tot
        - loo
            .iter()
            .zip(lengths)
            .map(|(t, m)| (1.0 - m / n) * t)
            .sum::<f64>();
    let u = loo
        .iter()
        .zip(lengths)
        .map(|(t, m)| {
            let h = n / m;
            let tau = h * tot - (h - 1.0) * t;
            (tau - est) / (h - 1.0).sqrt()
        })
        .collect();
    (u, est)
}

/// θ_J and the standard error sqrt((1/g) Σ u_k²).
fn estimate_and_se(loo: &[f64], lengths: &[f64], tot: f64) -> (f64, f64) {
    let (u, est) = deviations(loo, lengths, tot);
    let var = u.iter().map(|x| x * x).sum::<f64>() / u.len() as f64;
    (est, var.sqrt())
}

/// z of an estimate and its standard error: est / se, NaN when the SE is 0
/// or not finite (no information, not "no effect").
pub fn z_score(est: f64, se: f64) -> f64 {
    if se > 0.0 && se.is_finite() {
        est / se
    } else {
        f64::NAN
    }
}

/// Two-sided normal p-value of `z`, as admixtools2's `ztop`; NaN for NaN.
pub fn p_from_z(z: f64) -> f64 {
    libm::erfc(z.abs() / std::f64::consts::SQRT_2)
}

/// Jackknife estimate and standard error of a block mean.
///
/// Blocks whose value is not finite are left out entirely; `None` if none
/// remains.
pub fn block_jackknife(values: &[f64], lengths: &[f64]) -> Option<(f64, f64)> {
    let (x, m) = finite_blocks(values, lengths);
    if x.is_empty() {
        return None;
    }
    let (tot, _) = weighted_mean(&x, &m);
    Some(estimate_and_se(&leave_one_out(&x, &m), &m, tot))
}

/// Leave-one-block-out values θ_(k) of a block mean.
///
/// θ and n use only the finite blocks; a non-finite block's own value stays
/// non-finite.
pub fn leave_one_out(values: &[f64], lengths: &[f64]) -> Vec<f64> {
    let (x, m) = finite_blocks(values, lengths);
    let (tot, n) = weighted_mean(&x, &m);
    values
        .iter()
        .zip(lengths)
        .map(|(x, m)| (tot - x * m / n) / (1.0 - m / n))
        .collect()
}

/// Leave-one-block-out values with every block given a value.
///
/// θ is the mean over the finite blocks; each θ_(k) then takes a non-finite
/// block's value as 0, with m_k / n relative to the total length of all
/// blocks.
pub fn leave_one_out_nafix(values: &[f64], lengths: &[f64]) -> Vec<f64> {
    let (x, m) = finite_blocks(values, lengths);
    let (tot, _) = weighted_mean(&x, &m);
    let n: f64 = lengths.iter().sum();
    values
        .iter()
        .zip(lengths)
        .map(|(&x, m)| {
            let x = if x.is_finite() { x } else { 0.0 };
            (tot - x * m / n) / (1.0 - m / n)
        })
        .collect()
}

/// Jackknife estimate and standard error from given leave-one-out values
/// and a given full-data value `tot`, for statistics that are not block
/// means (such as ratios).
///
/// Blocks with a non-finite θ_(k) are left out; `None` if none remains.
pub fn jackknife_with_total(loo: &[f64], lengths: &[f64], tot: f64) -> Option<(f64, f64)> {
    let (t, m) = finite_blocks(loo, lengths);
    if t.is_empty() {
        return None;
    }
    Some(estimate_and_se(&t, &m, tot))
}

/// Jackknife estimate and standard error from leave-one-out values alone,
/// for block means and ratios of block means whose full-data value is not
/// given (the genotype f-statistics).
///
/// The total θ is the mean of the finite θ_(k) weighted by 1 − m_k / n,
/// n the total length of **all** blocks; then as [`jackknife_with_total`].
/// From `docs/cleanroom/small-helpers.md`, section 2.
///
/// That θ is the full-data value only for an additive statistic (a block
/// mean, or a linear combination of block means such as f3 and f4). For a
/// ratio of block means, such as D, it is the weighted mean of the
/// leave-one-out ratios, which differs from the full-data ratio by a term of
/// order 1/g: this is admixtools2's estimate for such statistics (as its
/// `jack_pairarr_stats` and `qpdstat`), matched for parity, not a bias
/// correction. Give the full-data value to [`jackknife_with_total`] where it
/// is known.
pub fn jackknife_from_loo(loo: &[f64], lengths: &[f64]) -> Option<(f64, f64)> {
    let n: f64 = lengths.iter().sum();
    let (s, w) = loo.iter().zip(lengths).filter(|(t, _)| t.is_finite()).fold(
        (0.0, 0.0),
        |(s, w), (t, m)| {
            let wk = 1.0 - m / n;
            (s + wk * t, w + wk)
        },
    );
    jackknife_with_total(loo, lengths, s / w)
}

/// Jackknife estimates and covariance of several block-mean statistics
/// (`stats[s][k]`).
///
/// Only the blocks where every statistic is finite are used (`None` if
/// fewer than two). The covariance of statistics s and t is
/// (1/g) Σ_k u_sk u_tk with u_sk = (τ_sk − θ_J,s) / sqrt(h_k − 1).
pub fn block_jackknife_cov(stats: &[Vec<f64>], lengths: &[f64]) -> Option<(Vec<f64>, Array2<f64>)> {
    let keep: Vec<usize> = (0..lengths.len())
        .filter(|&k| stats.iter().all(|s| s[k].is_finite()))
        .collect();
    if keep.len() < 2 {
        return None;
    }
    let m: Vec<f64> = keep.iter().map(|&k| lengths[k]).collect();
    let g = keep.len() as f64;

    let (u, ests): (Vec<Vec<f64>>, Vec<f64>) = stats
        .iter()
        .map(|s| {
            let x: Vec<f64> = keep.iter().map(|&k| s[k]).collect();
            let (tot, _) = weighted_mean(&x, &m);
            deviations(&leave_one_out(&x, &m), &m, tot)
        })
        .unzip();

    let p = stats.len();
    let mut cov = Array2::zeros((p, p));
    for a in 0..p {
        for b in a..p {
            let c = u[a].iter().zip(&u[b]).map(|(x, y)| x * y).sum::<f64>() / g;
            cov[[a, b]] = c;
            cov[[b, a]] = c;
        }
    }
    Some((ests, cov))
}

/// Jackknife estimates and covariance of several block-mean statistics
/// from their leave-one-out values alone (`loo[s][k]`), such as linear
/// combinations of the leave-one-out values of other block means.
///
/// Only the blocks where every θ_(k) is finite are used (`None` if fewer
/// than two). Each statistic's θ is the mean of its θ_(k) over those blocks
/// weighted by 1 − m_k / n, n the total length of **all** blocks: for a
/// block mean this is the full-data mean exactly, as in
/// [`jackknife_from_loo`]; everything below holds only for additive
/// statistics (block means and their linear combinations), not ratios. For a block mean θ_J = θ, so the estimate is θ,
/// and τ_k − θ_J = (h_k − 1)(θ − θ_(k)): the covariance of statistics s and
/// t is (1/g) Σ_k u_sk u_tk with u_sk = sqrt(h_k − 1)(θ_s − θ_(k),s), h_k
/// and g over the blocks used. That is [`block_jackknife_cov`]'s covariance,
/// written without the cancellation in g θ − Σ_k (1 − m_k / n) θ_(k) and in
/// h_k θ − (h_k − 1) θ_(k), which loses digits when there are many blocks.
pub fn jackknife_cov_from_loo(
    loo: &[Vec<f64>],
    lengths: &[f64],
) -> Option<(Vec<f64>, Array2<f64>)> {
    let keep: Vec<usize> = (0..lengths.len())
        .filter(|&k| loo.iter().all(|s| s[k].is_finite()))
        .collect();
    if keep.len() < 2 {
        return None;
    }
    let n_all: f64 = lengths.iter().sum();
    let m: Vec<f64> = keep.iter().map(|&k| lengths[k]).collect();
    let n: f64 = m.iter().sum();
    let g = keep.len() as f64;
    let (u, ests): (Vec<Vec<f64>>, Vec<f64>) = loo
        .iter()
        .map(|s| {
            let t: Vec<f64> = keep.iter().map(|&k| s[k]).collect();
            let (sw, w) = t.iter().zip(&m).fold((0.0, 0.0), |(sw, w), (t, mk)| {
                let wk = 1.0 - mk / n_all;
                (sw + wk * t, w + wk)
            });
            let tot = sw / w;
            let u = t
                .iter()
                .zip(&m)
                .map(|(tk, mk)| (n / mk - 1.0).sqrt() * (tot - tk))
                .collect();
            (u, tot)
        })
        .unzip();
    // statistics x blocks; the covariance is U Uᵀ / g
    let u = Array2::from_shape_fn((loo.len(), keep.len()), |(s, k)| u[s][k]);
    Some((ests, u.dot(&u.t()) / g))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn the_total_from_leave_one_out_values_is_exact_only_for_additive_statistics() {
        // G12: θ rebuilt from the θ_(k) is the full-data value for a block
        // mean, but not for a ratio of block means
        let lengths = [10.0, 30.0, 20.0, 40.0];
        let num = [1.0, 3.0, 2.0, 5.0];
        let den = [2.0, 4.0, 7.0, 6.0];
        let (full_num, _) = weighted_mean(&num, &lengths);
        let (full_den, _) = weighted_mean(&den, &lengths);
        let loo_num = leave_one_out(&num, &lengths);
        let (est, _) = jackknife_from_loo(&loo_num, &lengths).unwrap();
        assert_relative_eq!(est, full_num, max_relative = 1e-14);

        let loo_den = leave_one_out(&den, &lengths);
        let loo_ratio: Vec<f64> = loo_num.iter().zip(&loo_den).map(|(a, b)| a / b).collect();
        let (est, _) = jackknife_from_loo(&loo_ratio, &lengths).unwrap();
        let full_ratio = full_num / full_den;
        assert!(
            (est - full_ratio).abs() > 1e-4 * full_ratio,
            "{est} {full_ratio}"
        );
        assert!(
            (est - full_ratio).abs() < 1e-2 * full_ratio,
            "{est} {full_ratio}"
        );
    }

    #[test]
    fn z_is_nan_without_a_standard_error_and_p_is_two_sided() {
        assert_eq!(z_score(1.0, 0.5), 2.0);
        assert!(z_score(1.0, 0.0).is_nan());
        assert!(z_score(1.0, f64::NAN).is_nan());
        // R 4.5: print(2 * pnorm(-1.96), digits = 17)
        let p = p_from_z(1.96);
        assert!((p - 0.049_995_790_296_440_856).abs() < 1e-16, "{p:.17}");
        assert!((p_from_z(-1.96) - p_from_z(1.96)).abs() < 1e-18);
        assert!(p_from_z(f64::NAN).is_nan());
    }

    #[test]
    fn test_equal_blocks_give_the_standard_jackknife() {
        // With equal lengths the weighted jackknife is the delete-one
        // jackknife: var = (g - 1)/g * sum (loo - mean loo)^2.
        let x = [1.0, 2.0, 4.0, 7.0];
        let g = x.len() as f64;
        let mean = x.iter().sum::<f64>() / g;
        let loo: Vec<f64> = x.iter().map(|v| (mean * g - v) / (g - 1.0)).collect();
        let var = (g - 1.0) / g * loo.iter().map(|l| (l - mean).powi(2)).sum::<f64>();
        let (est, se) = block_jackknife(&x, &[5.0; 4]).unwrap();
        assert_relative_eq!(est, mean, max_relative = 1e-14);
        assert_relative_eq!(se, var.sqrt(), max_relative = 1e-14);
    }

    #[test]
    fn test_covariance_diagonal_is_the_single_statistic_variance() {
        let x = vec![0.3, -0.1, 0.8, 0.2, 0.5];
        let y = vec![1.0, 0.7, 0.2, 0.9, 0.4];
        let lengths = [120.0, 80.0, 95.0, 130.0, 60.0];
        let (ests, cov) = block_jackknife_cov(&[x.clone(), y.clone()], &lengths).unwrap();
        for (s, v) in [x, y].iter().enumerate() {
            let (est, se) = block_jackknife(v, &lengths).unwrap();
            assert_relative_eq!(ests[s], est, max_relative = 1e-12);
            assert_relative_eq!(cov[[s, s]], se * se, max_relative = 1e-12);
        }
        assert_eq!(cov[[0, 1]], cov[[1, 0]]);
        assert!(cov[[0, 1]] < 0.0, "x and y move in opposite directions");
    }

    #[test]
    fn test_jackknife_with_total_of_a_block_mean_is_the_plain_jackknife() {
        let x = [0.3, -0.1, 0.8, 0.2, 0.5];
        let lengths = [120.0, 80.0, 95.0, 130.0, 60.0];
        let tot =
            x.iter().zip(&lengths).map(|(v, l)| v * l).sum::<f64>() / lengths.iter().sum::<f64>();
        let loo = leave_one_out(&x, &lengths);
        let (est, se) = jackknife_with_total(&loo, &lengths, tot).unwrap();
        let (est2, se2) = block_jackknife(&x, &lengths).unwrap();
        assert_relative_eq!(est, est2, max_relative = 1e-12);
        assert_relative_eq!(se, se2, max_relative = 1e-12);
    }

    #[test]
    fn test_nafix_counts_missing_blocks_as_zero() {
        let loo = leave_one_out_nafix(&[1.0, f64::NAN, 3.0], &[1.0, 1.0, 2.0]);
        // tot = (1 + 6) / 3 over the finite blocks; rel = l / 4
        let tot = 7.0 / 3.0;
        assert_relative_eq!(loo[1], tot / 0.75, max_relative = 1e-14);
        assert_relative_eq!(loo[2], (tot - 3.0 * 0.5) / 0.5, max_relative = 1e-14);
    }

    #[test]
    fn test_jackknife_from_loo_of_a_block_mean_is_the_plain_jackknife() {
        let x = [0.3, -0.1, f64::NAN, 0.8, 0.2, 0.5];
        let lengths = [120.0, 80.0, 0.0, 95.0, 130.0, 60.0];
        let (est, se) = jackknife_from_loo(&leave_one_out(&x, &lengths), &lengths).unwrap();
        let (est2, se2) = block_jackknife(&x, &lengths).unwrap();
        assert_relative_eq!(est, est2, max_relative = 1e-12);
        assert_relative_eq!(se, se2, max_relative = 1e-12);
    }

    #[test]
    fn test_cov_from_loo_of_block_means_is_the_block_jackknife_cov() {
        let x = vec![0.3, -0.1, 0.8, 0.2, 0.5, f64::NAN];
        let y = vec![1.0, 0.7, 0.2, 0.9, 0.4, f64::NAN];
        let lengths = [120.0, 80.0, 95.0, 130.0, 60.0, 0.0];
        // a linear combination of leave-one-out values is the leave-one-out
        // value of the same combination of block means
        let (lx, ly) = (leave_one_out(&x, &lengths), leave_one_out(&y, &lengths));
        let lz: Vec<f64> = lx.iter().zip(&ly).map(|(a, b)| 2.0 * a - b).collect();
        let z: Vec<f64> = x.iter().zip(&y).map(|(a, b)| 2.0 * a - b).collect();
        let (ests, cov) = jackknife_cov_from_loo(&[lx, ly, lz], &lengths).unwrap();
        let (ests2, cov2) = block_jackknife_cov(&[x, y, z], &lengths).unwrap();
        for s in 0..3 {
            assert_relative_eq!(ests[s], ests2[s], max_relative = 1e-12);
            for t in 0..3 {
                assert_relative_eq!(cov[[s, t]], cov2[[s, t]], max_relative = 1e-12);
            }
        }
        assert!(jackknife_cov_from_loo(&[vec![1.0, f64::NAN]], &[1.0, 1.0]).is_none());
    }

    #[test]
    fn test_non_finite_blocks_are_left_out() {
        let with_nan = block_jackknife(&[1.0, f64::NAN, 3.0], &[2.0, 9.0, 1.0]).unwrap();
        let without = block_jackknife(&[1.0, 3.0], &[2.0, 1.0]).unwrap();
        assert_eq!(with_nan, without);
        assert!(block_jackknife(&[f64::NAN], &[1.0]).is_none());
    }
}
