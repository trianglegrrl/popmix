//! Weighted block jackknife, as in admixtools2.
//!
//! admixtools2 turns per-block estimates into leave-one-block-out estimates
//! (`est_to_loo`), weighting each block by its length in SNPs, and then
//! computes the jackknife estimate and variance from them
//! (`cpp_jack_vec_stats`, which without a total falls through to
//! `cpp_jack_vec_stats2`). The steps are reproduced literally, rather than
//! simplified algebraically, so results agree with admixtools2 to rounding.
//!
//! [`block_jackknife_cov`] is the same for several statistics at once, giving
//! their covariance (`est_to_loo_nafix` + `jack_pairarr_stats`, as used by
//! qpWave and qpAdm).

use ndarray::Array2;

/// Weighted block jackknife estimate and standard error.
///
/// `values[k]` is the statistic in block k and `lengths[k]` its length in
/// SNPs. Blocks with a non-finite value are left out, as admixtools2 does
/// for NA blocks. Returns `None` if no block is left.
pub fn block_jackknife(values: &[f64], lengths: &[f64]) -> Option<(f64, f64)> {
    assert_eq!(values.len(), lengths.len());
    let (x, bl): (Vec<f64>, Vec<f64>) = values
        .iter()
        .zip(lengths)
        .filter(|(v, _)| v.is_finite())
        .map(|(&v, &l)| (v, l))
        .unzip();
    if x.is_empty() {
        return None;
    }

    // est_to_loo: tot = weighted.mean(x, bl); loo = (tot - x rel) / (1 - rel)
    let n: f64 = bl.iter().sum();
    let tot = x.iter().zip(&bl).map(|(v, l)| v * l).sum::<f64>() / n;
    let loo: Vec<f64> = x
        .iter()
        .zip(&bl)
        .map(|(v, l)| {
            let rel = l / n;
            (tot - v * rel) / (1.0 - rel)
        })
        .collect();

    // cpp_jack_vec_stats2
    let h: Vec<f64> = bl.iter().map(|l| n / l).collect();
    let w_sum: f64 = h.iter().map(|h| 1.0 - 1.0 / h).sum();
    let est = loo
        .iter()
        .zip(&h)
        .map(|(l, h)| l * (1.0 - 1.0 / h))
        .sum::<f64>()
        / w_sum;
    let var = loo
        .iter()
        .zip(&h)
        .map(|(l, h)| (est - l).powi(2) * (h - 1.0))
        .sum::<f64>()
        / loo.len() as f64;
    Some((est, var.sqrt()))
}

/// Weighted block jackknife estimates and covariance of several statistics.
///
/// `stats[s][k]` is statistic s in block k and `lengths[k]` the block's
/// length in SNPs. Blocks where any statistic is non-finite are left out, as
/// admixtools2 `jack_pairarr_stats` does. (admixtools2's `est_to_loo_nafix`
/// would first set missing block values to 0; popmix drops the block
/// instead. The two agree whenever no block is missing, which the default
/// `maxmiss = 0` guarantees.) Returns `None` if fewer than two blocks are left.
pub fn block_jackknife_cov(stats: &[Vec<f64>], lengths: &[f64]) -> Option<(Vec<f64>, Array2<f64>)> {
    let keep: Vec<usize> = (0..lengths.len())
        .filter(|&k| stats.iter().all(|x| x[k].is_finite()))
        .collect();
    if keep.len() < 2 {
        return None;
    }
    let bl: Vec<f64> = keep.iter().map(|&k| lengths[k]).collect();
    let n: f64 = bl.iter().sum();
    let numblocks = bl.len() as f64;
    let h: Vec<f64> = bl.iter().map(|l| n / l).collect();

    let mut ests = Vec::with_capacity(stats.len());
    let mut xtau: Vec<Vec<f64>> = Vec::with_capacity(stats.len());
    for x in stats {
        let x: Vec<f64> = keep.iter().map(|&k| x[k]).collect();
        // est_to_loo_nafix
        let tot0 = x.iter().zip(&bl).map(|(v, l)| v * l).sum::<f64>() / n;
        let loo: Vec<f64> = x
            .iter()
            .zip(&bl)
            .map(|(v, l)| {
                let rel = l / n;
                (tot0 - v * rel) / (1.0 - rel)
            })
            .collect();
        // jack_pairarr_stats
        let w: Vec<f64> = bl.iter().map(|l| 1.0 - l / n).collect();
        let tot = loo.iter().zip(&w).map(|(v, w)| v * w).sum::<f64>() / w.iter().sum::<f64>();
        let est = loo.iter().map(|v| tot - v).sum::<f64>() / numblocks * numblocks
            + loo.iter().zip(&bl).map(|(v, l)| v * l).sum::<f64>() / n;
        xtau.push(
            loo.iter()
                .zip(&h)
                .map(|(v, h)| (est * h - v * (h - 1.0) - tot) / (h - 1.0).sqrt())
                .collect(),
        );
        ests.push(est);
    }

    let m = stats.len();
    let mut cov = Array2::zeros((m, m));
    for i in 0..m {
        for j in i..m {
            let c = xtau[i]
                .iter()
                .zip(&xtau[j])
                .map(|(a, b)| a * b)
                .sum::<f64>()
                / numblocks;
            cov[[i, j]] = c;
            cov[[j, i]] = c;
        }
    }
    Some((ests, cov))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
    fn test_non_finite_blocks_are_left_out() {
        let with_nan = block_jackknife(&[1.0, f64::NAN, 3.0], &[2.0, 9.0, 1.0]).unwrap();
        let without = block_jackknife(&[1.0, 3.0], &[2.0, 1.0]).unwrap();
        assert_eq!(with_nan, without);
        assert!(block_jackknife(&[f64::NAN], &[1.0]).is_none());
    }
}
