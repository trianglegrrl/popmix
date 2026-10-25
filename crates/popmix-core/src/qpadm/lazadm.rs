//! admixtools2 `lazadm`: admixture weights of a target from left
//! populations by least squares on f4-statistics with the right populations,
//! without qpAdm's covariance weighting, and their block jackknife.
//!
//! For every triple of right populations (r1, r2, r3) with r1 not r2 or r3
//! and r2 < r3 (in population order), the target's f4(target, r1; r2, r3)
//! is regressed on the left populations' f4(left, r1; r2, r3), from the f2
//! blocks with each block left out in turn (`est_to_loo_nafix`: missing
//! blocks count as 0) and from all blocks. The weights are fitted on each
//! (on the simplex when constrained, as quadprog does), scaled to sum to 1
//! when their sum is positive, and jackknifed (`jack_mat_stats` with the
//! all-block weights as the total).

use nalgebra::{DMatrix, DVector};

use super::fit;
use crate::fstats::F2Blocks;
use crate::{PopmixError, Result};

/// Result of [`lazadm`]: one weight per left population.
#[derive(Debug, Clone)]
pub struct LazAdmResult {
    pub left: Vec<String>,
    pub weights: Vec<f64>,
    pub std_errors: Vec<f64>,
}

/// admixtools2 `lazadm(f2_blocks, left, right, target, constrained)`.
pub fn lazadm(
    f2_blocks: &F2Blocks,
    target: &str,
    left: &[&str],
    right: &[&str],
    constrained: bool,
) -> Result<LazAdmResult> {
    let pops: Vec<&str> = std::iter::once(target)
        .chain(left.iter().copied())
        .chain(right.iter().copied())
        .collect();
    for (i, p) in pops.iter().enumerate() {
        if pops[..i].contains(p) {
            return Err(PopmixError::Validation(format!(
                "lazadm: population {p:?} given twice"
            )));
        }
    }
    let n = pops.len();
    let nb = f2_blocks.n_blocks;
    let lengths: Vec<f64> = f2_blocks.block_lengths.iter().map(|&l| l as f64).collect();
    let total_len: f64 = lengths.iter().sum();
    // block f2 means, [i][j][b]
    let mut blocks = vec![vec![vec![0.0; nb]; n]; n];
    for i in 0..n {
        for j in 0..n {
            let v = f2_blocks
                .get_block_f2_values(pops[i], pops[j])
                .ok_or_else(|| {
                    PopmixError::Validation(format!(
                        "lazadm: no f2 for {} and {}",
                        pops[i], pops[j]
                    ))
                })?;
            blocks[i][j] = v.to_vec();
        }
    }
    // est_to_loo_nafix, then the all-block weighted mean as the last block
    let mut loo = vec![vec![vec![0.0; nb + 1]; n]; n];
    for i in 0..n {
        for j in 0..n {
            let x = &blocks[i][j];
            let (mut s, mut w) = (0.0, 0.0);
            for b in 0..nb {
                if !x[b].is_nan() {
                    s += x[b] * lengths[b];
                    w += lengths[b];
                }
            }
            let tot = s / w;
            for b in 0..nb {
                let rel = lengths[b] / total_len;
                let v = if x[b].is_nan() { 0.0 } else { x[b] };
                loo[i][j][b] = (tot - v * rel) / (1.0 - rel);
            }
            loo[i][j][nb] = tot;
        }
    }
    // expand.grid(r, r, r): the first index varies fastest
    let r: Vec<usize> = (1 + left.len()..n).collect();
    let mut triples = Vec::new();
    for &v3 in &r {
        for &v2 in &r {
            for &v1 in &r {
                if v1 != v2 && v1 != v3 && v2 < v3 {
                    triples.push((v1, v2, v3));
                }
            }
        }
    }
    let nl = left.len();
    let mut wmat = vec![vec![0.0; nb + 1]; nl];
    for b in 0..=nb {
        let l = |i: usize, j: usize| loo[i][j][b];
        let y = DVector::from_iterator(
            triples.len(),
            triples
                .iter()
                .map(|&(p2, p3, p4)| l(0, p4) + l(p2, p3) - l(0, p3) - l(p2, p4)),
        );
        let x = DMatrix::from_fn(triples.len(), nl, |c, k| {
            let (p2, p3, p4) = triples[c];
            l(p4, 1 + k) + l(p2, p3) - l(p3, 1 + k) - l(p2, p4)
        });
        let lhs = x.transpose() * &y;
        let rhs = x.transpose() * &x;
        let mut w = if constrained {
            fit::simplex_qp(&rhs, &lhs)?
        } else {
            rhs.lu()
                .solve(&lhs)
                .ok_or_else(|| PopmixError::Numerical("lazadm: singular system".to_string()))?
        };
        let sum = w.sum();
        if sum > 1e-10 {
            w /= sum;
        }
        for k in 0..nl {
            wmat[k][b] = w[k];
        }
    }
    // jack_mat_stats(loo, block_lengths, tot)
    let nbf = nb as f64;
    let (mut weights, mut std_errors) = (Vec::new(), Vec::new());
    for row in &wmat {
        let tot = row[nb];
        let loo = &row[..nb];
        // rowMeans(tot - loo) * numblocks + weighted_row_means(loo, block_lengths)
        let est = loo.iter().map(|&v| tot - v).sum::<f64>() / nbf * nbf
            + loo.iter().zip(&lengths).map(|(v, l)| v * l).sum::<f64>() / total_len;
        let var = loo
            .iter()
            .zip(&lengths)
            .map(|(&v, &bl)| {
                let h = total_len / bl;
                let tau = h * tot - (h - 1.0) * v;
                ((tau - est) / (h - 1.0).sqrt()).powi(2)
            })
            .sum::<f64>()
            / nbf;
        weights.push(est);
        std_errors.push(var.sqrt());
    }
    Ok(LazAdmResult {
        left: left.iter().map(|s| s.to_string()).collect(),
        weights,
        std_errors,
    })
}
