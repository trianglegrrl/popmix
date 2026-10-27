//! Subspace minimization over the free variables at the Cauchy point (Byrd
//! et al. 1995, section 5), truncated at the first bound, as version 2.3.

use super::*;

/// The free and active sets at the Cauchy point.
pub(super) struct FreeSet {
    pub(super) nfree: usize,
    /// free variables first (1..nfree), then the active ones
    pub(super) index: Vec<usize>,
    /// variables entering the free set first (1..nenter), those leaving
    /// last (ileave..n)
    pub(super) indx2: Vec<usize>,
    pub(super) nenter: usize,
    pub(super) ileave: usize,
}

/// Update the free set from `iwhere`; returns whether the subspace matrix
/// must be formed again.
pub(super) fn freev(
    fs: &mut FreeSet,
    iwhere: &[i32],
    updatd: bool,
    cnstnd: bool,
    iter: usize,
) -> bool {
    let n = iwhere.len();
    fs.nenter = 0;
    fs.ileave = n + 1;
    if iter > 0 && cnstnd {
        for i in 1..=fs.nfree {
            let k = fs.index[i - 1];
            if iwhere[k - 1] > 0 {
                fs.ileave -= 1;
                fs.indx2[fs.ileave - 1] = k;
            }
        }
        for i in (fs.nfree + 1)..=n {
            let k = fs.index[i - 1];
            if iwhere[k - 1] <= 0 {
                fs.nenter += 1;
                fs.indx2[fs.nenter - 1] = k;
            }
        }
    }
    let wrk = fs.ileave < n + 1 || fs.nenter > 0 || updatd;
    fs.nfree = 0;
    let mut iact = n + 1;
    for i in 1..=n {
        if iwhere[i - 1] <= 0 {
            fs.nfree += 1;
            fs.index[fs.nfree - 1] = i;
        } else {
            iact -= 1;
            fs.index[iact - 1] = i;
        }
    }
    wrk
}

/// The reduced gradient of the model at the Cauchy point z over the free
/// variables: r = -Z'(B(z - x) + g). Returns a nonzero info on failure.
#[allow(clippy::too_many_arguments)]
pub(super) fn cmprlb(
    x: &[f64],
    g: &[f64],
    z: &[f64],
    r: &mut [f64],
    mem: &Memory,
    w: &mut CauchyWork,
    fs: &FreeSet,
    cnstnd: bool,
) -> i32 {
    let (m, col, head, theta) = (mem.m, mem.col, mem.head, mem.theta);
    if !cnstnd && col > 0 {
        for i in 0..x.len() {
            r[i] = -g[i];
        }
        return 0;
    }
    for i in 1..=fs.nfree {
        let k = fs.index[i - 1];
        r[i - 1] = -theta * (z[k - 1] - x[k - 1]) - g[k - 1];
    }
    // M c, into the p vector
    if bmv(mem, &w.c, &mut w.p) != 0 {
        return -8;
    }
    let mut pointr = head;
    for j in 1..=col {
        let a1 = w.p[j - 1];
        let a2 = theta * w.p[col + j - 1];
        for i in 1..=fs.nfree {
            let k = fs.index[i - 1];
            r[i - 1] = r[i - 1] + mem.wy(k, pointr) * a1 + mem.ws(k, pointr) * a2;
        }
        pointr = next(pointr, m);
    }
    0
}

/// Subspace minimisation (version 2.3): the Newton direction of the model
/// over the free variables, d = (1/theta) r + (1/theta^2) Z'W K^-1 W'Z r,
/// then the step from the Cauchy point `x` along it, truncated at the first
/// bound it meets. `d` holds r on entry; `x` is the Cauchy point on entry
/// and the subspace minimiser on exit. Returns a nonzero info if K is
/// singular.
pub(super) fn subsm(
    prob: &Problem,
    fs: &FreeSet,
    x: &mut [f64],
    d: &mut [f64],
    wv: &mut [f64],
    mem: &Memory,
) -> i32 {
    let (m, col, head, theta) = (mem.m, mem.col, mem.head, mem.theta);
    let nsub = fs.nfree;
    let ind = &fs.index;
    if nsub == 0 {
        return 0;
    }
    // wv = W'Z d
    let mut pointr = head;
    for i in 1..=col {
        let mut temp1 = 0.0;
        let mut temp2 = 0.0;
        for j in 1..=nsub {
            let k = ind[j - 1];
            temp1 += mem.wy(k, pointr) * d[j - 1];
            temp2 += mem.ws(k, pointr) * d[j - 1];
        }
        wv[i - 1] = temp1;
        wv[col + i - 1] = theta * temp2;
        pointr = next(pointr, m);
    }
    // wv = K^-1 wv
    let m2 = 2 * m;
    let col2 = 2 * col;
    if dtrsl(&mem.wn, 0, m2, col2, &mut wv[..col2], true) != 0 {
        return 1;
    }
    for v in &mut wv[..col] {
        *v = -*v;
    }
    if dtrsl(&mem.wn, 0, m2, col2, &mut wv[..col2], false) != 0 {
        return 1;
    }
    // d = (1/theta) d + (1/theta^2) Z'W wv
    let mut pointr = head;
    for jy in 1..=col {
        let js = col + jy;
        for i in 1..=nsub {
            let k = ind[i - 1];
            d[i - 1] =
                d[i - 1] + mem.wy(k, pointr) * wv[jy - 1] / theta + mem.ws(k, pointr) * wv[js - 1];
        }
        pointr = next(pointr, m);
    }
    let s = 1.0 / theta;
    for v in &mut d[..nsub] {
        *v *= s;
    }

    // the largest step alpha in [0, 1] keeping the free variables in their
    // bounds; the variable that limits it is put on its bound
    let (l, u, nbd) = (prob.l, prob.u, prob.nbd);
    let mut alpha = 1.0;
    let mut temp1 = alpha;
    let mut ibd = 0;
    for i in 1..=nsub {
        let k = ind[i - 1];
        let dk = d[i - 1];
        let b = nbd[k - 1];
        if b != 0 {
            if dk < 0.0 && b <= 2 {
                let temp2 = l[k - 1] - x[k - 1];
                if temp2 >= 0.0 {
                    temp1 = 0.0;
                } else if dk * alpha < temp2 {
                    temp1 = temp2 / dk;
                }
            } else if dk > 0.0 && b >= 2 {
                let temp2 = u[k - 1] - x[k - 1];
                if temp2 <= 0.0 {
                    temp1 = 0.0;
                } else if dk * alpha > temp2 {
                    temp1 = temp2 / dk;
                }
            }
            if temp1 < alpha {
                alpha = temp1;
                ibd = i;
            }
        }
    }
    if alpha < 1.0 {
        let dk = d[ibd - 1];
        let k = ind[ibd - 1];
        if dk > 0.0 {
            x[k - 1] = u[k - 1];
            d[ibd - 1] = 0.0;
        } else if dk < 0.0 {
            x[k - 1] = l[k - 1];
            d[ibd - 1] = 0.0;
        }
    }
    for i in 1..=nsub {
        let k = ind[i - 1];
        x[k - 1] += alpha * d[i - 1];
    }
    0
}
