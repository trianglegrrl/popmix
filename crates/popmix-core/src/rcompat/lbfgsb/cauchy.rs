//! The generalized Cauchy point: the first local minimizer of the quadratic
//! model along the projected steepest-descent path (Byrd et al. 1995,
//! section 4).

use super::*;

/// Max norm of the projected gradient.
pub(super) fn projgr(prob: &Problem, x: &[f64], g: &[f64]) -> f64 {
    let mut sbgnrm: f64 = 0.0;
    for i in 0..x.len() {
        let mut gi = g[i];
        let nbd = prob.nbd[i];
        if nbd != 0 {
            if gi < 0.0 {
                if nbd >= 2 {
                    gi = (x[i] - prob.u[i]).max(gi);
                }
            } else if nbd <= 2 {
                gi = (x[i] - prob.l[i]).min(gi);
            }
        }
        sbgnrm = sbgnrm.max(gi.abs());
    }
    sbgnrm
}

/// Work vectors of the Cauchy point computation, of length 2m: p = W'd,
/// c = W'(xcp - x), and scratch.
pub(super) struct CauchyWork {
    pub(super) p: Vec<f64>,
    pub(super) c: Vec<f64>,
    pub(super) wbp: Vec<f64>,
    pub(super) v: Vec<f64>,
}

/// The generalised Cauchy point `xcp`: the first local minimiser of the
/// quadratic model along the projected steepest descent path x(t) =
/// P(x - t g), examined segment by segment between breakpoints taken from a
/// heap. Marks the variables it fixes at bounds in `iwhere`. Returns a
/// nonzero info on a singular system.
#[allow(clippy::too_many_arguments)]
pub(super) fn cauchy(
    prob: &Problem,
    x: &[f64],
    g: &[f64],
    iorder: &mut [usize],
    iwhere: &mut [i32],
    t: &mut [f64],
    d: &mut [f64],
    xcp: &mut [f64],
    mem: &Memory,
    w: &mut CauchyWork,
    sbgnrm: f64,
) -> usize {
    let n = x.len();
    let (m, col, head, theta) = (mem.m, mem.col, mem.head, mem.theta);
    let (l, u, nbd) = (prob.l, prob.u, prob.nbd);
    if sbgnrm <= 0.0 {
        xcp.copy_from_slice(x);
        return 0;
    }
    let mut bnded = true;
    let mut nfree = n + 1;
    let mut nbreak = 0usize;
    let mut ibkmin = 0usize;
    let mut bkmin = 0.0;
    let col2 = 2 * col;
    let mut f1 = 0.0;
    for v in &mut w.p[..col2] {
        *v = 0.0;
    }
    for i in 1..=n {
        let neggi = -g[i - 1];
        let mut tl = 0.0;
        let mut tu = 0.0;
        let b = nbd[i - 1];
        if iwhere[i - 1] != 3 && iwhere[i - 1] != -1 {
            if b <= 2 {
                tl = x[i - 1] - l[i - 1];
            }
            if b >= 2 {
                tu = u[i - 1] - x[i - 1];
            }
            let xlower = b <= 2 && tl <= 0.0;
            let xupper = b >= 2 && tu <= 0.0;
            iwhere[i - 1] = 0;
            if xlower {
                if neggi <= 0.0 {
                    iwhere[i - 1] = 1;
                }
            } else if xupper {
                if neggi >= 0.0 {
                    iwhere[i - 1] = 2;
                }
            } else if neggi.abs() <= 0.0 {
                iwhere[i - 1] = -3;
            }
        }
        let mut pointr = head;
        if iwhere[i - 1] != 0 && iwhere[i - 1] != -1 {
            d[i - 1] = 0.0;
        } else {
            d[i - 1] = neggi;
            f1 -= neggi * neggi;
            for j in 1..=col {
                w.p[j - 1] += mem.wy(i, pointr) * neggi;
                w.p[col + j - 1] += mem.ws(i, pointr) * neggi;
                pointr = next(pointr, m);
            }
            if b <= 2 && b != 0 && neggi < 0.0 {
                nbreak += 1;
                iorder[nbreak - 1] = i;
                t[nbreak - 1] = tl / (-neggi);
                if nbreak == 1 || t[nbreak - 1] < bkmin {
                    bkmin = t[nbreak - 1];
                    ibkmin = nbreak;
                }
            } else if b >= 2 && neggi > 0.0 {
                nbreak += 1;
                iorder[nbreak - 1] = i;
                t[nbreak - 1] = tu / neggi;
                if nbreak == 1 || t[nbreak - 1] < bkmin {
                    bkmin = t[nbreak - 1];
                    ibkmin = nbreak;
                }
            } else {
                nfree -= 1;
                iorder[nfree - 1] = i;
                if neggi.abs() > 0.0 {
                    bnded = false;
                }
            }
        }
    }
    if theta != 1.0 {
        for v in &mut w.p[col..col2] {
            *v *= theta;
        }
    }
    xcp.copy_from_slice(x);
    if nbreak == 0 && nfree == n + 1 {
        return 0;
    }
    for v in &mut w.c[..col2] {
        *v = 0.0;
    }
    let mut f2 = -theta * f1;
    let f2_org = f2;
    if col > 0 {
        let info = bmv(mem, &w.p, &mut w.v);
        if info != 0 {
            return info;
        }
        f2 -= ddot(&w.v[..col2], &w.p[..col2]);
    }
    let mut dtm = -f1 / f2;
    let mut tsum = 0.0;
    let mut reached_last = false;
    if nbreak > 0 {
        let mut nleft = nbreak;
        let mut iter = 1usize;
        let mut tj = 0.0;
        loop {
            let tj0 = tj;
            let ibp;
            if iter == 1 {
                tj = bkmin;
                ibp = iorder[ibkmin - 1];
            } else {
                if iter == 2 && ibkmin != nbreak {
                    t[ibkmin - 1] = t[nbreak - 1];
                    iorder[ibkmin - 1] = iorder[nbreak - 1];
                }
                hpsolb(nleft, t, iorder, iter - 2);
                tj = t[nleft - 1];
                ibp = iorder[nleft - 1];
            }
            let dt = tj - tj0;
            if dtm < dt {
                break;
            }
            tsum += dt;
            nleft -= 1;
            iter += 1;
            let dibp = d[ibp - 1];
            d[ibp - 1] = 0.0;
            let zibp;
            if dibp > 0.0 {
                zibp = u[ibp - 1] - x[ibp - 1];
                xcp[ibp - 1] = u[ibp - 1];
                iwhere[ibp - 1] = 2;
            } else {
                zibp = l[ibp - 1] - x[ibp - 1];
                xcp[ibp - 1] = l[ibp - 1];
                iwhere[ibp - 1] = 1;
            }
            if nleft == 0 && nbreak == n {
                dtm = dt;
                reached_last = true;
                break;
            }
            let dibp2 = dibp * dibp;
            f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
            f2 -= theta * dibp2;
            if col > 0 {
                for j in 0..col2 {
                    w.c[j] += dt * w.p[j];
                }
                let mut pointr = head;
                for j in 1..=col {
                    w.wbp[j - 1] = mem.wy(ibp, pointr);
                    w.wbp[col + j - 1] = theta * mem.ws(ibp, pointr);
                    pointr = next(pointr, m);
                }
                let info = bmv(mem, &w.wbp, &mut w.v);
                if info != 0 {
                    return info;
                }
                let wmc = ddot(&w.c[..col2], &w.v[..col2]);
                let wmp = ddot(&w.p[..col2], &w.v[..col2]);
                let wmw = ddot(&w.wbp[..col2], &w.v[..col2]);
                for j in 0..col2 {
                    w.p[j] += -dibp * w.wbp[j];
                }
                f1 += dibp * wmc;
                f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
            }
            f2 = (EPSMCH * f2_org).max(f2);
            if nleft > 0 {
                dtm = -f1 / f2;
                continue;
            } else if bnded {
                dtm = 0.0;
            } else {
                dtm = -f1 / f2;
            }
            break;
        }
    }
    if !reached_last {
        if dtm <= 0.0 {
            dtm = 0.0;
        }
        tsum += dtm;
        for i in 0..n {
            xcp[i] += tsum * d[i];
        }
    }
    if col > 0 {
        for j in 0..col2 {
            w.c[j] += dtm * w.p[j];
        }
    }
    0
}

/// Heap of breakpoints: with `iheap` 0 first arrange t(1..n) as a heap,
/// then move its least element to t(n).
pub(super) fn hpsolb(n: usize, t: &mut [f64], iorder: &mut [usize], iheap: usize) {
    if iheap == 0 {
        for k in 2..=n {
            let ddum = t[k - 1];
            let indxin = iorder[k - 1];
            let mut i = k;
            while i > 1 {
                let j = i / 2;
                if ddum < t[j - 1] {
                    t[i - 1] = t[j - 1];
                    iorder[i - 1] = iorder[j - 1];
                    i = j;
                } else {
                    break;
                }
            }
            t[i - 1] = ddum;
            iorder[i - 1] = indxin;
        }
    }
    if n > 1 {
        let mut i = 1;
        let out = t[0];
        let indxou = iorder[0];
        let ddum = t[n - 1];
        let indxin = iorder[n - 1];
        loop {
            let mut j = i + i;
            if j < n {
                if t[j] < t[j - 1] {
                    j += 1;
                }
                if t[j - 1] < ddum {
                    t[i - 1] = t[j - 1];
                    iorder[i - 1] = iorder[j - 1];
                    i = j;
                    continue;
                }
            }
            break;
        }
        t[i - 1] = ddum;
        iorder[i - 1] = indxin;
        t[n - 1] = out;
        iorder[n - 1] = indxou;
    }
}
