//! R's `optim(method = "L-BFGS-B")`.
//!
//! L-BFGS-B as published by Byrd, Lu, Nocedal and Zhu (1995, "A limited
//! memory algorithm for bound constrained optimization", SIAM J. Sci.
//! Comput. 16:1190), in the form R uses (its version 2.x): each iteration
//! finds the generalized Cauchy point along the projected steepest descent
//! path of the limited-memory quadratic model B = theta I - W M W', minimizes
//! the model over the variables still free there (the direct primal method,
//! backtracked into the box), and searches along the result with the line
//! search of Moré and Thuente (1994, ACM TOMS 20:286; MINPACK-2 `dcsrch`,
//! ftol 1e-3, gtol 0.9, xtol 0.1). It stops when the projected gradient is
//! at most `pgtol` or the relative reduction of f is at most
//! `factr * epsilon`.
//!
//! Around it, as R's `optim`: the objective is divided by `fnscale`; without
//! a gradient, each gradient takes central differences with step `ndeps`,
//! shortened to stay within the bounds; and `maxit` counts iterations, the
//! run stopping (convergence code 1) when one more would exceed it.

use std::collections::VecDeque;

use crate::{PopmixError, Result};

const EPSMCH: f64 = f64::EPSILON;
const FTOL: f64 = 1e-3;
const GTOL: f64 = 0.9;
const XTOL: f64 = 0.1;
const BIG: f64 = 1e10;

/// Controls of R's `optim` that L-BFGS-B uses; defaults are R's.
#[derive(Debug, Clone)]
pub struct LbfgsbOptions {
    /// Number of correction pairs kept (`lmm`, default 5).
    pub lmm: usize,
    /// Relative reduction tolerance, in units of machine epsilon (default 1e7).
    pub factr: f64,
    /// Projected gradient tolerance (default 0: never the reason to stop
    /// unless the projected gradient vanishes).
    pub pgtol: f64,
    /// Iteration limit (default 100).
    pub maxit: usize,
    /// The objective is divided by this (default 1).
    pub fnscale: f64,
    /// Finite-difference step for each variable (default 1e-3).
    pub ndeps: f64,
}

impl Default for LbfgsbOptions {
    fn default() -> Self {
        Self {
            lmm: 5,
            factr: 1e7,
            pgtol: 0.0,
            maxit: 100,
            fnscale: 1.0,
            ndeps: 1e-3,
        }
    }
}

/// What R's `optim` returns.
#[derive(Debug, Clone)]
pub struct LbfgsbResult {
    pub par: Vec<f64>,
    /// The objective at `par` (not divided by `fnscale`).
    pub value: f64,
    /// Objective-and-gradient evaluations (R's `counts`, both entries).
    pub evaluations: usize,
    /// 0 converged, 1 iteration limit, 51 warning, 52 error (R's codes).
    pub convergence: i32,
    pub message: String,
}

/// Bound type of a variable: 0 none, 1 lower, 2 both, 3 upper (`nbd`).
fn bound_type(l: f64, u: f64) -> u8 {
    match (l.is_finite(), u.is_finite()) {
        (true, true) => 2,
        (true, false) => 1,
        (false, true) => 3,
        (false, false) => 0,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| s + x * y)
}

/// LINPACK `dpofa` on the leading n x n block of `a` at offset `o`: the
/// upper Cholesky factor R (A = R'R) in place. Err(j) when not positive
/// definite at column j.
fn dpofa(a: &mut [Vec<f64>], o: usize, n: usize) -> std::result::Result<(), usize> {
    for j in 0..n {
        let mut s = 0.0;
        for k in 0..j {
            let mut t = a[o + k][o + j];
            for i in 0..k {
                t -= a[o + i][o + k] * a[o + i][o + j];
            }
            t /= a[o + k][o + k];
            a[o + k][o + j] = t;
            s += t * t;
        }
        let s = a[o + j][o + j] - s;
        if s <= 0.0 {
            return Err(j + 1);
        }
        a[o + j][o + j] = s.sqrt();
    }
    Ok(())
}

/// LINPACK `dtrsl` with an upper triangular T (n x n at offset `o`):
/// `trans` solves T'x = b, otherwise Tx = b; b is overwritten.
fn dtrsl(
    t: &[Vec<f64>],
    o: usize,
    n: usize,
    b: &mut [f64],
    trans: bool,
) -> std::result::Result<(), usize> {
    if let Some(j) = (0..n).find(|&j| t[o + j][o + j] == 0.0) {
        return Err(j + 1);
    }
    if n == 0 {
        return Ok(());
    }
    if trans {
        b[0] /= t[o][o];
        for j in 1..n {
            let mut s = 0.0;
            for i in 0..j {
                s += t[o + i][o + j] * b[i];
            }
            b[j] = (b[j] - s) / t[o + j][o + j];
        }
    } else {
        b[n - 1] /= t[o + n - 1][o + n - 1];
        for j in (0..n - 1).rev() {
            let temp = -b[j + 1];
            for i in 0..=j {
                b[i] += temp * t[o + i][o + j + 1];
            }
            b[j] /= t[o + j][o + j];
        }
    }
    Ok(())
}

/// The limited-memory BFGS matrix B = theta I - W M W', W = [Y, theta S].
struct Memory {
    m: usize,
    s: VecDeque<Vec<f64>>,
    y: VecDeque<Vec<f64>>,
    /// S'Y: entry (i, j) is s_i'y_j; only the lower triangle is kept
    sy: Vec<Vec<f64>>,
    /// S'S: upper triangle
    ss: Vec<Vec<f64>>,
    /// upper Cholesky factor of theta S'S + L D^-1 L'
    wt: Vec<Vec<f64>>,
    theta: f64,
}

impl Memory {
    fn new(m: usize) -> Self {
        Self {
            m,
            s: VecDeque::new(),
            y: VecDeque::new(),
            sy: vec![vec![0.0; m]; m],
            ss: vec![vec![0.0; m]; m],
            wt: vec![vec![0.0; m]; m],
            theta: 1.0,
        }
    }

    fn col(&self) -> usize {
        self.s.len()
    }

    fn reset(&mut self) {
        self.s.clear();
        self.y.clear();
        self.theta = 1.0;
    }

    /// Add the pair (s, y) (`matupd`): s = stp d, with d'd = `dtd`, s'y = `dr`
    /// and y'y = `rr`.
    fn update(&mut self, s: Vec<f64>, y: Vec<f64>, stp: f64, dtd: f64, dr: f64, rr: f64) {
        self.theta = rr / dr;
        if self.col() == self.m {
            self.s.pop_front();
            self.y.pop_front();
            for i in 0..self.m - 1 {
                for j in 0..self.m - 1 {
                    self.sy[i][j] = self.sy[i + 1][j + 1];
                    self.ss[i][j] = self.ss[i + 1][j + 1];
                }
            }
        }
        let c = self.col();
        for j in 0..c {
            self.sy[c][j] = dot(&s, &self.y[j]);
            self.ss[j][c] = dot(&self.s[j], &s);
        }
        self.ss[c][c] = if stp == 1.0 { dtd } else { stp * stp * dtd };
        self.sy[c][c] = dr;
        self.s.push_back(s);
        self.y.push_back(y);
    }

    /// Factor theta S'S + L D^-1 L' (`formt`).
    fn form_t(&mut self) -> std::result::Result<(), usize> {
        let col = self.col();
        for j in 0..col {
            self.wt[0][j] = self.theta * self.ss[0][j];
        }
        for i in 1..col {
            for j in i..col {
                let mut ddum = 0.0;
                for k in 0..i {
                    ddum += self.sy[i][k] * self.sy[j][k] / self.sy[k][k];
                }
                self.wt[i][j] = ddum + self.theta * self.ss[i][j];
            }
        }
        dpofa(&mut self.wt, 0, col)
    }

    /// The product M v for a vector of length 2 col (`bmv`).
    fn bmv(&self, v: &[f64]) -> std::result::Result<Vec<f64>, usize> {
        let col = self.col();
        let mut p = vec![0.0; 2 * col];
        if col == 0 {
            return Ok(p);
        }
        let sy = &self.sy;
        p[col] = v[col];
        for i in 1..col {
            let mut sum = 0.0;
            for k in 0..i {
                sum += sy[i][k] * v[k] / sy[k][k];
            }
            p[col + i] = v[col + i] + sum;
        }
        dtrsl(&self.wt, 0, col, &mut p[col..], true)?;
        for i in 0..col {
            p[i] = v[i] / sy[i][i].sqrt();
        }
        dtrsl(&self.wt, 0, col, &mut p[col..], false)?;
        for i in 0..col {
            p[i] = -p[i] / sy[i][i].sqrt();
        }
        for i in 0..col {
            let mut sum = 0.0;
            for k in i + 1..col {
                sum += sy[k][i] * p[col + k] / sy[i][i];
            }
            p[i] += sum;
        }
        Ok(p)
    }
}

struct Problem<'a> {
    l: &'a [f64],
    u: &'a [f64],
    nbd: Vec<u8>,
}

impl Problem<'_> {
    /// Infinity norm of the projected gradient (`projgr`).
    fn projected_gradient_norm(&self, x: &[f64], g: &[f64]) -> f64 {
        let mut norm: f64 = 0.0;
        for i in 0..x.len() {
            let mut gi = g[i];
            if self.nbd[i] != 0 {
                if gi < 0.0 {
                    if self.nbd[i] >= 2 {
                        gi = gi.max(x[i] - self.u[i]);
                    }
                } else if self.nbd[i] <= 2 {
                    gi = gi.min(x[i] - self.l[i]);
                }
            }
            norm = norm.max(gi.abs());
        }
        norm
    }
}

/// Heap of breakpoints (`hpsolb`): with `build`, make t[0..n] a heap; then
/// move its least element to t[n-1] and re-heap t[0..n-1].
fn hpsolb(n: usize, t: &mut [f64], order: &mut [usize], build: bool) {
    if build {
        for k in 1..n {
            let ddum = t[k];
            let indxin = order[k];
            let mut i = k + 1;
            while i > 1 {
                let j = i / 2;
                if ddum < t[j - 1] {
                    t[i - 1] = t[j - 1];
                    order[i - 1] = order[j - 1];
                    i = j;
                } else {
                    break;
                }
            }
            t[i - 1] = ddum;
            order[i - 1] = indxin;
        }
    }
    if n > 1 {
        let mut i = 1;
        let out = t[0];
        let indxou = order[0];
        let ddum = t[n - 1];
        let indxin = order[n - 1];
        loop {
            let mut j = i + i;
            if j > n - 1 {
                break;
            }
            if j < n - 1 && t[j] < t[j - 1] {
                j += 1;
            }
            if t[j - 1] < ddum {
                t[i - 1] = t[j - 1];
                order[i - 1] = order[j - 1];
                i = j;
            } else {
                break;
            }
        }
        t[i - 1] = ddum;
        order[i - 1] = indxin;
        t[n - 1] = out;
        order[n - 1] = indxou;
    }
}

/// Where each variable is (`iwhere`): -1 unbounded, 0 free, 1 at its lower
/// bound, 2 at its upper bound, 3 fixed (l = u), -3 free with zero gradient.
type Where = Vec<i8>;

/// The generalized Cauchy point (`cauchy`): xcp, and c = W'(xcp - x).
fn cauchy(
    pr: &Problem,
    x: &[f64],
    g: &[f64],
    iwhere: &mut Where,
    mem: &Memory,
    sbgnrm: f64,
) -> std::result::Result<(Vec<f64>, Vec<f64>), usize> {
    let n = x.len();
    let col = mem.col();
    let col2 = 2 * col;
    let theta = mem.theta;
    let (l, u, nbd) = (pr.l, pr.u, &pr.nbd);
    if sbgnrm <= 0.0 {
        return Ok((x.to_vec(), vec![0.0; col2]));
    }
    let mut bnded = true;
    let mut nfree = n; // free variables fill order[nfree..n], from the end
    let mut nbreak = 0;
    let mut ibkmin = 0;
    let mut bkmin = 0.0;
    let mut f1 = 0.0;
    let mut p = vec![0.0; col2];
    let mut c = vec![0.0; col2];
    let mut d = vec![0.0; n];
    let mut t = vec![0.0; n];
    let mut order = vec![0usize; n];
    for i in 0..n {
        let neggi = -g[i];
        let mut tl = 0.0;
        let mut tu = 0.0;
        if iwhere[i] != 3 && iwhere[i] != -1 {
            if nbd[i] <= 2 {
                tl = x[i] - l[i];
            }
            if nbd[i] >= 2 {
                tu = u[i] - x[i];
            }
            let xlower = nbd[i] <= 2 && tl <= 0.0;
            let xupper = nbd[i] >= 2 && tu <= 0.0;
            iwhere[i] = 0;
            if xlower {
                if neggi <= 0.0 {
                    iwhere[i] = 1;
                }
            } else if xupper {
                if neggi >= 0.0 {
                    iwhere[i] = 2;
                }
            } else if neggi.abs() <= 0.0 {
                iwhere[i] = -3;
            }
        }
        if iwhere[i] != 0 && iwhere[i] != -1 {
            d[i] = 0.0;
        } else {
            d[i] = neggi;
            f1 -= neggi * neggi;
            for j in 0..col {
                p[j] += mem.y[j][i] * neggi;
                p[col + j] += mem.s[j][i] * neggi;
            }
            if nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.0 {
                order[nbreak] = i;
                t[nbreak] = tl / (-neggi);
                if nbreak == 0 || t[nbreak] < bkmin {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
                nbreak += 1;
            } else if nbd[i] >= 2 && neggi > 0.0 {
                order[nbreak] = i;
                t[nbreak] = tu / neggi;
                if nbreak == 0 || t[nbreak] < bkmin {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
                nbreak += 1;
            } else {
                nfree -= 1;
                order[nfree] = i;
                if neggi.abs() > 0.0 {
                    bnded = false;
                }
            }
        }
    }
    if theta != 1.0 {
        for v in p[col..].iter_mut() {
            *v *= theta;
        }
    }
    let mut xcp = x.to_vec();
    if nbreak == 0 && nfree == n {
        return Ok((xcp, c));
    }
    let mut f2 = -theta * f1;
    let f2_org = f2;
    if col > 0 {
        let v = mem.bmv(&p)?;
        f2 -= dot(&v, &p);
    }
    let mut dtm = -f1 / f2;
    let mut tsum = 0.0;
    let mut all_fixed = false;
    if nbreak > 0 {
        let mut nleft = nbreak;
        let mut iter = 1;
        let mut tj = 0.0;
        loop {
            let tj0 = tj;
            let ibp;
            if iter == 1 {
                tj = bkmin;
                ibp = order[ibkmin];
            } else {
                if iter == 2 && ibkmin != nbreak - 1 {
                    t[ibkmin] = t[nbreak - 1];
                    order[ibkmin] = order[nbreak - 1];
                }
                hpsolb(nleft, &mut t, &mut order, iter == 2);
                tj = t[nleft - 1];
                ibp = order[nleft - 1];
            }
            let dt = tj - tj0;
            if dtm < dt {
                break;
            }
            tsum += dt;
            nleft -= 1;
            iter += 1;
            let dibp = d[ibp];
            d[ibp] = 0.0;
            let zibp;
            if dibp > 0.0 {
                zibp = u[ibp] - x[ibp];
                xcp[ibp] = u[ibp];
                iwhere[ibp] = 2;
            } else {
                zibp = l[ibp] - x[ibp];
                xcp[ibp] = l[ibp];
                iwhere[ibp] = 1;
            }
            if nleft == 0 && nbreak == n {
                dtm = dt;
                all_fixed = true;
                break;
            }
            let dibp2 = dibp * dibp;
            f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
            f2 -= theta * dibp2;
            if col > 0 {
                for j in 0..col2 {
                    c[j] += dt * p[j];
                }
                let mut wbp = vec![0.0; col2];
                for j in 0..col {
                    wbp[j] = mem.y[j][ibp];
                    wbp[col + j] = theta * mem.s[j][ibp];
                }
                let v = mem.bmv(&wbp)?;
                let wmc = dot(&c, &v);
                let wmp = dot(&p, &v);
                let wmw = dot(&wbp, &v);
                for j in 0..col2 {
                    p[j] -= dibp * wbp[j];
                }
                f1 += dibp * wmc;
                f2 = f2 + dibp * 2.0 * wmp - dibp2 * wmw;
            }
            f2 = f2.max(EPSMCH * f2_org);
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
    if !all_fixed {
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
            c[j] += dtm * p[j];
        }
    }
    Ok((xcp, c))
}

/// Minimize the model over the variables free at the Cauchy point z
/// (`formk`, `cmprlb`, `subsm`): z is moved along the subspace Newton
/// direction, backtracked to stay within the bounds.
#[allow(clippy::too_many_arguments)]
fn subspace_min(
    pr: &Problem,
    x: &[f64],
    g: &[f64],
    z: &mut [f64],
    c: &[f64],
    free: &[usize],
    active: &[usize],
    mem: &Memory,
) -> std::result::Result<(), usize> {
    let col = mem.col();
    let col2 = 2 * col;
    let theta = mem.theta;
    let (s, y) = (&mem.s, &mem.y);

    // K = [-D - Y'ZZ'Y/theta, L_a' - R_z'; L_a - R_z, theta S'AA'S], factored
    // as J'EJ with E = diag(-I, I), J upper triangular, stored in wn
    let mut wn = vec![vec![0.0; col2]; col2];
    let over =
        |a: &Vec<f64>, b: &Vec<f64>, set: &[usize]| set.iter().fold(0.0, |t, &k| t + a[k] * b[k]);
    for iy in 0..col {
        for jy in 0..=iy {
            wn[jy][iy] = over(&y[iy], &y[jy], free) / theta;
            wn[col + jy][col + iy] = over(&s[iy], &s[jy], active) * theta;
        }
        // (1,2) block, -L_a' + R_z': s_iy'y_jy over the free variables when
        // iy <= jy, minus it over the active ones when iy > jy
        for jy in 0..iy {
            wn[jy][col + iy] = -over(&s[iy], &y[jy], active);
        }
        for jy in iy..col {
            wn[jy][col + iy] = over(&s[iy], &y[jy], free);
        }
        wn[iy][iy] += mem.sy[iy][iy];
    }
    dpofa(&mut wn, 0, col)?;
    for js in col..col2 {
        let mut b: Vec<f64> = (0..col).map(|i| wn[i][js]).collect();
        dtrsl(&wn, 0, col, &mut b, true)?;
        for i in 0..col {
            wn[i][js] = b[i];
        }
    }
    for is in col..col2 {
        for js in is..col2 {
            let sum = wn[..col].iter().fold(0.0, |t, row| t + row[is] * row[js]);
            wn[is][js] += sum;
        }
    }
    dpofa(&mut wn, col, col)?;

    // r = -Z'(B(z - x) + g)
    let v = mem.bmv(c)?;
    let mut d: Vec<f64> = free
        .iter()
        .map(|&k| -theta * (z[k] - x[k]) - g[k])
        .collect();
    for j in 0..col {
        let a1 = v[j];
        let a2 = theta * v[col + j];
        for (i, &k) in free.iter().enumerate() {
            d[i] = d[i] + y[j][k] * a1 + s[j][k] * a2;
        }
    }

    // the Newton direction (1/theta) r + (1/theta^2) Z'W K^-1 W'Z r
    let mut wv = vec![0.0; col2];
    for i in 0..col {
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        for (j, &k) in free.iter().enumerate() {
            t1 += y[i][k] * d[j];
            t2 += s[i][k] * d[j];
        }
        wv[i] = t1;
        wv[col + i] = theta * t2;
    }
    let full: Vec<Vec<f64>> = wn;
    dtrsl(&full, 0, col2, &mut wv, true)?;
    for w in wv[..col].iter_mut() {
        *w = -*w;
    }
    dtrsl(&full, 0, col2, &mut wv, false)?;
    for jy in 0..col {
        let js = col + jy;
        for (i, &k) in free.iter().enumerate() {
            d[i] = d[i] + y[jy][k] * wv[jy] / theta + s[jy][k] * wv[js];
        }
    }
    for di in d.iter_mut() {
        *di /= theta;
    }

    // backtrack along it into the box
    let mut alpha = 1.0;
    let mut temp1 = alpha;
    let mut ibd = 0;
    for (i, &k) in free.iter().enumerate() {
        let dk = d[i];
        if pr.nbd[k] != 0 {
            if dk < 0.0 && pr.nbd[k] <= 2 {
                let temp2 = pr.l[k] - z[k];
                if temp2 >= 0.0 {
                    temp1 = 0.0;
                } else if dk * alpha < temp2 {
                    temp1 = temp2 / dk;
                }
            } else if dk > 0.0 && pr.nbd[k] >= 2 {
                let temp2 = pr.u[k] - z[k];
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
        let dk = d[ibd];
        let k = free[ibd];
        if dk > 0.0 {
            z[k] = pr.u[k];
            d[ibd] = 0.0;
        } else if dk < 0.0 {
            z[k] = pr.l[k];
            d[ibd] = 0.0;
        }
    }
    for (i, &k) in free.iter().enumerate() {
        z[k] += alpha * d[i];
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Search {
    Evaluate,
    Converged,
    Warning,
}

/// The Moré-Thuente line search (MINPACK-2 `dcsrch`), as a state machine:
/// `start` then `step` with f and its directional derivative at each trial.
struct Dcsrch {
    brackt: bool,
    stage: u8,
    ginit: f64,
    gtest: f64,
    gx: f64,
    gy: f64,
    finit: f64,
    fx: f64,
    fy: f64,
    stx: f64,
    sty: f64,
    stmin: f64,
    stmax: f64,
    width: f64,
    width1: f64,
    stpmin: f64,
    stpmax: f64,
}

impl Dcsrch {
    fn start(f: f64, g: f64, stp: f64, stpmin: f64, stpmax: f64) -> Self {
        Self {
            brackt: false,
            stage: 1,
            finit: f,
            ginit: g,
            gtest: FTOL * g,
            width: stpmax - stpmin,
            width1: (stpmax - stpmin) / 0.5,
            stx: 0.0,
            fx: f,
            gx: g,
            sty: 0.0,
            fy: f,
            gy: g,
            stmin: 0.0,
            stmax: stp + 4.0 * stp,
            stpmin,
            stpmax,
        }
    }

    fn step(&mut self, f: f64, g: f64, stp: &mut f64) -> Search {
        let ftest = self.finit + *stp * self.gtest;
        if self.stage == 1 && f <= ftest && g >= 0.0 {
            self.stage = 2;
        }
        let mut status = Search::Evaluate;
        if self.brackt && (*stp <= self.stmin || *stp >= self.stmax) {
            status = Search::Warning;
        }
        if self.brackt && self.stmax - self.stmin <= XTOL * self.stmax {
            status = Search::Warning;
        }
        if *stp == self.stpmax && f <= ftest && g <= self.gtest {
            status = Search::Warning;
        }
        if *stp == self.stpmin && (f > ftest || g >= self.gtest) {
            status = Search::Warning;
        }
        if f <= ftest && g.abs() <= GTOL * (-self.ginit) {
            status = Search::Converged;
        }
        if status != Search::Evaluate {
            return status;
        }
        if self.stage == 1 && f <= self.fx && f > ftest {
            let gt = self.gtest;
            let fm = f - *stp * gt;
            let mut fxm = self.fx - self.stx * gt;
            let mut fym = self.fy - self.sty * gt;
            let gm = g - gt;
            let mut gxm = self.gx - gt;
            let mut gym = self.gy - gt;
            dcstep(
                &mut self.stx,
                &mut fxm,
                &mut gxm,
                &mut self.sty,
                &mut fym,
                &mut gym,
                stp,
                fm,
                gm,
                &mut self.brackt,
                self.stmin,
                self.stmax,
            );
            self.fx = fxm + self.stx * gt;
            self.fy = fym + self.sty * gt;
            self.gx = gxm + gt;
            self.gy = gym + gt;
        } else {
            dcstep(
                &mut self.stx,
                &mut self.fx,
                &mut self.gx,
                &mut self.sty,
                &mut self.fy,
                &mut self.gy,
                stp,
                f,
                g,
                &mut self.brackt,
                self.stmin,
                self.stmax,
            );
        }
        if self.brackt {
            if (self.sty - self.stx).abs() >= 0.66 * self.width1 {
                *stp = self.stx + 0.5 * (self.sty - self.stx);
            }
            self.width1 = self.width;
            self.width = (self.sty - self.stx).abs();
        }
        if self.brackt {
            self.stmin = self.stx.min(self.sty);
            self.stmax = self.stx.max(self.sty);
        } else {
            self.stmin = *stp + 1.1 * (*stp - self.stx);
            self.stmax = *stp + 4.0 * (*stp - self.stx);
        }
        *stp = stp.max(self.stpmin).min(self.stpmax);
        if self.brackt
            && (*stp <= self.stmin
                || *stp >= self.stmax
                || self.stmax - self.stmin <= XTOL * self.stmax)
        {
            *stp = self.stx;
        }
        Search::Evaluate
    }
}

/// One step of the safeguarded interval update (MINPACK-2 `dcstep`).
#[allow(clippy::too_many_arguments)]
fn dcstep(
    stx: &mut f64,
    fx: &mut f64,
    dx: &mut f64,
    sty: &mut f64,
    fy: &mut f64,
    dy: &mut f64,
    stp: &mut f64,
    fp: f64,
    dp: f64,
    brackt: &mut bool,
    stpmin: f64,
    stpmax: f64,
) {
    let sgnd = dp * (*dx / dx.abs());
    let stpf;
    if fp > *fx {
        let theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        let s = theta.abs().max(dx.abs()).max(dp.abs());
        let mut gamma = s * ((theta / s).powi(2) - (*dx / s) * (dp / s)).sqrt();
        if *stp < *stx {
            gamma = -gamma;
        }
        let p = (gamma - *dx) + theta;
        let q = ((gamma - *dx) + gamma) + dp;
        let r = p / q;
        let stpc = *stx + r * (*stp - *stx);
        let stpq = *stx + ((*dx / ((*fx - fp) / (*stp - *stx) + *dx)) / 2.0) * (*stp - *stx);
        stpf = if (stpc - *stx).abs() < (stpq - *stx).abs() {
            stpc
        } else {
            stpc + (stpq - stpc) / 2.0
        };
        *brackt = true;
    } else if sgnd < 0.0 {
        let theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        let s = theta.abs().max(dx.abs()).max(dp.abs());
        let mut gamma = s * ((theta / s).powi(2) - (*dx / s) * (dp / s)).sqrt();
        if *stp > *stx {
            gamma = -gamma;
        }
        let p = (gamma - dp) + theta;
        let q = ((gamma - dp) + gamma) + *dx;
        let r = p / q;
        let stpc = *stp + r * (*stx - *stp);
        let stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
        stpf = if (stpc - *stp).abs() > (stpq - *stp).abs() {
            stpc
        } else {
            stpq
        };
        *brackt = true;
    } else if dp.abs() < dx.abs() {
        let theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        let s = theta.abs().max(dx.abs()).max(dp.abs());
        let mut gamma = s * 0f64.max((theta / s).powi(2) - (*dx / s) * (dp / s)).sqrt();
        if *stp > *stx {
            gamma = -gamma;
        }
        let p = (gamma - dp) + theta;
        let q = (gamma + (*dx - dp)) + gamma;
        let r = p / q;
        let stpc = if r < 0.0 && gamma != 0.0 {
            *stp + r * (*stx - *stp)
        } else if *stp > *stx {
            stpmax
        } else {
            stpmin
        };
        let stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
        if *brackt {
            let mut f = if (stpc - *stp).abs() < (stpq - *stp).abs() {
                stpc
            } else {
                stpq
            };
            if *stp > *stx {
                f = f.min(*stp + 0.66 * (*sty - *stp));
            } else {
                f = f.max(*stp + 0.66 * (*sty - *stp));
            }
            stpf = f;
        } else {
            let f = if (stpc - *stp).abs() > (stpq - *stp).abs() {
                stpc
            } else {
                stpq
            };
            stpf = f.min(stpmax).max(stpmin);
        }
    } else if *brackt {
        let theta = 3.0 * (fp - *fy) / (*sty - *stp) + *dy + dp;
        let s = theta.abs().max(dy.abs()).max(dp.abs());
        let mut gamma = s * ((theta / s).powi(2) - (*dy / s) * (dp / s)).sqrt();
        if *stp > *sty {
            gamma = -gamma;
        }
        let p = (gamma - dp) + theta;
        let q = ((gamma - dp) + gamma) + *dy;
        let r = p / q;
        stpf = *stp + r * (*sty - *stp);
    } else if *stp > *stx {
        stpf = stpmax;
    } else {
        stpf = stpmin;
    }
    if fp > *fx {
        *sty = *stp;
        *fy = fp;
        *dy = dp;
    } else {
        if sgnd < 0.0 {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = fp;
        *dx = dp;
    }
    *stp = stpf;
}

/// A gradient function for [`optim_lbfgsb`].
pub type Gradient<'a> = &'a mut dyn FnMut(&[f64]) -> Result<Vec<f64>>;

/// Minimize `f` over the box [lower, upper] from `par`, as R's
/// `optim(par, f, gr, method = "L-BFGS-B", lower, upper, control)`.
///
/// `gr` is the gradient; without one, R's finite differences are used.
/// Infinite bounds leave a side unbounded.
pub fn optim_lbfgsb(
    f: &mut dyn FnMut(&[f64]) -> Result<f64>,
    mut gr: Option<Gradient>,
    par: &[f64],
    lower: &[f64],
    upper: &[f64],
    opts: &LbfgsbOptions,
) -> Result<LbfgsbResult> {
    let n = par.len();
    if lower.len() != n || upper.len() != n {
        return Err(PopmixError::Validation(
            "optim_lbfgsb: bounds must have one entry per parameter".to_string(),
        ));
    }
    let pr = Problem {
        l: lower,
        u: upper,
        nbd: (0..n).map(|i| bound_type(lower[i], upper[i])).collect(),
    };
    let fnscale = opts.fnscale;
    let evals = std::cell::Cell::new(0usize);
    let mut fg = |x: &[f64]| -> Result<(f64, Vec<f64>)> {
        evals.set(evals.get() + 1);
        let fx = f(x)? / fnscale;
        if !fx.is_finite() {
            return Err(PopmixError::Numerical(
                "L-BFGS-B needs finite values of the objective".to_string(),
            ));
        }
        let g = match gr.as_mut() {
            Some(gr) => gr(x)?.iter().map(|v| v / fnscale).collect(),
            None => {
                let mut g = vec![0.0; n];
                let mut xt = x.to_vec();
                for i in 0..n {
                    let mut eps = opts.ndeps;
                    let mut epsused = eps;
                    let mut tmp = x[i] + eps;
                    if tmp > upper[i] {
                        tmp = upper[i];
                        epsused = tmp - x[i];
                    }
                    xt[i] = tmp;
                    let val1 = f(&xt)? / fnscale;
                    tmp = x[i] - eps;
                    if tmp < lower[i] {
                        tmp = lower[i];
                        eps = x[i] - tmp;
                    }
                    xt[i] = tmp;
                    let val2 = f(&xt)? / fnscale;
                    g[i] = (val1 - val2) / (epsused + eps);
                    if !g[i].is_finite() {
                        return Err(PopmixError::Numerical(format!(
                            "non-finite finite-difference value [{}]",
                            i + 1
                        )));
                    }
                    xt[i] = x[i];
                }
                g
            }
        };
        Ok((fx, g))
    };

    // project the start into the box (`active`)
    let mut x = par.to_vec();
    let mut iwhere: Where = vec![0; n];
    let mut cnstnd = false;
    let mut boxed = true;
    for i in 0..n {
        if pr.nbd[i] > 0 {
            if pr.nbd[i] <= 2 && x[i] <= lower[i] {
                x[i] = lower[i];
            } else if pr.nbd[i] >= 2 && x[i] >= upper[i] {
                x[i] = upper[i];
            }
        }
        if pr.nbd[i] != 2 {
            boxed = false;
        }
        if pr.nbd[i] == 0 {
            iwhere[i] = -1;
        } else {
            cnstnd = true;
            iwhere[i] = if pr.nbd[i] == 2 && upper[i] - lower[i] <= 0.0 {
                3
            } else {
                0
            };
        }
    }

    let tol = opts.factr * EPSMCH;
    let mut mem = Memory::new(opts.lmm.max(1));
    let (mut fx, mut g) = fg(&x)?;
    let mut iter = 0usize;
    let finish = |x: Vec<f64>, fx: f64, convergence: i32, message: &str| LbfgsbResult {
        par: x,
        value: fx * fnscale,
        evaluations: evals.get(),
        convergence,
        message: message.to_string(),
    };
    let mut sbgnrm = pr.projected_gradient_norm(&x, &g);
    if sbgnrm <= opts.pgtol {
        return Ok(finish(
            x,
            fx,
            0,
            "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",
        ));
    }
    loop {
        // the search target z: the Cauchy point, then the subspace minimizer
        let (mut z, c) = if !cnstnd && mem.col() > 0 {
            (x.clone(), vec![0.0; 2 * mem.col()])
        } else {
            match cauchy(&pr, &x, &g, &mut iwhere, &mem, sbgnrm) {
                Ok(r) => r,
                Err(_) => {
                    mem.reset();
                    continue;
                }
            }
        };
        let free: Vec<usize> = (0..n).filter(|&i| iwhere[i] <= 0).collect();
        let active: Vec<usize> = (0..n).rev().filter(|&i| iwhere[i] > 0).collect();
        if !free.is_empty()
            && mem.col() > 0
            && subspace_min(&pr, &x, &g, &mut z, &c, &free, &active, &mem).is_err()
        {
            mem.reset();
            continue;
        }

        // line search along d = z - x (`lnsrlb`)
        let d: Vec<f64> = (0..n).map(|i| z[i] - x[i]).collect();
        let dtd = dot(&d, &d);
        let dnorm = dtd.sqrt();
        let mut stpmx = BIG;
        if cnstnd {
            if iter == 0 {
                stpmx = 1.0;
            } else {
                for i in 0..n {
                    let a1 = d[i];
                    if pr.nbd[i] != 0 {
                        if a1 < 0.0 && pr.nbd[i] <= 2 {
                            let a2 = lower[i] - x[i];
                            if a2 >= 0.0 {
                                stpmx = 0.0;
                            } else if a1 * stpmx < a2 {
                                stpmx = a2 / a1;
                            }
                        } else if a1 > 0.0 && pr.nbd[i] >= 2 {
                            let a2 = upper[i] - x[i];
                            if a2 <= 0.0 {
                                stpmx = 0.0;
                            } else if a1 * stpmx > a2 {
                                stpmx = a2 / a1;
                            }
                        }
                    }
                }
            }
        }
        let mut stp = if iter == 0 && !boxed {
            (1.0 / dnorm).min(stpmx)
        } else {
            1.0
        };
        let t = x.clone();
        let r = g.clone();
        let fold = fx;
        let mut ifun = 0usize;
        let mut gdold = 0.0;
        let mut search: Option<Dcsrch> = None;
        let mut failed = false;
        let gd = loop {
            let gd = dot(&g, &d);
            if ifun == 0 {
                gdold = gd;
                if gd >= 0.0 {
                    // not a descent direction
                    failed = true;
                    break gd;
                }
            }
            let status = match search.as_mut() {
                None => {
                    search = Some(Dcsrch::start(fx, gd, stp, 0.0, stpmx));
                    Search::Evaluate
                }
                Some(s) => s.step(fx, gd, &mut stp),
            };
            if status != Search::Evaluate {
                break gd;
            }
            ifun += 1;
            // the 21st trial is not evaluated
            if ifun > 20 {
                failed = true;
                break gd;
            }
            x = if stp == 1.0 {
                z.clone()
            } else {
                (0..n).map(|i| stp * d[i] + t[i]).collect()
            };
            (fx, g) = fg(&x)?;
        };
        if failed {
            x = t;
            g = r;
            fx = fold;
            if mem.col() == 0 {
                return Ok(finish(x, fx, 52, "ABNORMAL_TERMINATION_IN_LNSRCH"));
            }
            mem.reset();
            continue;
        }
        iter += 1;
        sbgnrm = pr.projected_gradient_norm(&x, &g);
        if iter > opts.maxit {
            return Ok(finish(x, fx, 1, "NEW_X"));
        }
        if sbgnrm <= opts.pgtol {
            return Ok(finish(
                x,
                fx,
                0,
                "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",
            ));
        }
        let ddum = fold.abs().max(fx.abs()).max(1.0);
        if fold - fx <= tol * ddum {
            return Ok(finish(
                x,
                fx,
                0,
                "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",
            ));
        }
        // the correction pair
        let yv: Vec<f64> = (0..n).map(|i| g[i] - r[i]).collect();
        let rr = dot(&yv, &yv);
        let (dr, ddum, s) = if stp == 1.0 {
            (gd - gdold, -gdold, d)
        } else {
            (
                (gd - gdold) * stp,
                -gdold * stp,
                d.iter().map(|v| stp * v).collect(),
            )
        };
        if dr <= EPSMCH * ddum {
            continue;
        }
        mem.update(s, yv, stp, dtd, dr, rr);
        if mem.form_t().is_err() {
            mem.reset();
        }
    }
}
