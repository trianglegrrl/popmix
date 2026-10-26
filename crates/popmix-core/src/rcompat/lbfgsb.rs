//! L-BFGS-B as R's `optim(method = "L-BFGS-B")` runs it.
//!
//! qpgraph optimises its admixture weights with this so it takes the same
//! steps as admixtools2, which calls R's `optim`. The objective is evaluated
//! the same number of times, at the same points, and the run stops with the
//! same convergence code and message as R (4.5.2; checked against R itself by
//! `tests/rcompat_golden.rs`).
//!
//! # Provenance
//!
//! The optimiser is derived from L-BFGS-B 3.0 (Fortran; Ciyou Zhu, Richard
//! Byrd, Peihuang Lu and Jorge Nocedal, with Jose Luis Morales and Jorge
//! Nocedal's 2011 correction), BSD-3, Copyright (c) 2011 Ciyou Zhu, Richard
//! Byrd, Jorge Nocedal and Jose Luis Morales; its licence is in
//! `LICENSE.lbfgsb.txt` beside this file. The algorithm is
//! described in Byrd, Lu, Nocedal and Zhu (1995), SIAM J. Sci. Comput.
//! 16:1190; Zhu, Byrd, Lu and Nocedal (1997), ACM TOMS 23:550; Morales and
//! Nocedal (2011), ACM TOMS 38:7. The line search is Moré and Thuente (1994),
//! ACM TOMS 20:286 (MINPACK-2's `dcsrch` and `dcstep`).
//!
//! R runs version 2.3 of the algorithm, which differs from 3.0 where it
//! matters for the steps taken: the subspace minimisation truncates the
//! Newton step at the first bound it meets (the largest step in [0, 1]
//! keeping the free variables feasible; Byrd et al. 1995, section 5) instead
//! of 3.0's projection, and machine epsilon is the double precision one.
//! This module does the same.
//!
//! The driver's behaviour (finite differences, counting, iteration limit,
//! messages and errors) was established from `?optim` and by running R as a
//! black box; see `tests/golden/rlang/generate_optim.R`.

use crate::error::{PopmixError, Result};

/// Settings of `optim(method = "L-BFGS-B")`; [`Default`] gives R's defaults.
#[derive(Debug, Clone)]
pub struct LbfgsbOptions {
    /// Number of correction pairs kept (R's `lmm`, 5).
    pub lmm: usize,
    /// Relative reduction tolerance in units of machine epsilon (1e7).
    pub factr: f64,
    /// Projected gradient tolerance (0).
    pub pgtol: f64,
    /// Iteration limit (100).
    pub maxit: usize,
    /// The objective is divided by this (1).
    pub fnscale: f64,
    /// Finite-difference step for every variable (1e-3).
    pub ndeps: f64,
}

impl Default for LbfgsbOptions {
    fn default() -> Self {
        LbfgsbOptions {
            lmm: 5,
            factr: 1e7,
            pgtol: 0.0,
            maxit: 100,
            fnscale: 1.0,
            ndeps: 1e-3,
        }
    }
}

/// What `optim` returns.
#[derive(Debug, Clone)]
pub struct LbfgsbResult {
    pub par: Vec<f64>,
    /// The objective at `par`, not divided by `fnscale`.
    pub value: f64,
    /// Objective (and gradient) evaluations, R's `counts`.
    pub evaluations: usize,
    /// 0 converged, 1 iteration limit reached, 51 warning, 52 error.
    pub convergence: i32,
    pub message: String,
}

/// An analytic gradient of the objective.
pub type Gradient<'a> = &'a mut dyn FnMut(&[f64]) -> Result<Vec<f64>>;

/// R's `optim(par, f, gr, method = "L-BFGS-B", lower, upper, control)`.
///
/// Without `gr` the gradient is R's central differences with step `ndeps`,
/// each side cut at its bound. Infinite bounds leave that side unbounded.
/// A start outside the box is projected into it.
pub fn optim_lbfgsb(
    f: &mut dyn FnMut(&[f64]) -> Result<f64>,
    gr: Option<Gradient>,
    par: &[f64],
    lower: &[f64],
    upper: &[f64],
    opts: &LbfgsbOptions,
) -> Result<LbfgsbResult> {
    let n = par.len();
    if lower.len() != n || upper.len() != n {
        return Err(PopmixError::Validation(format!(
            "optim_lbfgsb: {n} parameters but {} lower and {} upper bounds",
            lower.len(),
            upper.len()
        )));
    }
    let nbd: Vec<u8> = lower
        .iter()
        .zip(upper)
        .map(|(&l, &u)| match (l.is_finite(), u.is_finite()) {
            (false, false) => 0,
            (true, false) => 1,
            (true, true) => 2,
            (false, true) => 3,
        })
        .collect();
    let mut objective = Objective {
        f,
        gr,
        lower,
        upper,
        fnscale: opts.fnscale,
        ndeps: opts.ndeps,
        count: 0,
    };
    let mut x = par.to_vec();
    let problem = Problem {
        l: lower,
        u: upper,
        nbd: &nbd,
        m: opts.lmm,
        factr: opts.factr,
        pgtol: opts.pgtol,
        maxit: opts.maxit,
    };
    let (stop, value) = minimize(&problem, &mut x, &mut objective)?;
    let (convergence, message) = match stop {
        Stop::Pgtol => (
            0,
            "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL".to_string(),
        ),
        Stop::RelReduction => (
            0,
            "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH".to_string(),
        ),
        Stop::MaxIt => (1, "NEW_X".to_string()),
        Stop::Abnormal => (52, "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH".to_string()),
        Stop::Error(msg) => (52, msg.to_string()),
    };
    Ok(LbfgsbResult {
        par: x,
        value,
        evaluations: objective.count,
        convergence,
        message,
    })
}

/// The objective as the optimiser sees it: divided by `fnscale`, with its
/// gradient, analytic or by R's finite differences.
struct Objective<'a, 'g> {
    f: &'a mut dyn FnMut(&[f64]) -> Result<f64>,
    gr: Option<Gradient<'g>>,
    lower: &'a [f64],
    upper: &'a [f64],
    fnscale: f64,
    ndeps: f64,
    count: usize,
}

impl Objective<'_, '_> {
    /// The scaled objective, its unscaled value and the scaled gradient.
    fn eval(&mut self, x: &[f64]) -> Result<(f64, f64, Vec<f64>)> {
        let raw = (self.f)(x)?;
        self.count += 1;
        if !raw.is_finite() {
            return Err(PopmixError::Numerical(
                "L-BFGS-B needs finite values of 'fn'".to_string(),
            ));
        }
        let fx = raw / self.fnscale;
        let n = x.len();
        let g = match self.gr.as_mut() {
            Some(gr) => {
                let g = gr(x)?;
                if g.len() != n {
                    return Err(PopmixError::Validation(format!(
                        "gradient in optim evaluated to length {} not {n}",
                        g.len()
                    )));
                }
                if g.iter().any(|v| !v.is_finite()) {
                    return Err(PopmixError::Numerical(
                        "non-finite value supplied by optim".to_string(),
                    ));
                }
                g.iter().map(|v| v / self.fnscale).collect()
            }
            None => {
                let eps = self.ndeps;
                let mut xs = x.to_vec();
                let mut g = vec![0.0; n];
                for i in 0..n {
                    let xi = x[i];
                    let mut hi = xi + eps;
                    let mut e1 = eps;
                    if hi > self.upper[i] {
                        hi = self.upper[i];
                        e1 = hi - xi;
                    }
                    xs[i] = hi;
                    let f1 = (self.f)(&xs)? / self.fnscale;
                    let mut lo = xi - eps;
                    let mut e2 = eps;
                    if lo < self.lower[i] {
                        lo = self.lower[i];
                        e2 = xi - lo;
                    }
                    xs[i] = lo;
                    let f2 = (self.f)(&xs)? / self.fnscale;
                    xs[i] = xi;
                    g[i] = (f1 - f2) / (e1 + e2);
                    if !g[i].is_finite() {
                        return Err(PopmixError::Numerical(format!(
                            "non-finite finite-difference value [{}]",
                            i + 1
                        )));
                    }
                }
                g
            }
        };
        Ok((fx, raw, g))
    }
}

struct Problem<'a> {
    l: &'a [f64],
    u: &'a [f64],
    /// 0 unbounded, 1 lower bound only, 2 both, 3 upper bound only
    nbd: &'a [u8],
    m: usize,
    factr: f64,
    pgtol: f64,
    maxit: usize,
}

enum Stop {
    Pgtol,
    RelReduction,
    MaxIt,
    Abnormal,
    Error(&'static str),
}

const EPSMCH: f64 = f64::EPSILON;

// Line search parameters (L-BFGS-B's).
const FTOL: f64 = 1e-3;
const GTOL: f64 = 0.9;
const XTOL: f64 = 0.1;
const STPMAX_FREE: f64 = 1e10;

// Indices below are 1-based, as in the papers and the Fortran; `at` maps a
// 1-based (row, column) of a column-major matrix to its offset.
#[inline]
fn at(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

#[inline]
fn next(p: usize, m: usize) -> usize {
    p % m + 1
}

fn ddot(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for (a, b) in x.iter().zip(y) {
        s += a * b;
    }
    s
}

/// The limited-memory BFGS matrix in compact form: S and Y (n x m, a
/// circular buffer of columns from `head`), S'Y and S'S (m x m), the
/// Cholesky factor of theta S'S + L D^-1 L' (`wt`), and the factored middle
/// matrix of the subspace problem (`wn`, with its pieces kept in `wn1`).
struct Memory {
    n: usize,
    m: usize,
    ws: Vec<f64>,
    wy: Vec<f64>,
    sy: Vec<f64>,
    ss: Vec<f64>,
    wt: Vec<f64>,
    wn: Vec<f64>,
    wn1: Vec<f64>,
    theta: f64,
    col: usize,
    head: usize,
    itail: usize,
    iupdat: usize,
    updatd: bool,
}

impl Memory {
    fn new(n: usize, m: usize) -> Self {
        Memory {
            n,
            m,
            ws: vec![0.0; n * m],
            wy: vec![0.0; n * m],
            sy: vec![0.0; m * m],
            ss: vec![0.0; m * m],
            wt: vec![0.0; m * m],
            wn: vec![0.0; 4 * m * m],
            wn1: vec![0.0; 4 * m * m],
            theta: 1.0,
            col: 0,
            head: 1,
            itail: 0,
            iupdat: 0,
            updatd: false,
        }
    }

    /// Forget the corrections (after a numerical failure).
    fn reset(&mut self) {
        self.col = 0;
        self.head = 1;
        self.theta = 1.0;
        self.iupdat = 0;
        self.updatd = false;
    }

    #[inline]
    fn wy(&self, i: usize, j: usize) -> f64 {
        self.wy[at(i, j, self.n)]
    }

    #[inline]
    fn ws(&self, i: usize, j: usize) -> f64 {
        self.ws[at(i, j, self.n)]
    }
}

/// Solve with an upper triangular matrix `t` (column-major, leading
/// dimension `ldt`, its (1, 1) at offset `off`), of order `n`: T x = b
/// (`transpose` false) or T' x = b. Returns the index of a zero diagonal
/// element, or 0.
fn dtrsl(t: &[f64], off: usize, ldt: usize, n: usize, b: &mut [f64], transpose: bool) -> usize {
    let tt = |i: usize, j: usize| t[off + at(i, j, ldt)];
    for info in 1..=n {
        if tt(info, info) == 0.0 {
            return info;
        }
    }
    if n == 0 {
        return 0;
    }
    if transpose {
        b[0] /= tt(1, 1);
        for j in 2..=n {
            let mut s = 0.0;
            for k in 1..j {
                s += tt(k, j) * b[k - 1];
            }
            b[j - 1] -= s;
            b[j - 1] /= tt(j, j);
        }
    } else {
        b[n - 1] /= tt(n, n);
        for jj in 2..=n {
            let j = n - jj + 1;
            let temp = -b[j];
            for i in 1..=j {
                b[i - 1] += temp * tt(i, j + 1);
            }
            b[j - 1] /= tt(j, j);
        }
    }
    0
}

/// Cholesky factorisation A = R'R of the upper triangle of a symmetric
/// matrix (order `n`, at offset `off`, leading dimension `lda`), in place.
/// Returns 0, or the order of the leading minor that is not positive
/// definite.
fn dpofa(a: &mut [f64], off: usize, lda: usize, n: usize) -> usize {
    let ix = |i: usize, j: usize| off + at(i, j, lda);
    for j in 1..=n {
        let mut s = 0.0;
        for k in 1..j {
            let mut dot = 0.0;
            for i in 1..k {
                dot += a[ix(i, k)] * a[ix(i, j)];
            }
            let mut t = a[ix(k, j)] - dot;
            t /= a[ix(k, k)];
            a[ix(k, j)] = t;
            s += t * t;
        }
        s = a[ix(j, j)] - s;
        if s <= 0.0 {
            return j;
        }
        a[ix(j, j)] = s.sqrt();
    }
    0
}

/// p = M v for the 2col x 2col middle matrix M of the compact form, using
/// S'Y and the factor in `wt`. Returns a nonzero info when `wt` is singular.
fn bmv(mem: &Memory, v: &[f64], p: &mut [f64]) -> usize {
    let (m, col) = (mem.m, mem.col);
    if col == 0 {
        return 0;
    }
    let sy = |i: usize, j: usize| mem.sy[at(i, j, m)];
    p[col] = v[col];
    for i in 2..=col {
        let i2 = col + i;
        let mut sum = 0.0;
        for k in 1..i {
            sum += sy(i, k) * v[k - 1] / sy(k, k);
        }
        p[i2 - 1] = v[i2 - 1] + sum;
    }
    let info = dtrsl(&mem.wt, 0, m, col, &mut p[col..2 * col], true);
    if info != 0 {
        return info;
    }
    for i in 1..=col {
        p[i - 1] = v[i - 1] / sy(i, i).sqrt();
    }
    let info = dtrsl(&mem.wt, 0, m, col, &mut p[col..2 * col], false);
    if info != 0 {
        return info;
    }
    for i in 1..=col {
        p[i - 1] = -p[i - 1] / sy(i, i).sqrt();
    }
    for i in 1..=col {
        let mut sum = 0.0;
        for k in (i + 1)..=col {
            sum += sy(k, i) * p[col + k - 1] / sy(i, i);
        }
        p[i - 1] += sum;
    }
    0
}

/// Max norm of the projected gradient.
fn projgr(prob: &Problem, x: &[f64], g: &[f64]) -> f64 {
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
struct CauchyWork {
    p: Vec<f64>,
    c: Vec<f64>,
    wbp: Vec<f64>,
    v: Vec<f64>,
}

/// The generalised Cauchy point `xcp`: the first local minimiser of the
/// quadratic model along the projected steepest descent path x(t) =
/// P(x - t g), examined segment by segment between breakpoints taken from a
/// heap. Marks the variables it fixes at bounds in `iwhere`. Returns a
/// nonzero info on a singular system.
#[allow(clippy::too_many_arguments)]
fn cauchy(
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
    for v in w.p[..col2].iter_mut() {
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
        for v in w.p[col..col2].iter_mut() {
            *v *= theta;
        }
    }
    xcp.copy_from_slice(x);
    if nbreak == 0 && nfree == n + 1 {
        return 0;
    }
    for v in w.c[..col2].iter_mut() {
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
fn hpsolb(n: usize, t: &mut [f64], iorder: &mut [usize], iheap: usize) {
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

/// The free and active sets at the Cauchy point.
struct FreeSet {
    nfree: usize,
    /// free variables first (1..nfree), then the active ones
    index: Vec<usize>,
    /// variables entering the free set first (1..nenter), those leaving
    /// last (ileave..n)
    indx2: Vec<usize>,
    nenter: usize,
    ileave: usize,
}

/// Update the free set from `iwhere`; returns whether the subspace matrix
/// must be formed again.
fn freev(fs: &mut FreeSet, iwhere: &[i32], updatd: bool, cnstnd: bool, iter: usize) -> bool {
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

/// Form and factor the middle matrix K of the subspace problem,
///   K = [-D - Y'ZZ'Y/theta     L_a' - R_z'  ]
///       [L_a - R_z            theta S'AA'S ],
/// updating its pieces in `wn1` for the variables that entered or left the
/// free set. Returns a nonzero info if a factorisation fails.
fn formk(mem: &mut Memory, fs: &FreeSet) -> i32 {
    let (n, m, col, head, theta) = (mem.n, mem.m, mem.col, mem.head, mem.theta);
    let (nsub, nenter, ileave) = (fs.nfree, fs.nenter, fs.ileave);
    let (ind, indx2) = (&fs.index, &fs.indx2);
    let m2 = 2 * m;
    let w1 = |i: usize, j: usize| at(i, j, m2);
    let wyv = |mem: &Memory, i: usize, j: usize| mem.wy[at(i, j, n)];
    let wsv = |mem: &Memory, i: usize, j: usize| mem.ws[at(i, j, n)];

    let upcl = if mem.updatd {
        if mem.iupdat > m {
            // shift the old part of wn1 up-left
            for jy in 1..m {
                let js = m + jy;
                for k in 0..(m - jy) {
                    mem.wn1[w1(jy + k, jy)] = mem.wn1[w1(jy + 1 + k, jy + 1)];
                }
                for k in 0..(m - jy) {
                    mem.wn1[w1(js + k, js)] = mem.wn1[w1(js + 1 + k, js + 1)];
                }
                for k in 0..(m - 1) {
                    mem.wn1[w1(m + 1 + k, jy)] = mem.wn1[w1(m + 2 + k, jy + 1)];
                }
            }
        }
        // new row and column of the blocks from the free (pbegin..pend)
        // and active (dbegin..dend) variables
        let (pbegin, pend, dbegin, dend) = (1, nsub, nsub + 1, n);
        let iy = col;
        let is = m + col;
        let mut ipntr = head + col - 1;
        if ipntr > m {
            ipntr -= m;
        }
        let mut jpntr = head;
        for jy in 1..=col {
            let js = m + jy;
            let mut temp1 = 0.0;
            let mut temp2 = 0.0;
            let mut temp3 = 0.0;
            for k in pbegin..=pend {
                let k1 = ind[k - 1];
                temp1 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            for k in dbegin..=dend {
                let k1 = ind[k - 1];
                temp2 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            mem.wn1[w1(iy, jy)] = temp1;
            mem.wn1[w1(is, js)] = temp2;
            mem.wn1[w1(is, jy)] = temp3;
            jpntr = next(jpntr, m);
        }
        let jy = col;
        let mut jpntr = head + col - 1;
        if jpntr > m {
            jpntr -= m;
        }
        let mut ipntr = head;
        for i in 1..=col {
            let is = m + i;
            let mut temp3 = 0.0;
            for k in pbegin..=pend {
                let k1 = ind[k - 1];
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            ipntr = next(ipntr, m);
            mem.wn1[w1(is, jy)] = temp3;
        }
        col - 1
    } else {
        col
    };

    // the old part, for the variables that entered or left
    let mut ipntr = head;
    for iy in 1..=upcl {
        let is = m + iy;
        let mut jpntr = head;
        for jy in 1..=iy {
            let js = m + jy;
            let mut temp1 = 0.0;
            let mut temp2 = 0.0;
            let mut temp3 = 0.0;
            let mut temp4 = 0.0;
            for k in 1..=nenter {
                let k1 = indx2[k - 1];
                temp1 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
                temp2 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
            }
            for k in ileave..=n {
                let k1 = indx2[k - 1];
                temp3 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
                temp4 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
            }
            mem.wn1[w1(iy, jy)] = mem.wn1[w1(iy, jy)] + temp1 - temp3;
            mem.wn1[w1(is, js)] = mem.wn1[w1(is, js)] - temp2 + temp4;
            jpntr = next(jpntr, m);
        }
        ipntr = next(ipntr, m);
    }
    let mut ipntr = head;
    for is in (m + 1)..=(m + upcl) {
        let mut jpntr = head;
        for jy in 1..=upcl {
            let mut temp1 = 0.0;
            let mut temp3 = 0.0;
            for k in 1..=nenter {
                let k1 = indx2[k - 1];
                temp1 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            for k in ileave..=n {
                let k1 = indx2[k - 1];
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            if is <= jy + m {
                mem.wn1[w1(is, jy)] = mem.wn1[w1(is, jy)] + temp1 - temp3;
            } else {
                mem.wn1[w1(is, jy)] = mem.wn1[w1(is, jy)] - temp1 + temp3;
            }
            jpntr = next(jpntr, m);
        }
        ipntr = next(ipntr, m);
    }

    // the upper triangle of K, in wn
    for iy in 1..=col {
        let is = col + iy;
        let is1 = m + iy;
        for jy in 1..=iy {
            let js = col + jy;
            let js1 = m + jy;
            mem.wn[w1(jy, iy)] = mem.wn1[w1(iy, jy)] / theta;
            mem.wn[w1(js, is)] = mem.wn1[w1(is1, js1)] * theta;
        }
        for jy in 1..iy {
            mem.wn[w1(jy, is)] = -mem.wn1[w1(is1, jy)];
        }
        for jy in iy..=col {
            mem.wn[w1(jy, is)] = mem.wn1[w1(is1, jy)];
        }
        mem.wn[w1(iy, iy)] += mem.sy[at(iy, iy, m)];
    }

    // LEL' factorisation: factor the (1,1) block, solve for the (1,2)
    // block, update and factor the (2,2) block
    if dpofa(&mut mem.wn, 0, m2, col) != 0 {
        return -1;
    }
    let col2 = 2 * col;
    let mut b = vec![0.0; col];
    for js in (col + 1)..=col2 {
        for i in 1..=col {
            b[i - 1] = mem.wn[w1(i, js)];
        }
        dtrsl(&mem.wn, 0, m2, col, &mut b, true);
        for i in 1..=col {
            mem.wn[w1(i, js)] = b[i - 1];
        }
    }
    for is in (col + 1)..=col2 {
        for js in is..=col2 {
            let mut dot = 0.0;
            for k in 1..=col {
                dot += mem.wn[w1(k, is)] * mem.wn[w1(k, js)];
            }
            mem.wn[w1(is, js)] += dot;
        }
    }
    if dpofa(&mut mem.wn, w1(col + 1, col + 1), m2, col) != 0 {
        return -2;
    }
    0
}

/// Form and factor T = theta S'S + L D^-1 L' (upper triangle) into `wt`.
fn formt(mem: &mut Memory) -> i32 {
    let (m, col, theta) = (mem.m, mem.col, mem.theta);
    let sy = |mem: &Memory, i: usize, j: usize| mem.sy[at(i, j, m)];
    for j in 1..=col {
        mem.wt[at(1, j, m)] = theta * mem.ss[at(1, j, m)];
    }
    for i in 2..=col {
        for j in i..=col {
            let k1 = i.min(j) - 1;
            let mut ddum = 0.0;
            for k in 1..=k1 {
                ddum += sy(mem, i, k) * sy(mem, j, k) / sy(mem, k, k);
            }
            mem.wt[at(i, j, m)] = ddum + theta * mem.ss[at(i, j, m)];
        }
    }
    if dpofa(&mut mem.wt, 0, m, col) != 0 {
        return -3;
    }
    0
}

/// Add the correction pair (s, y) = (d, r) to the memory.
fn matupd(mem: &mut Memory, d: &[f64], r: &[f64], rr: f64, dr: f64, stp: f64, dtd: f64) {
    let (n, m) = (mem.n, mem.m);
    if mem.iupdat <= m {
        mem.col = mem.iupdat;
        mem.itail = (mem.head + mem.iupdat - 2) % m + 1;
    } else {
        mem.itail = next(mem.itail, m);
        mem.head = next(mem.head, m);
    }
    let it = mem.itail;
    mem.ws[at(1, it, n)..at(1, it, n) + n].copy_from_slice(d);
    mem.wy[at(1, it, n)..at(1, it, n) + n].copy_from_slice(r);
    mem.theta = rr / dr;
    let col = mem.col;
    if mem.iupdat > m {
        // drop the oldest pair from S'S and S'Y
        for j in 1..col {
            for k in 0..j {
                mem.ss[at(1 + k, j, m)] = mem.ss[at(2 + k, j + 1, m)];
            }
            for k in 0..(col - j) {
                mem.sy[at(j + k, j, m)] = mem.sy[at(j + 1 + k, j + 1, m)];
            }
        }
    }
    let mut pointr = mem.head;
    for j in 1..col {
        let wy = &mem.wy[at(1, pointr, n)..at(1, pointr, n) + n];
        let ws = &mem.ws[at(1, pointr, n)..at(1, pointr, n) + n];
        let syv = ddot(d, wy);
        let ssv = ddot(ws, d);
        mem.sy[at(col, j, m)] = syv;
        mem.ss[at(j, col, m)] = ssv;
        pointr = next(pointr, m);
    }
    mem.ss[at(col, col, m)] = if stp == 1.0 { dtd } else { stp * stp * dtd };
    mem.sy[at(col, col, m)] = dr;
}

/// The reduced gradient of the model at the Cauchy point z over the free
/// variables: r = -Z'(B(z - x) + g). Returns a nonzero info on failure.
#[allow(clippy::too_many_arguments)]
fn cmprlb(
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
fn subsm(
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
    for v in wv[..col].iter_mut() {
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
    for v in d[..nsub].iter_mut() {
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

#[derive(Debug, Clone, Copy, PartialEq)]
enum SearchTask {
    Start,
    Fg,
    Convergence,
    Warning,
    Error,
}

/// The state dcsrch keeps between calls.
#[derive(Debug, Clone, Default)]
struct SearchState {
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
}

/// Moré-Thuente line search (MINPACK-2 `dcsrch`): called with the function
/// value `f` and derivative `g` at `stp`, it updates `stp` and sets `task`
/// to `Fg` to ask for the next evaluation, or to say it has finished.
#[allow(clippy::too_many_arguments)]
fn dcsrch(
    f: f64,
    g: f64,
    stp: &mut f64,
    ftol: f64,
    gtol: f64,
    xtol: f64,
    stpmin: f64,
    stpmax: f64,
    task: &mut SearchTask,
    s: &mut SearchState,
) {
    const P5: f64 = 0.5;
    const P66: f64 = 0.66;
    const XTRAPL: f64 = 1.1;
    const XTRAPU: f64 = 4.0;
    if *task == SearchTask::Start {
        if *stp < stpmin
            || *stp > stpmax
            || g >= 0.0
            || ftol < 0.0
            || gtol < 0.0
            || xtol < 0.0
            || stpmin < 0.0
            || stpmax < stpmin
        {
            *task = SearchTask::Error;
            return;
        }
        s.brackt = false;
        s.stage = 1;
        s.finit = f;
        s.ginit = g;
        s.gtest = ftol * s.ginit;
        s.width = stpmax - stpmin;
        s.width1 = s.width / P5;
        s.stx = 0.0;
        s.fx = s.finit;
        s.gx = s.ginit;
        s.sty = 0.0;
        s.fy = s.finit;
        s.gy = s.ginit;
        s.stmin = 0.0;
        s.stmax = *stp + XTRAPU * *stp;
        *task = SearchTask::Fg;
        return;
    }
    let ftest = s.finit + *stp * s.gtest;
    if s.stage == 1 && f <= ftest && g >= 0.0 {
        s.stage = 2;
    }
    let mut done = None;
    if s.brackt && (*stp <= s.stmin || *stp >= s.stmax) {
        done = Some(SearchTask::Warning); // rounding errors prevent progress
    }
    if s.brackt && s.stmax - s.stmin <= xtol * s.stmax {
        done = Some(SearchTask::Warning); // xtol test satisfied
    }
    if *stp == stpmax && f <= ftest && g <= s.gtest {
        done = Some(SearchTask::Warning); // stp = stpmax
    }
    if *stp == stpmin && (f > ftest || g >= s.gtest) {
        done = Some(SearchTask::Warning); // stp = stpmin
    }
    if f <= ftest && g.abs() <= gtol * (-s.ginit) {
        done = Some(SearchTask::Convergence);
    }
    if let Some(t) = done {
        *task = t;
        return;
    }
    if s.stage == 1 && f <= s.fx && f > ftest {
        // the modified function
        let fm = f - *stp * s.gtest;
        let mut fxm = s.fx - s.stx * s.gtest;
        let mut fym = s.fy - s.sty * s.gtest;
        let gm = g - s.gtest;
        let mut gxm = s.gx - s.gtest;
        let mut gym = s.gy - s.gtest;
        dcstep(
            &mut s.stx,
            &mut fxm,
            &mut gxm,
            &mut s.sty,
            &mut fym,
            &mut gym,
            stp,
            fm,
            gm,
            &mut s.brackt,
            s.stmin,
            s.stmax,
        );
        s.fx = fxm + s.stx * s.gtest;
        s.fy = fym + s.sty * s.gtest;
        s.gx = gxm + s.gtest;
        s.gy = gym + s.gtest;
    } else {
        dcstep(
            &mut s.stx,
            &mut s.fx,
            &mut s.gx,
            &mut s.sty,
            &mut s.fy,
            &mut s.gy,
            stp,
            f,
            g,
            &mut s.brackt,
            s.stmin,
            s.stmax,
        );
    }
    if s.brackt {
        if (s.sty - s.stx).abs() >= P66 * s.width1 {
            *stp = s.stx + P5 * (s.sty - s.stx);
        }
        s.width1 = s.width;
        s.width = (s.sty - s.stx).abs();
    }
    if s.brackt {
        s.stmin = s.stx.min(s.sty);
        s.stmax = s.stx.max(s.sty);
    } else {
        s.stmin = *stp + XTRAPL * (*stp - s.stx);
        s.stmax = *stp + XTRAPU * (*stp - s.stx);
    }
    *stp = stp.max(stpmin);
    *stp = stp.min(stpmax);
    if s.brackt && (*stp <= s.stmin || *stp >= s.stmax || s.stmax - s.stmin <= xtol * s.stmax) {
        *stp = s.stx;
    }
    *task = SearchTask::Fg;
}

/// One safeguarded step of the Moré-Thuente search: a new trial step from
/// the best step `stx`, the other end of the interval `sty` and the current
/// trial `stp` (with values and derivatives), by cubic and quadratic
/// interpolation; updates the interval.
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
        // higher function value: the minimum is bracketed
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
        if (stpc - *stx).abs() < (stpq - *stx).abs() {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.0;
        }
        *brackt = true;
    } else if sgnd < 0.0 {
        // derivatives of opposite sign: bracketed
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
        if (stpc - *stp).abs() > (stpq - *stp).abs() {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = true;
    } else if dp.abs() < dx.abs() {
        // same sign, the derivative's magnitude decreases
        let theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        let s = theta.abs().max(dx.abs()).max(dp.abs());
        let mut gamma = s
            * (0.0f64)
                .max((theta / s).powi(2) - (*dx / s) * (dp / s))
                .sqrt();
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
            let mut sf = if (stpc - *stp).abs() < (stpq - *stp).abs() {
                stpc
            } else {
                stpq
            };
            if *stp > *stx {
                sf = sf.min(*stp + 0.66 * (*sty - *stp));
            } else {
                sf = sf.max(*stp + 0.66 * (*sty - *stp));
            }
            stpf = sf;
        } else {
            let mut sf = if (stpc - *stp).abs() > (stpq - *stp).abs() {
                stpc
            } else {
                stpq
            };
            sf = stpmax.min(sf);
            sf = stpmin.max(sf);
            stpf = sf;
        }
    } else {
        // same sign, the derivative's magnitude does not decrease
        if *brackt {
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

/// Input checks of the Fortran (`errclb`): the message of the last failing
/// check, if any.
fn check_input(prob: &Problem, n: usize) -> Option<&'static str> {
    let mut task = None;
    if n == 0 {
        task = Some("ERROR: N .LE. 0");
    }
    if prob.m == 0 {
        task = Some("ERROR: M .LE. 0");
    }
    if prob.factr < 0.0 {
        task = Some("ERROR: FACTR .LT. 0");
    }
    for i in 0..n {
        if prob.nbd[i] == 2 && prob.l[i] > prob.u[i] {
            task = Some("ERROR: NO FEASIBLE SOLUTION");
        }
    }
    task
}

/// Project the start into the box; returns (constrained, boxed) and sets
/// `iwhere` (-1 unbounded, 3 fixed by equal bounds, 0 otherwise).
fn active(prob: &Problem, x: &mut [f64], iwhere: &mut [i32]) -> (bool, bool) {
    let (l, u, nbd) = (prob.l, prob.u, prob.nbd);
    let mut cnstnd = false;
    let mut boxed = true;
    for i in 0..x.len() {
        if nbd[i] > 0 {
            if nbd[i] <= 2 && x[i] <= l[i] {
                if x[i] < l[i] {
                    x[i] = l[i];
                }
            } else if nbd[i] >= 2 && x[i] > u[i] {
                x[i] = u[i];
            }
        }
    }
    for i in 0..x.len() {
        if nbd[i] != 2 {
            boxed = false;
        }
        if nbd[i] == 0 {
            iwhere[i] = -1;
        } else {
            cnstnd = true;
            iwhere[i] = if nbd[i] == 2 && u[i] - l[i] <= 0.0 {
                3
            } else {
                0
            };
        }
    }
    (cnstnd, boxed)
}

/// The L-BFGS-B iteration with R's driver around it. Returns why it stopped
/// and the unscaled objective at `x`.
fn minimize(prob: &Problem, x: &mut [f64], obj: &mut Objective) -> Result<(Stop, f64)> {
    let n = x.len();
    if let Some(msg) = check_input(prob, n) {
        return Ok((Stop::Error(msg), 0.0));
    }
    let m = prob.m;
    let tol = prob.factr * EPSMCH;
    let mut mem = Memory::new(n, m);
    let mut iwhere = vec![0i32; n];
    let (cnstnd, boxed) = active(prob, x, &mut iwhere);
    let mut fs = FreeSet {
        nfree: n,
        index: vec![0; n],
        indx2: vec![0; n],
        nenter: 0,
        ileave: n + 1,
    };
    let mut work = CauchyWork {
        p: vec![0.0; 2 * m],
        c: vec![0.0; 2 * m],
        wbp: vec![0.0; 2 * m],
        v: vec![0.0; 2 * m],
    };
    let mut z = vec![0.0; n]; // Cauchy point, then the subspace minimiser
    let mut r = vec![0.0; n]; // reduced gradient; the previous gradient
    let mut d = vec![0.0; n]; // search direction
    let mut t = vec![0.0; n]; // breakpoints; the previous x
    let mut search = SearchState::default();

    let (mut f, mut f_raw, mut g) = obj.eval(x)?;
    let mut sbgnrm = projgr(prob, x, &g);
    if sbgnrm <= prob.pgtol {
        return Ok((Stop::Pgtol, f_raw));
    }
    let mut iter = 0usize;
    loop {
        // the generalised Cauchy point and the free set at it
        let wrk = if !cnstnd && mem.col > 0 {
            z.copy_from_slice(x);
            mem.updatd
        } else {
            let info = cauchy(
                prob,
                x,
                &g,
                &mut fs.indx2,
                &mut iwhere,
                &mut t,
                &mut d,
                &mut z,
                &mem,
                &mut work,
                sbgnrm,
            );
            if info != 0 {
                mem.reset();
                continue;
            }
            freev(&mut fs, &iwhere, mem.updatd, cnstnd, iter)
        };

        // subspace minimisation from the Cauchy point
        if fs.nfree != 0 && mem.col != 0 {
            if wrk && formk(&mut mem, &fs) != 0 {
                mem.reset();
                continue;
            }
            let mut info = cmprlb(x, &g, &z, &mut r, &mem, &mut work, &fs, cnstnd);
            if info == 0 {
                info = subsm(prob, &fs, &mut z, &mut r, &mut work.p, &mem);
            }
            if info != 0 {
                mem.reset();
                continue;
            }
        }

        // line search along d = z - x
        for i in 0..n {
            d[i] = z[i] - x[i];
        }
        let dtd = ddot(&d, &d);
        let dnorm = dtd.sqrt();
        let mut stpmx = STPMAX_FREE;
        if cnstnd {
            if iter == 0 {
                stpmx = 1.0;
            } else {
                for i in 0..n {
                    let a1 = d[i];
                    let b = prob.nbd[i];
                    if b != 0 {
                        if a1 < 0.0 && b <= 2 {
                            let a2 = prob.l[i] - x[i];
                            if a2 >= 0.0 {
                                stpmx = 0.0;
                            } else if a1 * stpmx < a2 {
                                stpmx = a2 / a1;
                            }
                        } else if a1 > 0.0 && b >= 2 {
                            let a2 = prob.u[i] - x[i];
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
        t.copy_from_slice(x);
        r.copy_from_slice(&g);
        let fold = f;
        let fold_raw = f_raw;
        let mut ifun = 0usize;
        let mut task = SearchTask::Start;
        let mut gd;
        let mut gdold = 0.0;
        let failed = loop {
            gd = ddot(&g, &d);
            if ifun == 0 {
                gdold = gd;
                if gd >= 0.0 {
                    // not a descent direction
                    break true;
                }
            }
            dcsrch(
                f,
                gd,
                &mut stp,
                FTOL,
                GTOL,
                XTOL,
                0.0,
                stpmx,
                &mut task,
                &mut search,
            );
            if matches!(task, SearchTask::Convergence | SearchTask::Warning) {
                break false;
            }
            ifun += 1;
            let iback = ifun - 1;
            if stp == 1.0 {
                x.copy_from_slice(&z);
            } else {
                for i in 0..n {
                    x[i] = stp * d[i] + t[i];
                }
            }
            if iback >= 20 {
                break true;
            }
            let (fv, fr, gv) = obj.eval(x)?;
            f = fv;
            f_raw = fr;
            g = gv;
        };
        if failed {
            x.copy_from_slice(&t);
            g.copy_from_slice(&r);
            f = fold;
            f_raw = fold_raw;
            if mem.col == 0 {
                return Ok((Stop::Abnormal, f_raw));
            }
            mem.reset();
            continue;
        }
        iter += 1;
        sbgnrm = projgr(prob, x, &g);

        // R's iteration limit
        if iter > prob.maxit {
            return Ok((Stop::MaxIt, f_raw));
        }

        // convergence tests
        if sbgnrm <= prob.pgtol {
            return Ok((Stop::Pgtol, f_raw));
        }
        let ddum = fold.abs().max(f.abs()).max(1.0);
        if fold - f <= tol * ddum {
            return Ok((Stop::RelReduction, f_raw));
        }

        // BFGS update, skipped without enough curvature
        for i in 0..n {
            r[i] = g[i] - r[i];
        }
        let rr = ddot(&r, &r);
        let dr;
        let ddum;
        if stp == 1.0 {
            dr = gd - gdold;
            ddum = -gdold;
        } else {
            dr = (gd - gdold) * stp;
            for v in d.iter_mut() {
                *v *= stp;
            }
            ddum = -gdold * stp;
        }
        if dr <= EPSMCH * ddum {
            mem.updatd = false;
            continue;
        }
        mem.updatd = true;
        mem.iupdat += 1;
        matupd(&mut mem, &d, &r, rr, dr, stp, dtd);
        if formt(&mut mem) != 0 {
            mem.reset();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn problem<'a>(l: &'a [f64], u: &'a [f64], nbd: &'a [u8], m: usize) -> Problem<'a> {
        Problem {
            l,
            u,
            nbd,
            m,
            factr: 1e7,
            pgtol: 0.0,
            maxit: 100,
        }
    }

    #[test]
    fn test_dcsrch_quadratic_takes_exact_minimiser() {
        // phi(a) = (a - 2)^2 from a = 0: phi'(0) = -4; the first trial at 1
        // has phi' = -2, so the cubic step interpolates to the minimiser 2.
        let phi = |a: f64| ((a - 2.0).powi(2), 2.0 * (a - 2.0));
        let mut s = SearchState::default();
        let mut task = SearchTask::Start;
        let mut stp = 1.0;
        let (f0, g0) = phi(0.0);
        dcsrch(
            f0, g0, &mut stp, FTOL, GTOL, XTOL, 0.0, 1e10, &mut task, &mut s,
        );
        assert_eq!(task, SearchTask::Fg);
        assert_eq!(stp, 1.0);
        let mut trials = vec![];
        for _ in 0..10 {
            let (f, g) = phi(stp);
            dcsrch(
                f, g, &mut stp, FTOL, GTOL, XTOL, 0.0, 1e10, &mut task, &mut s,
            );
            trials.push(stp);
            if task != SearchTask::Fg {
                break;
            }
        }
        // first trial at 1 satisfies sufficient decrease (1 <= 4 - 0.004)
        // but not the curvature condition (|-2| > 0.9 * 4 is false: 2 <=
        // 3.6), so it converges at once
        assert_eq!(task, SearchTask::Convergence);
        assert_eq!(trials, vec![1.0]);
    }

    #[test]
    fn test_dcsrch_rejects_ascent() {
        let mut s = SearchState::default();
        let mut task = SearchTask::Start;
        let mut stp = 1.0;
        dcsrch(
            1.0, 0.5, &mut stp, FTOL, GTOL, XTOL, 0.0, 10.0, &mut task, &mut s,
        );
        assert_eq!(task, SearchTask::Error);
    }

    #[test]
    fn test_dcsrch_backtracks_on_increase() {
        // phi(a) = a^4 - a from 0 with a trial step 10: f rises, so the
        // minimum is bracketed and the next trial lies inside (0, 10)
        let phi = |a: f64| (a.powi(4) - a, 4.0 * a.powi(3) - 1.0);
        let mut s = SearchState::default();
        let mut task = SearchTask::Start;
        let mut stp = 10.0;
        let (f0, g0) = phi(0.0);
        dcsrch(
            f0, g0, &mut stp, FTOL, GTOL, XTOL, 0.0, 1e10, &mut task, &mut s,
        );
        let (f, g) = phi(stp);
        dcsrch(
            f, g, &mut stp, FTOL, GTOL, XTOL, 0.0, 1e10, &mut task, &mut s,
        );
        assert_eq!(task, SearchTask::Fg);
        assert!(s.brackt);
        assert!(stp > 0.0 && stp < 10.0, "{stp}");
        let mut k = 0;
        while task == SearchTask::Fg && k < 30 {
            let (f, g) = phi(stp);
            dcsrch(
                f, g, &mut stp, FTOL, GTOL, XTOL, 0.0, 1e10, &mut task, &mut s,
            );
            k += 1;
        }
        assert_eq!(task, SearchTask::Convergence);
        // strong Wolfe conditions at the accepted step
        let (f, g) = phi(stp);
        assert!(f <= f0 + FTOL * stp * g0);
        assert!(g.abs() <= GTOL * g0.abs());
    }

    #[test]
    fn test_cauchy_without_memory_is_projected_steepest_descent() {
        // B = I (no corrections): along x - t g the model is minimised at
        // t = 1 unless a bound is met first. x = (0.5, 0.5), g = (1, -0.2),
        // box [0, 1]^2: variable 1 hits 0 at t = 0.5; then only variable 2
        // moves, to its unconstrained minimiser x2 + 0.2 = 0.7.
        let (l, u, nbd) = ([0.0, 0.0], [1.0, 1.0], [2u8, 2]);
        let prob = problem(&l, &u, &nbd, 5);
        let mem = Memory::new(2, 5);
        let x = [0.5, 0.5];
        let g = [1.0, -0.2];
        let mut iorder = [0usize; 2];
        let mut iwhere = [0i32; 2];
        let (mut t, mut d, mut xcp) = ([0.0; 2], [0.0; 2], [0.0; 2]);
        let mut w = CauchyWork {
            p: vec![0.0; 10],
            c: vec![0.0; 10],
            wbp: vec![0.0; 10],
            v: vec![0.0; 10],
        };
        let sb = projgr(&prob, &x, &g);
        let info = cauchy(
            &prob,
            &x,
            &g,
            &mut iorder,
            &mut iwhere,
            &mut t,
            &mut d,
            &mut xcp,
            &mem,
            &mut w,
            sb,
        );
        assert_eq!(info, 0);
        assert_eq!(xcp[0], 0.0);
        assert!((xcp[1] - 0.7).abs() < 1e-15, "{xcp:?}");
        assert_eq!(iwhere, [1, 0]);
    }

    #[test]
    fn test_cauchy_stops_at_bound_held_variable() {
        // at the lower bound with the gradient pushing outward the variable
        // is held; the other one moves freely
        let (l, u, nbd) = ([0.0, f64::NEG_INFINITY], [1.0, f64::INFINITY], [2u8, 0]);
        let prob = problem(&l, &u, &nbd, 5);
        let mem = Memory::new(2, 5);
        let x = [0.0, 3.0];
        let g = [2.0, 1.5];
        let mut iorder = [0usize; 2];
        let mut iwhere = [0i32, -1];
        let (mut t, mut d, mut xcp) = ([0.0; 2], [0.0; 2], [0.0; 2]);
        let mut w = CauchyWork {
            p: vec![0.0; 10],
            c: vec![0.0; 10],
            wbp: vec![0.0; 10],
            v: vec![0.0; 10],
        };
        let sb = projgr(&prob, &x, &g);
        cauchy(
            &prob,
            &x,
            &g,
            &mut iorder,
            &mut iwhere,
            &mut t,
            &mut d,
            &mut xcp,
            &mem,
            &mut w,
            sb,
        );
        assert_eq!(xcp, [0.0, 1.5]);
        assert_eq!(iwhere[0], 1);
    }

    #[test]
    fn test_subsm_truncates_at_first_bound() {
        // theta = 1 and no corrections: the Newton direction is r itself.
        // From xcp = (0.5, 0.5) along r = (1, 0.25) in [0, 1]^2 variable 1
        // reaches 1 at alpha = 0.5, so the step is cut there (2.3's rule;
        // 3.0 would project to (1, 0.75)).
        let (l, u, nbd) = ([0.0, 0.0], [1.0, 1.0], [2u8, 2]);
        let prob = problem(&l, &u, &nbd, 5);
        let mem = Memory::new(2, 5);
        let fs = FreeSet {
            nfree: 2,
            index: vec![1, 2],
            indx2: vec![0, 0],
            nenter: 0,
            ileave: 3,
        };
        let mut x = [0.5, 0.5];
        let mut d = [1.0, 0.25];
        let mut wv = vec![0.0; 10];
        assert_eq!(subsm(&prob, &fs, &mut x, &mut d, &mut wv, &mem), 0);
        assert_eq!(x, [1.0, 0.625]);
    }

    #[test]
    fn test_subsm_inside_box_takes_full_step() {
        let (l, u, nbd) = ([0.0, 0.0], [1.0, 1.0], [2u8, 2]);
        let prob = problem(&l, &u, &nbd, 5);
        let mem = Memory::new(2, 5);
        let fs = FreeSet {
            nfree: 2,
            index: vec![1, 2],
            indx2: vec![0, 0],
            nenter: 0,
            ileave: 3,
        };
        let mut x = [0.5, 0.5];
        let mut d = [0.25, -0.125];
        let mut wv = vec![0.0; 10];
        subsm(&prob, &fs, &mut x, &mut d, &mut wv, &mem);
        assert_eq!(x, [0.75, 0.375]);
    }

    #[test]
    fn test_hpsolb_yields_breakpoints_in_order() {
        let mut t = vec![0.7, 0.2, 0.9, 0.1, 0.5];
        let mut iorder = vec![1, 2, 3, 4, 5];
        let mut got = vec![];
        let n = t.len();
        for (k, nleft) in (1..=n).rev().enumerate() {
            hpsolb(nleft, &mut t, &mut iorder, k);
            got.push((t[nleft - 1], iorder[nleft - 1]));
        }
        assert_eq!(got, vec![(0.1, 4), (0.2, 2), (0.5, 5), (0.7, 1), (0.9, 3)]);
    }

    #[test]
    fn test_dpofa_and_dtrsl_solve_spd_system() {
        // A = [[4, 2], [2, 3]] = R'R; solve A x = b by R' y = b, R x = y
        let mut a = vec![4.0, 2.0, 2.0, 3.0];
        assert_eq!(dpofa(&mut a, 0, 2, 2), 0);
        let mut b = vec![2.0, 1.0]; // x = (0.5, 0)
        dtrsl(&a, 0, 2, 2, &mut b, true);
        dtrsl(&a, 0, 2, 2, &mut b, false);
        assert!((b[0] - 0.5).abs() < 1e-15 && b[1].abs() < 1e-15, "{b:?}");
        let mut c = vec![1.0, 2.0, 2.0, 1.0];
        assert_eq!(dpofa(&mut c, 0, 2, 2), 2);
    }

    fn quad(x: &[f64]) -> f64 {
        (x[0] - 0.3).powi(2) + (x[1] - 0.7).powi(2)
    }

    #[test]
    fn test_optim_projects_start_and_converges() {
        let mut f = |x: &[f64]| Ok(quad(x));
        let r = optim_lbfgsb(
            &mut f,
            None,
            &[2.0, -1.0],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap();
        assert_eq!(r.convergence, 0);
        assert!((r.par[0] - 0.3).abs() < 1e-6 && (r.par[1] - 0.7).abs() < 1e-6);
    }

    #[test]
    fn test_optim_maxit_zero_stops_after_first_iteration() {
        // R: maxit = 0 still takes one iteration, then reports "NEW_X"
        let mut f = |x: &[f64]| Ok(quad(x));
        let opts = LbfgsbOptions {
            maxit: 0,
            ..LbfgsbOptions::default()
        };
        let r = optim_lbfgsb(&mut f, None, &[0.1, 0.1], &[0.0, 0.0], &[1.0, 1.0], &opts).unwrap();
        assert_eq!((r.convergence, r.message.as_str()), (1, "NEW_X"));
        assert_eq!(r.evaluations, 2);
    }

    #[test]
    fn test_optim_errors_like_r() {
        let mut f = |x: &[f64]| Ok(quad(x));
        let r = optim_lbfgsb(
            &mut f,
            None,
            &[0.1, 0.1],
            &[1.0, 0.0],
            &[0.0, 1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap();
        assert_eq!(r.convergence, 52);
        assert_eq!(r.message, "ERROR: NO FEASIBLE SOLUTION");
        assert_eq!(r.evaluations, 0);

        let mut inf = |_: &[f64]| Ok(f64::INFINITY);
        let e = optim_lbfgsb(
            &mut inf,
            None,
            &[0.1],
            &[0.0],
            &[1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap_err();
        assert!(e
            .to_string()
            .contains("L-BFGS-B needs finite values of 'fn'"));

        // equal bounds: the finite difference is 0/0, as in R
        let e = optim_lbfgsb(
            &mut f,
            None,
            &[0.5, 0.1],
            &[0.5, 0.0],
            &[0.5, 1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap_err();
        assert!(e
            .to_string()
            .contains("non-finite finite-difference value [1]"));

        let e = optim_lbfgsb(
            &mut f,
            None,
            &[0.1, 0.1],
            &[0.0],
            &[1.0, 1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(e, PopmixError::Validation(_)));
    }

    #[test]
    fn test_optim_wrong_gradient_is_abnormal_termination() {
        // R: a gradient of the wrong sign makes every line search fail;
        // after 20 trials it stops with code 52 at the start, 21 evaluations
        let mut f = |x: &[f64]| Ok(quad(x));
        let mut gr = |x: &[f64]| Ok(vec![-2.0 * (x[0] - 0.3), -2.0 * (x[1] - 0.7)]);
        let r = optim_lbfgsb(
            &mut f,
            Some(&mut gr),
            &[0.1, 0.1],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &LbfgsbOptions::default(),
        )
        .unwrap();
        assert_eq!(r.convergence, 52);
        assert_eq!(r.message, "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH");
        assert_eq!(r.evaluations, 21);
        assert_eq!(r.par, vec![0.1, 0.1]);
    }
}
