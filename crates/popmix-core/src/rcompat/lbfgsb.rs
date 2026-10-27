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
//!
//! This file is the driver (`mainlb`); the steps are beside it: the
//! limited-memory matrices (`memory`), the Cauchy point (`cauchy`), the
//! subspace minimization (`subspace`), the line search (`line_search`) and
//! the LINPACK routines they use (`linalg`).

mod cauchy;

mod linalg;
mod line_search;
mod memory;
mod subspace;

use cauchy::*;
use linalg::*;
use line_search::*;
use memory::*;
use subspace::*;

use crate::error::{PopmixError, Result};

/// Settings of `optim(method = "L-BFGS-B")`; [`Default`] gives R's defaults.
#[derive(Debug, Clone)]
#[non_exhaustive]
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
#[non_exhaustive]
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
        return Err(PopmixError::DimensionMismatch {
            context: "optim_lbfgsb bounds".to_string(),
            expected: format!("{n} lower and {n} upper bounds"),
            actual: format!("{} and {}", lower.len(), upper.len()),
        });
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
            // R's own words, here and below
            return Err(PopmixError::numerical(
                "optim",
                "L-BFGS-B needs finite values of 'fn'",
            ));
        }
        let fx = raw / self.fnscale;
        let n = x.len();
        let g = match self.gr.as_mut() {
            Some(gr) => {
                let g = gr(x)?;
                if g.len() != n {
                    return Err(PopmixError::DimensionMismatch {
                        context: "the gradient in optim".to_string(),
                        expected: format!("length {n}"),
                        actual: format!("length {}", g.len()),
                    });
                }
                if g.iter().any(|v| !v.is_finite()) {
                    return Err(PopmixError::numerical(
                        "optim",
                        "non-finite value supplied by optim",
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
                        return Err(PopmixError::numerical(
                            "optim",
                            format!("non-finite finite-difference value [{}]", i + 1),
                        ));
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
mod tests;
