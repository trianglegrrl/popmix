//! The line search of Moré and Thuente (1994), MINPACK-2's `dcsrch` and
//! `dcstep`.

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum SearchTask {
    Start,
    Fg,
    Convergence,
    Warning,
    Error,
}

/// The state dcsrch keeps between calls.
#[derive(Debug, Clone, Default)]
pub(super) struct SearchState {
    pub(super) brackt: bool,
    pub(super) stage: u8,
    pub(super) ginit: f64,
    pub(super) gtest: f64,
    pub(super) gx: f64,
    pub(super) gy: f64,
    pub(super) finit: f64,
    pub(super) fx: f64,
    pub(super) fy: f64,
    pub(super) stx: f64,
    pub(super) sty: f64,
    pub(super) stmin: f64,
    pub(super) stmax: f64,
    pub(super) width: f64,
    pub(super) width1: f64,
}

/// Moré-Thuente line search (MINPACK-2 `dcsrch`): called with the function
/// value `f` and derivative `g` at `stp`, it updates `stp` and sets `task`
/// to `Fg` to ask for the next evaluation, or to say it has finished.
#[allow(clippy::too_many_arguments)]
pub(super) fn dcsrch(
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
pub(super) fn dcstep(
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
