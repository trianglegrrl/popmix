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
    assert!(matches!(e, PopmixError::DimensionMismatch { .. }));
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
