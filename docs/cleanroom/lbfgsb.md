# Clean-room specification: L-BFGS-B with R's `optim` behaviour

`crates/popmix-core/src/rcompat/lbfgsb.rs`: `optim_lbfgsb`, which
qpgraph uses so its admixture-weight optimisation takes the same steps as
admixtools2's, which calls R's `optim(method = "L-BFGS-B")`.

The earlier implementation followed R's C code (GPL-2+). The new one is
built on:

- the L-BFGS-B 3.0 Fortran of Zhu, Byrd, Lu and Nocedal (with Morales and
  Nocedal's 2011 correction), BSD-3-licensed, from its distribution
  (https://users.iems.northwestern.edu/~nocedal/lbfgsb.html; SciPy's
  `scipy/optimize/lbfgsb_src` up to 1.14 also carries it), and its papers:
  Byrd, Lu, Nocedal, Zhu 1995 (SIAM J. Sci. Comput. 16:1190); Zhu, Byrd,
  Lu, Nocedal 1997 (ACM TOMS 23:550); Morales, Nocedal 2011 (ACM TOMS
  38:7); the line search is Moré and Thuente 1994 (ACM TOMS 20:286), as in
  MINPACK-2's `dcsrch`/`dcstep`;
- R's documentation of `optim` (`?optim`) for the driver's behaviour;
- R 4.5.2 run as a black box: `tests/golden/rlang/generate_optim.R` writes
  `tests/golden/rlang/r-4.5.2/optim.json`, and
  `crates/popmix-core/tests/rcompat_golden.rs` checks that popmix evaluates
  the objective the same number of times, at the same points, and stops
  with the same convergence code and message. You may write more R
  experiments of the same kind (R's *behaviour* is fair to observe; its
  source is not to be read).

Do not read R's source (`lbfgsb.c`, `optim.c`) or the earlier
`lbfgsb.rs`.

## What R runs, behaviourally

R's L-BFGS-B is version 2.3 of the algorithm (its documentation says
"the code for method L-BFGS-B is based on Fortran code by Zhu, Byrd, Lu-Chen
and Nocedal obtained from Netlib", version 2.3). Version 3.0 differs from
2.3 in two ways that change the steps taken, so to match R:

1. **Subspace minimisation:** 3.0 (Morales-Nocedal) projects the subspace
   step back into the box and backtracks to the feasible region along the
   path; 2.3 instead truncates the subspace step at the first bound it
   meets (the largest step length in [0, 1] keeping every free variable in
   its bounds, as described in Byrd et al. 1995, section 5). Implement
   2.3's truncation.
2. **Machine precision:** use `f64::EPSILON` (2.3 and R use the double
   precision epsilon; 3.0 computes its own).

Everything else (generalised Cauchy point, the compact limited-memory
matrices, the BFGS update with its curvature check, the Moré-Thuente line
search with `ftol = 1e-3`, `gtol = 0.9`, `xtol = 0.1`, the initial step
length of 1/‖d‖ on the first iteration and 1 thereafter, bounded by the box)
is the same in both; follow 3.0's code and the papers.

## The driver (`optim_lbfgsb`)

The public API (keep it exactly; the callers are
`crates/popmix-core/src/qpgraph/qpgraph_fit.rs` and
`crates/popmix-python/src/rcompat.rs`):

```text
pub struct LbfgsbOptions {      // Debug, Clone; Default gives R's defaults
    pub lmm: usize,             // correction pairs kept (5)
    pub factr: f64,             // relative reduction tolerance in units of epsilon (1e7)
    pub pgtol: f64,             // projected gradient tolerance (0)
    pub maxit: usize,           // iteration limit (100)
    pub fnscale: f64,           // the objective is divided by this (1)
    pub ndeps: f64,             // finite-difference step, every variable (1e-3)
}
pub struct LbfgsbResult {       // Debug, Clone
    pub par: Vec<f64>,
    pub value: f64,             // objective at par, not divided by fnscale
    pub evaluations: usize,     // objective(-and-gradient) evaluations: R's counts
    pub convergence: i32,       // 0 converged, 1 iteration limit, 51 warning, 52 error
    pub message: String,
}
pub type Gradient<'a> = &'a mut dyn FnMut(&[f64]) -> Result<Vec<f64>>;
pub fn optim_lbfgsb(
    f: &mut dyn FnMut(&[f64]) -> Result<f64>,
    gr: Option<Gradient>,       // None: R's finite differences
    par: &[f64], lower: &[f64], upper: &[f64],   // infinite bounds: unbounded side
    opts: &LbfgsbOptions,
) -> Result<LbfgsbResult>
```

Bounds of the wrong length are a `PopmixError::Validation` error.

- Without an analytic gradient, R uses central differences with step
  `ndeps` (in the scaled parameters) — confirm the exact rule at the bounds
  by experiment (`?optim`: "for L-BFGS-B ... [the finite-difference
  gradient] ... respects the bounds"); the goldens include bounded starts.
- Convergence codes and messages as R returns them (the goldens record the
  messages, e.g. "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",
  "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"); code 1 when the
  iteration limit is reached; a non-finite objective or gradient value is
  an error with the message R gives (observe it).
- The relative-reduction test is (f_old − f) / max(|f_old|, |f|, 1) <=
  factr · epsilon.

The number and order of objective evaluations must match R exactly; the
goldens are the oracle.
