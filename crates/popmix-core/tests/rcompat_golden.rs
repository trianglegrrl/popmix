//! popmix's R-compatible random numbers and L-BFGS-B against R itself.
//!
//! Reference: tests/golden/rlang/r-4.5.2/optim.json, from
//! tests/golden/rlang/generate_optim.R:
//! - runif: set.seed(seed); runif(1500), which must match exactly
//! - sample: set.seed(seed), then sample(n, size, replace = TRUE) repeatedly,
//!   exactly
//!
//! Numbers are "%.17g" strings, so they are exact.
//! - lbfgsb: optim(method = "L-BFGS-B") on bounded test problems; popmix
//!   must evaluate the objective the same number of times at the same
//!   points (see `tolerances`) and stop where R stops, with the same
//!   convergence code and message

use popmix_core::rcompat::{optim_lbfgsb, LbfgsbOptions, RRng};
use serde_json::Value;
use std::fs;
use std::path::PathBuf;

fn reference() -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/rlang/r-4.5.2/optim.json");
    serde_json::from_str(&fs::read_to_string(p).unwrap()).unwrap()
}

fn num(v: &Value) -> f64 {
    match v {
        Value::String(s) => s.parse().unwrap(),
        _ => v.as_f64().unwrap(),
    }
}

fn nums(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(num).collect()
}

fn rosen(x: &[f64]) -> f64 {
    (0..x.len() - 1)
        .map(|i| 100.0 * (x[i + 1] - x[i] * x[i]).powi(2) + (1.0 - x[i]).powi(2))
        .sum()
}

fn rosen_gr(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut g = vec![0.0; n];
    for i in 0..n - 1 {
        g[i] = -400.0 * x[i] * (x[i + 1] - x[i] * x[i]) - 2.0 * (1.0 - x[i]);
    }
    for i in 1..n {
        g[i] += 200.0 * (x[i] - x[i - 1] * x[i - 1]);
    }
    g
}

const CENTRE: [f64; 4] = [2.0, -1.0, 0.5, 3.0];

fn quad(x: &[f64]) -> f64 {
    x.iter()
        .enumerate()
        .map(|(i, &v)| (i + 1) as f64 * (v - CENTRE[i]).powi(2))
        .sum::<f64>()
        + x[0] * x[1]
}

fn quad_gr(x: &[f64]) -> Vec<f64> {
    let mut g: Vec<f64> = x
        .iter()
        .enumerate()
        .map(|(i, &v)| 2.0 * (i + 1) as f64 * (v - CENTRE[i]))
        .collect();
    g[0] += x[1];
    g[1] += x[0];
    g
}

fn bumpy(x: &[f64]) -> f64 {
    1e-4 * ((7.0 * x[0]).sin() * (5.0 * x[1]).cos()
        + (x[0] - 0.3).powi(2)
        + 2.0 * (x[1] - 0.6).powi(2)
        + 1.0)
}

/// How far popmix may be from R (relative to max(|x|, 1)): over the
/// evaluation points, and at the end. R's dot products and axpys come from
/// its BLAS (here OpenBLAS's SSE2 kernels, whose summation order depends on
/// memory alignment), so the two differ in the last bit from the first
/// update on, and Rosenbrock's valley amplifies that along the path.
fn tolerances(name: &str) -> (f64, f64) {
    match name {
        // finite differences of a function scaled up 1e6, in a narrow valley
        "rosen5_fd_scaled" => (1e-4, 1e-8),
        _ => (1e-8, 1e-10),
    }
}

#[test]
fn test_runif_matches_r() {
    for case in reference()["runif"].as_array().unwrap() {
        let seed = case["seed"].as_i64().unwrap() as i32;
        let want = nums(&case["values"]);
        assert_eq!(RRng::set_seed(seed).runif(want.len()), want, "seed {seed}");
    }
}

#[test]
fn test_sample_with_replacement_matches_r() {
    for case in reference()["sample"].as_array().unwrap() {
        let seed = case["seed"].as_i64().unwrap() as i32;
        let n = case["n"].as_u64().unwrap() as usize;
        let mut rng = RRng::set_seed(seed);
        for (rep, values) in case["values"].as_array().unwrap().iter().enumerate() {
            let want: Vec<usize> = values
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_u64().unwrap() as usize)
                .collect();
            let got: Vec<usize> = rng
                .sample_replace(n, want.len())
                .iter()
                .map(|i| i + 1)
                .collect();
            assert_eq!(got, want, "seed {seed}, n {n}, replicate {rep}");
        }
    }
}

#[test]
fn test_lbfgsb_matches_r_optim() {
    for case in reference()["lbfgsb"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let (f, g): (fn(&[f64]) -> f64, fn(&[f64]) -> Vec<f64>) = if name.starts_with("rosen") {
            (rosen, rosen_gr)
        } else if name.starts_with("quad") {
            (quad, quad_gr)
        } else {
            (bumpy, |_| unreachable!())
        };
        let control = &case["control"];
        let mut opts = LbfgsbOptions::default();
        if let Some(v) = control.get("fnscale") {
            opts.fnscale = num(v);
        }
        if let Some(v) = control.get("maxit") {
            opts.maxit = num(v) as usize;
        }
        if let Some(v) = control.get("lmm") {
            opts.lmm = num(v) as usize;
        }
        if let Some(v) = control.get("factr") {
            opts.factr = num(v);
        }
        let mut evals: Vec<Vec<f64>> = Vec::new();
        let mut obj = |x: &[f64]| {
            evals.push(x.to_vec());
            Ok(f(x))
        };
        let mut grad = |x: &[f64]| Ok(g(x));
        let gr: Option<&mut dyn FnMut(&[f64]) -> popmix_core::Result<Vec<f64>>> =
            if case["gradient"] == "analytic" {
                Some(&mut grad)
            } else {
                None
            };
        let res = optim_lbfgsb(
            &mut obj,
            gr,
            &nums(&case["par0"]),
            &nums(&case["lower"]),
            &nums(&case["upper"]),
            &opts,
        )
        .unwrap();

        let want_evals: Vec<Vec<f64>> =
            case["evals"].as_array().unwrap().iter().map(nums).collect();
        assert_eq!(
            evals.len(),
            want_evals.len(),
            "{name}: number of evaluations"
        );
        assert_eq!(
            res.evaluations as f64,
            nums(&case["counts"])[0],
            "{name}: counts"
        );
        assert_eq!(
            res.convergence as f64,
            num(&case["convergence"]),
            "{name}: convergence"
        );
        assert_eq!(
            res.message,
            case["message"].as_str().unwrap(),
            "{name}: message"
        );
        assert_eq!(evals[0], want_evals[0], "{name}: start");
        let (path_tol, par_tol) = tolerances(name);
        let dev = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .map(|(x, y)| (x - y).abs() / y.abs().max(1.0))
                .fold(0.0, f64::max)
        };
        let path = evals
            .iter()
            .zip(&want_evals)
            .map(|(a, b)| dev(a, b))
            .fold(0.0, f64::max);
        assert!(
            path <= path_tol,
            "{name}: evaluation points differ by {path:e}"
        );
        let par = dev(&res.par, &nums(&case["par"]));
        assert!(
            par <= par_tol,
            "{name}: par {:?} differs by {par:e}",
            res.par
        );
        let v = num(&case["value"]);
        assert!(
            (res.value - v).abs() <= par_tol * v.abs().max(1.0),
            "{name}: value {} vs {v}",
            res.value
        );
    }
}
