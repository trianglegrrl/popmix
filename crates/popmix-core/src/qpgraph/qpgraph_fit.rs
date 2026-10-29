//! qpGraph as in admixtools2 (`qpgraph`, with its default `cpp = TRUE`).
//!
//! A graph is a list of (from, to) edges. Vertices are numbered in order of
//! first appearance, row by row; the leaves (no outgoing edge) in that order
//! are the populations, and the first is the f3 base population. A vertex
//! with two incoming edges is an admixture node: the incoming edge whose
//! source vertex comes first carries weight w, the other 1 - w (igraph's
//! `incident_edges` order).
//!
//! 1. f3(base; i, j) for every pair i <= j of the other leaves, from the
//!    leave-one-block-out f2 values, with their weighted block-jackknife
//!    covariance (Busing et al. 1999; `docs/cleanroom/jackknife.md`).
//!    `diag_f3 * trace` is added to the covariance diagonal before it is
//!    inverted.
//! 2. For given admixture weights, each drift edge's share of each leaf's
//!    ancestry (the weights of the paths from the root to the leaf through
//!    that edge) makes f3 linear in the edge lengths (Patterson et al. 2012).
//!    The lengths minimise the covariance-weighted f3 residual with
//!    `diag * mean(diag)` added to the normal equations, subject to
//!    lengths >= 0: a box-constrained quadratic program, solved exactly by
//!    popmix's own active-set solver (`docs/cleanroom/qpgraph-edge-lengths.md`).
//!    The score is the covariance-weighted squared residual.
//! 3. The admixture weights minimise the score over [0, 1] from `numstart`
//!    starts drawn as admixtools2 draws them after `set.seed(seed)` (R's
//!    uniforms, numstart × nadmix of them, start by start within each
//!    admixture weight), each with R's
//!    `optim(method = "L-BFGS-B")` and admixtools2's controls (`fnscale`
//!    1e-6, `maxit` 1e4, finite-difference gradients); the best start wins,
//!    the first among ties. Its weights are then polished past a zero-bound
//!    kink where L-BFGS-B stalls (`polish`, ADMIXTOOLS2-0009; not with
//!    `compat.qpgraph_stall_at_kink`). The starts' own results are kept.
//!
//! Provenance: Patterson et al. 2012 (Genetics, qpGraph), Maier et al. 2023
//! (eLife, admixtools2), admixtools2's documentation, the clean-room
//! specifications in `docs/cleanroom/`, and admixtools2's black-box outputs.

use nalgebra::{DMatrix, DVector};
use ndarray::Array2;

use crate::fstats::{jackknife_cov_from_loo, leave_one_out, F2Blocks};
use crate::qpadm::Admixtools2Compat;
use crate::rcompat::{optim_lbfgsb, LbfgsbOptions, RRng};
use crate::{PopmixError, Result};

/// R `.Machine$integer.max`: admixtools2's upper bound on drift edges.
const EDGE_UPPER: f64 = 2147483647.0;

/// Options for [`qpgraph`]; defaults are admixtools2's.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpGraphOptions {
    /// Regularization of the edge-length normal equations, scaled by the
    /// mean of their diagonal (admixtools2 `diag`, default 1e-4).
    pub diag: f64,
    /// Regularization of the f3 covariance, scaled by its trace (admixtools2
    /// `diag_f3`, default 1e-5).
    pub diag_f3: f64,
    /// Drift lengths >= 0 and admixture weights in [0, 1] (default true).
    pub constrained: bool,
    /// Random starts for the admixture weights (default 10).
    pub numstart: usize,
    /// Seed for the random starts, as R's `set.seed` (admixtools2 `seed`).
    pub seed: i32,
    /// Fit at these admixture weights (one per admixture node, the weight
    /// of its first incoming edge) instead of optimizing them.
    pub admix_weights: Option<Vec<f64>>,
    /// Also score the fit against these blocks' f3 statistics, weighted by
    /// the fitted blocks' inverse f3 covariance (admixtools2
    /// `f2_blocks_test`, for an out-of-sample score).
    pub f2_blocks_test: Option<F2Blocks>,
    /// admixtools2 2.0.10 quirks to reproduce (ADMIXTOOLS2-0009: stop where
    /// L-BFGS-B stalls, without polishing).
    pub compat: Admixtools2Compat,
}

impl Default for QpGraphOptions {
    fn default() -> Self {
        Self {
            diag: 1e-4,
            diag_f3: 1e-5,
            constrained: true,
            numstart: 10,
            seed: 1,
            admix_weights: None,
            f2_blocks_test: None,
            compat: Admixtools2Compat::CORRECTED,
        }
    }
}

/// One fitted edge, in the order the graph lists them.
#[derive(Debug, Clone)]
pub struct QpGraphEdge {
    pub from: String,
    pub to: String,
    /// Whether this is an admixture edge (into an admixture node).
    pub admix: bool,
    /// Drift length, or admixture weight for an admixture edge.
    pub weight: f64,
    /// Admixture edges: lowest and highest weight over the random starts.
    pub low: Option<f64>,
    pub high: Option<f64>,
}

/// One fitted f3(pop1; pop2, pop3).
#[derive(Debug, Clone)]
pub struct QpGraphF3 {
    pub pop1: String,
    pub pop2: String,
    pub pop3: String,
    pub est: f64,
    pub se: f64,
    pub fit: f64,
    /// est - fit
    pub diff: f64,
    /// diff / se
    pub z: f64,
    /// Two-sided normal p-value of z.
    pub p: f64,
}

/// One random start of the admixture weight optimization.
#[derive(Debug, Clone)]
pub struct QpGraphStart {
    pub start: Vec<f64>,
    pub weights: Vec<f64>,
    pub score: f64,
    /// Score evaluations with their gradient (admixtools2 `fevals`).
    pub evaluations: usize,
    /// R's `optim` convergence code: 0 converged, 1 iteration limit.
    pub convergence: i32,
}

/// Result of [`qpgraph`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct QpGraphResult {
    /// Leaves in graph order; the first is the f3 base population.
    pub pops: Vec<String>,
    pub edges: Vec<QpGraphEdge>,
    /// Weighted f3 residual; lower is better.
    pub score: f64,
    pub f3: Vec<QpGraphF3>,
    /// The random starts (empty without admixture or with fixed weights).
    pub opt: Vec<QpGraphStart>,
    /// Inverse of the regularized f3 covariance.
    pub ppinv: Array2<f64>,
    /// The score against `f2_blocks_test`'s f3 statistics, when given.
    pub score_test: Option<f64>,
}

/// An edge-list graph, numbered as igraph `graph_from_edgelist` numbers it.
#[derive(Debug)]
struct EdgeListGraph {
    vertices: Vec<String>,
    edges: Vec<(usize, usize)>,
    /// Leaf vertices, in vertex order.
    leaves: Vec<usize>,
    root: usize,
    /// Per admixture node (in vertex order): its (first, second) incoming edges.
    admix: Vec<[usize; 2]>,
    /// Edges that are not admixture edges, in edge order.
    normal: Vec<usize>,
}

impl EdgeListGraph {
    fn parse(edges: &[(String, String)]) -> Result<Self> {
        let invalid = |msg: String| PopmixError::InvalidGraph {
            details: format!("qpgraph: {msg}"),
        };
        let mut vertices: Vec<String> = Vec::new();
        let id = |name: &str, vertices: &mut Vec<String>| match vertices
            .iter()
            .position(|v| v == name)
        {
            Some(i) => i,
            None => {
                vertices.push(name.to_string());
                vertices.len() - 1
            }
        };
        let edges: Vec<(usize, usize)> = edges
            .iter()
            .map(|(a, b)| {
                let a = id(a, &mut vertices);
                (a, id(b, &mut vertices))
            })
            .collect();
        let n = vertices.len();
        let mut indeg = vec![0usize; n];
        let mut outdeg = vec![0usize; n];
        for &(a, b) in &edges {
            if a == b {
                return Err(invalid(format!("self-loop at {}", vertices[a])));
            }
            outdeg[a] += 1;
            indeg[b] += 1;
        }
        let roots: Vec<usize> = (0..n).filter(|&v| indeg[v] == 0).collect();
        if roots.len() != 1 {
            let names: Vec<&str> = roots.iter().map(|&v| vertices[v].as_str()).collect();
            return Err(invalid(format!(
                "graph needs exactly one root, found {names:?}"
            )));
        }
        if let Some(v) = (0..n).find(|&v| indeg[v] > 2) {
            return Err(invalid(format!(
                "{} has more than two incoming edges",
                vertices[v]
            )));
        }
        // Kahn's algorithm: every vertex must be reachable without a cycle
        let mut remaining = indeg.clone();
        let mut queue = roots.clone();
        let mut seen = 0;
        while let Some(v) = queue.pop() {
            seen += 1;
            for &(a, b) in &edges {
                if a == v {
                    remaining[b] -= 1;
                    if remaining[b] == 0 {
                        queue.push(b);
                    }
                }
            }
        }
        if seen != n {
            return Err(invalid("graph has a cycle".to_string()));
        }
        let admix: Vec<[usize; 2]> = (0..n)
            .filter(|&v| indeg[v] == 2)
            .map(|v| {
                let mut inc: Vec<usize> = (0..edges.len()).filter(|&e| edges[e].1 == v).collect();
                inc.sort_by_key(|&e| (edges[e].0, e));
                [inc[0], inc[1]]
            })
            .collect();
        let normal = (0..edges.len())
            .filter(|&e| indeg[edges[e].1] != 2)
            .collect();
        Ok(Self {
            leaves: (0..n).filter(|&v| outdeg[v] == 0).collect(),
            root: roots[0],
            vertices,
            edges,
            admix,
            normal,
        })
    }
}

/// One root-to-leaf path: its leaf (index into `leaves`), drift edges
/// (indices into `normal`) and admixture edges as (node, is first edge).
struct Path {
    leaf: usize,
    normal: Vec<usize>,
    admix: Vec<(usize, bool)>,
}

/// The paths of every leaf, and how many paths reach each leaf.
struct Paths {
    paths: Vec<Path>,
    numpaths: Vec<usize>,
    n_normal: usize,
}

impl Paths {
    fn new(g: &EdgeListGraph) -> Self {
        let normal_index: Vec<Option<usize>> = (0..g.edges.len())
            .map(|e| g.normal.iter().position(|&n| n == e))
            .collect();
        let admix_index: Vec<Option<(usize, bool)>> = (0..g.edges.len())
            .map(|e| {
                g.admix.iter().enumerate().find_map(|(k, pair)| {
                    (pair[0] == e || pair[1] == e).then_some((k, pair[0] == e))
                })
            })
            .collect();
        let mut paths = Vec::new();
        let mut stack: Vec<(usize, Vec<usize>)> = vec![(g.root, Vec::new())];
        while let Some((v, path)) = stack.pop() {
            if let Some(leaf) = g.leaves.iter().position(|&l| l == v) {
                paths.push(Path {
                    leaf,
                    normal: path.iter().filter_map(|&e| normal_index[e]).collect(),
                    admix: path.iter().filter_map(|&e| admix_index[e]).collect(),
                });
                continue;
            }
            for (e, &(a, b)) in g.edges.iter().enumerate().rev() {
                if a == v {
                    let mut p = path.clone();
                    p.push(e);
                    stack.push((b, p));
                }
            }
        }
        let mut numpaths = vec![0; g.leaves.len()];
        for p in &paths {
            numpaths[p.leaf] += 1;
        }
        Self {
            paths,
            numpaths,
            n_normal: g.normal.len(),
        }
    }

    /// Share of each leaf's ancestry through each drift edge: drift edges x
    /// leaves. Without admixture it is the fraction of the leaf's paths
    /// through the edge; with admixture, edges on only some of the paths to
    /// the leaf get the summed weights of those paths.
    fn pwts(&self, weights: &[f64]) -> DMatrix<f64> {
        let nleaf = self.numpaths.len();
        let mut base = DMatrix::zeros(self.n_normal, nleaf);
        let mut count = DMatrix::<usize>::zeros(self.n_normal, nleaf);
        let mut filled = DMatrix::zeros(self.n_normal, nleaf);
        for p in &self.paths {
            let w: f64 = p.admix.iter().fold(1.0, |acc, &(k, first)| {
                acc * if first { weights[k] } else { 1.0 - weights[k] }
            });
            for &e in &p.normal {
                base[(e, p.leaf)] += 1.0 / self.numpaths[p.leaf] as f64;
                count[(e, p.leaf)] += 1;
                filled[(e, p.leaf)] += w;
            }
        }
        if weights.is_empty() {
            return base;
        }
        // only the edges on some but not all of the paths to the leaf take
        // the path weights; the others keep the unweighted share
        DMatrix::from_fn(self.n_normal, nleaf, |e, l| {
            if count[(e, l)] > 0 && count[(e, l)] < self.numpaths[l] {
                filled[(e, l)]
            } else {
                base[(e, l)]
            }
        })
    }
}

/// The f3 pairs (i, j), 1 <= i <= j < n, in admixtools2 order.
fn f3_pairs(n: usize) -> Vec<(usize, usize)> {
    (1..n).flat_map(|i| (i..n).map(move |j| (i, j))).collect()
}

/// f3 rows x drift edges: f3(base; i, j) = sum_e len_e d_ei d_ej, with
/// d = pwts - pwts[, base].
fn f3_design(pwts: &DMatrix<f64>) -> DMatrix<f64> {
    let pairs = f3_pairs(pwts.ncols());
    DMatrix::from_fn(pairs.len(), pwts.nrows(), |r, e| {
        let (i, j) = pairs[r];
        (pwts[(e, i)] - pwts[(e, 0)]) * (pwts[(e, j)] - pwts[(e, 0)])
    })
}

/// f3 estimates, their covariance and the inverse of the regularized
/// covariance.
struct F3Stats {
    est: DVector<f64>,
    var: DMatrix<f64>,
    ppinv: DMatrix<f64>,
}

/// f3(P0; Pi, Pj) for the pairs of [`f3_pairs`], with `pops[0]` = P0.
///
/// Per block, f3(P0; Pi, Pj) = (f2(P0, Pi) + f2(P0, Pj) - f2(Pi, Pj)) / 2
/// (so f3(P0; Pi, Pi) = f2(P0, Pi)). f3 is linear in f2, so its
/// leave-one-block-out values are this combination of the f2
/// leave-one-block-out values. The estimates and covariance are the
/// weighted block jackknife of these block-mean statistics over the blocks
/// where all are finite, the covariance of `docs/cleanroom/jackknife.md`'s
/// `block_jackknife_cov` computed from the leave-one-out values
/// ([`jackknife_cov_from_loo`], which avoids cancellation). The covariance,
/// with `diag_f3` times its trace added to every diagonal element
/// (admixtools2's `?qpgraph`: "after scaling by the matrix trace"), is then
/// inverted; a singular matrix is an error.
fn precompute_f3(f2_blocks: &F2Blocks, pops: &[&str], diag_f3: f64) -> Result<F3Stats> {
    let idx: Vec<usize> = pops
        .iter()
        .map(|p| {
            f2_blocks
                .pop_index(p)
                .ok_or_else(|| PopmixError::PopulationNotFound(p.to_string()))
        })
        .collect::<Result<_>>()?;
    let lengths = f2_blocks.block_lengths_f64();
    let n = pops.len();
    // leave-one-block-out f2 of every pair of `pops`
    let mut f2_loo = vec![vec![Vec::new(); n]; n];
    for i in 0..n {
        for j in i..n {
            let means: Vec<f64> = (0..f2_blocks.n_blocks)
                .map(|b| f2_blocks.block_mean(b, idx[i], idx[j]))
                .collect();
            let loo = leave_one_out(&means, &lengths);
            f2_loo[j][i].clone_from(&loo);

            f2_loo[i][j] = loo;
        }
    }
    let loo: Vec<Vec<f64>> = f3_pairs(n)
        .into_iter()
        .map(|(i, j)| {
            (0..f2_blocks.n_blocks)
                .map(|b| (f2_loo[0][i][b] + f2_loo[0][j][b] - f2_loo[i][j][b]) / 2.0)
                .collect()
        })
        .collect();
    let (est, cov) = jackknife_cov_from_loo(&loo, &lengths).ok_or_else(|| {
        PopmixError::data("qpgraph: fewer than two blocks with finite f3 statistics")
    })?;
    let k = est.len();
    let var = DMatrix::from_fn(k, k, |a, b| cov[[a, b]]);
    let mut reg = var.clone();
    let add = diag_f3 * reg.trace();
    for a in 0..k {
        reg[(a, a)] += add;
    }
    let ppinv = reg
        .try_inverse()
        .ok_or_else(|| PopmixError::SingularMatrix {
            context: "qpgraph f3 covariance".to_string(),
        })?;
    Ok(F3Stats {
        est: DVector::from_vec(est),
        var,
        ppinv,
    })
}

/// Strictly convex box-constrained QP: minimise ½ x'Dx - q'x subject to
/// lo <= x <= hi, by a primal active-set method. Exact (up to rounding) in
/// finitely many steps; quadprog solves the same problem in admixtools2.
fn box_qp(
    d: &DMatrix<f64>,
    q: &DVector<f64>,
    lo: &DVector<f64>,
    hi: &DVector<f64>,
) -> Result<DVector<f64>> {
    #[derive(Clone, Copy, PartialEq)]
    enum State {
        Free,
        Lower,
        Upper,
    }
    let n = q.len();
    let singular = || PopmixError::SingularMatrix {
        context: "qpgraph edge-length system".to_string(),
    };
    let solve_free = |state: &[State], x: &DVector<f64>| -> Result<DVector<f64>> {
        let free: Vec<usize> = (0..n).filter(|&i| state[i] == State::Free).collect();
        let mut out = x.clone();
        if free.is_empty() {
            return Ok(out);
        }
        let dff = DMatrix::from_fn(free.len(), free.len(), |a, b| d[(free[a], free[b])]);
        let rhs = DVector::from_fn(free.len(), |a, _| {
            q[free[a]]
                - (0..n)
                    .filter(|&j| state[j] != State::Free)
                    .map(|j| d[(free[a], j)] * x[j])
                    .sum::<f64>()
        });
        let sol = dff.lu().solve(&rhs).ok_or_else(singular)?;
        for (a, &i) in free.iter().enumerate() {
            out[i] = sol[a];
        }
        Ok(out)
    };

    // Start from the unconstrained optimum, clamped into the box.
    let mut state = vec![State::Free; n];
    let mut x = solve_free(&state, &DVector::zeros(n))?;
    for i in 0..n {
        if x[i] <= lo[i] {
            x[i] = lo[i];
            state[i] = State::Lower;
        } else if x[i] >= hi[i] {
            x[i] = hi[i];
            state[i] = State::Upper;
        }
    }
    let gtol = 1e-13 * (1.0 + q.amax());
    let max_steps = 20 * n + 100;
    for _ in 0..max_steps {
        let target = solve_free(&state, &x)?;
        // largest step towards the subproblem optimum that stays in the box
        let mut alpha = 1.0;
        let mut blocking = None;
        for i in (0..n).filter(|&i| state[i] == State::Free) {
            let step = target[i] - x[i];
            let (bound, s) = if step < 0.0 {
                (lo[i], State::Lower)
            } else {
                (hi[i], State::Upper)
            };
            if (step < 0.0 && target[i] < lo[i]) || (step > 0.0 && target[i] > hi[i]) {
                let a = (bound - x[i]) / step;
                if a < alpha {
                    alpha = a;
                    blocking = Some((i, s));
                }
            }
        }
        if let Some((i, s)) = blocking {
            x += (target - &x) * alpha;
            x[i] = if s == State::Lower { lo[i] } else { hi[i] };
            state[i] = s;
            continue;
        }
        x = target;
        // KKT: a variable at a bound must be held there by its gradient
        let g = d * &x - q;
        let release = (0..n)
            .filter_map(|i| match state[i] {
                State::Lower if g[i] < -gtol => Some((i, -g[i])),
                State::Upper if g[i] > gtol => Some((i, g[i])),
                _ => None,
            })
            .max_by(|a, b| a.1.total_cmp(&b.1));
        match release {
            Some((i, _)) => state[i] = State::Free,
            None => return Ok(x),
        }
    }
    Err(PopmixError::ConvergenceFailure {
        iterations: max_steps,
        context: "qpgraph edge-length optimization".to_string(),
    })
}

/// Everything that is fixed once the f3 statistics are known.
struct Model {
    paths: Paths,
    f3: F3Stats,
    diag: f64,
    constrained: bool,
}

impl Model {
    /// Drift lengths at fixed admixture weights, given the design D
    /// (fitted f3 = D l).
    ///
    /// Minimises (y − Dl)ᵀP(y − Dl) + λ‖l‖² with λ = `diag` · tr(DᵀPD) / k
    /// (k edges): by the normal equations (DᵀPD + λI) l = DᵀPy when
    /// unconstrained, otherwise as a box QP with 0 <= l <= `EDGE_UPPER`,
    /// solved by [`box_qp`] in the variables z_i = s_i l_i, s_i the square
    /// root of the i-th diagonal element of DᵀPD + λI.
    ///
    /// Provenance: written from the clean-room specification
    /// `docs/cleanroom/qpgraph-edge-lengths.md`, after Patterson et al. 2012
    /// (Genetics, qpGraph) and Maier et al. 2023 (eLife, admixtools2), and
    /// checked against admixtools2's black-box outputs at fixed weights.
    fn edge_lengths(&self, design: &DMatrix<f64>) -> Result<DVector<f64>> {
        let k = design.ncols();
        let dtp = design.transpose() * &self.f3.ppinv;
        let mut a = &dtp * design;
        let b = &dtp * &self.f3.est;
        let lambda = self.diag * a.trace() / k as f64;
        for i in 0..k {
            a[(i, i)] += lambda;
        }
        if !self.constrained {
            return a.lu().solve(&b).ok_or_else(|| PopmixError::SingularMatrix {
                context: "qpgraph edge-length system".to_string(),
            });
        }
        let s = DVector::from_fn(k, |i, _| a[(i, i)].sqrt());
        let d = DMatrix::from_fn(k, k, |i, j| a[(i, j)] / (s[i] * s[j]));
        let q = DVector::from_fn(k, |i, _| b[i] / s[i]);
        let lo = DVector::zeros(k);
        let hi = s.map(|si| EDGE_UPPER * si);
        let z = box_qp(&d, &q, &lo, &hi)?;
        Ok(z.component_div(&s))
    }

    /// Edge lengths, fitted f3 and score at the given admixture weights.
    fn fit(&self, weights: &[f64]) -> Result<(DVector<f64>, DVector<f64>, f64)> {
        let design = f3_design(&self.paths.pwts(weights));
        let lengths = self.edge_lengths(&design)?;
        let fit = &design * &lengths;
        let res = &self.f3.est - &fit;
        let score = (res.transpose() * &self.f3.ppinv * &res)[(0, 0)];
        Ok((lengths, fit, score))
    }

    fn score(&self, weights: &[f64]) -> Result<f64> {
        Ok(self.fit(weights)?.2)
    }
}

/// Fit an admixture graph to f2 blocks, as admixtools2 `qpgraph` does.
///
/// `edges` lists the graph's (from, to) edges; see the module docs for how
/// populations and admixture weights are ordered.
pub fn qpgraph(
    f2_blocks: &F2Blocks,
    edges: &[(String, String)],
    options: &QpGraphOptions,
) -> Result<QpGraphResult> {
    let opts = options;
    let g = EdgeListGraph::parse(edges)?;
    let pops: Vec<String> = g.leaves.iter().map(|&v| g.vertices[v].clone()).collect();
    if pops.len() < 3 {
        return Err(PopmixError::too_few(
            "qpgraph: the graph",
            "leaves",
            3,
            pops.len(),
        ));
    }
    let pop_refs: Vec<&str> = pops.iter().map(std::string::String::as_str).collect();
    let model = Model {
        paths: Paths::new(&g),
        f3: precompute_f3(f2_blocks, &pop_refs, opts.diag_f3)?,
        diag: opts.diag,
        constrained: opts.constrained,
    };
    let nadmix = g.admix.len();

    let mut opt = Vec::new();
    let weights = match &opts.admix_weights {
        Some(w) => {
            if w.len() != nadmix {
                return Err(PopmixError::DimensionMismatch {
                    context: "qpgraph admixture weights".to_string(),
                    expected: format!("{nadmix}, one per admixture node"),
                    actual: w.len().to_string(),
                });
            }
            w.clone()
        }
        None if nadmix == 0 => Vec::new(),
        None => {
            let (lo, hi) = if opts.constrained {
                (vec![0.0; nadmix], vec![1.0; nadmix])
            } else {
                (vec![f64::NEG_INFINITY; nadmix], vec![f64::INFINITY; nadmix])
            };
            let numstart = opts.numstart.max(1);
            let u = RRng::set_seed(opts.seed).runif(numstart * nadmix);
            let control = LbfgsbOptions {
                fnscale: 1e-6,
                maxit: 10_000,
                ..LbfgsbOptions::default()
            };
            for i in 0..numstart {
                // R fills the start matrix column by column
                let start: Vec<f64> = (0..nadmix).map(|j| u[j * numstart + i]).collect();
                let mut score = |w: &[f64]| model.score(w);
                let r = optim_lbfgsb(&mut score, None, &start, &lo, &hi, &control)?;
                opt.push(QpGraphStart {
                    start,
                    weights: r.par,
                    score: r.value,
                    evaluations: r.evaluations,
                    convergence: r.convergence,
                });
            }
            let best = opt
                .iter()
                .fold(None::<&QpGraphStart>, |best, s| match best {
                    Some(b) if b.score <= s.score => Some(b),
                    _ => Some(s),
                })
                .map(|s| s.weights.clone())
                .unwrap_or_default();
            if opts.compat.qpgraph_stall_at_kink {
                best
            } else {
                let mut score = |w: &[f64]| model.score(w);
                super::polish::polish(&mut score, &best, &lo, &hi)?.0
            }
        }
    };

    let (lengths, fit, score) = model.fit(&weights)?;
    let mut out_edges: Vec<QpGraphEdge> = g
        .edges
        .iter()
        .map(|&(a, b)| QpGraphEdge {
            from: g.vertices[a].clone(),
            to: g.vertices[b].clone(),
            admix: false,
            weight: f64::NAN,
            low: None,
            high: None,
        })
        .collect();
    for (k, &e) in g.normal.iter().enumerate() {
        out_edges[e].weight = lengths[k];
    }
    for (k, pair) in g.admix.iter().enumerate() {
        let range = (!opt.is_empty()).then(|| {
            opt.iter()
                .map(|s| s.weights[k])
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), w| {
                    (lo.min(w), hi.max(w))
                })
        });
        for (j, &e) in pair.iter().enumerate() {
            let edge = &mut out_edges[e];
            edge.admix = true;
            edge.weight = if j == 0 { weights[k] } else { 1.0 - weights[k] };
            if let Some((lo, hi)) = range {
                (edge.low, edge.high) = if j == 0 {
                    (Some(lo), Some(hi))
                } else {
                    (Some(1.0 - hi), Some(1.0 - lo))
                };
            }
        }
    }

    let normal = statrs::distribution::Normal::new(0.0, 1.0).expect("standard normal");
    let f3 = f3_pairs(pops.len())
        .iter()
        .enumerate()
        .map(|(r, &(i, j))| {
            use statrs::distribution::ContinuousCDF;
            let est = model.f3.est[r];
            let se = model.f3.var[(r, r)].sqrt();
            let diff = est - fit[r];
            let z = diff / se;
            QpGraphF3 {
                pop1: pops[0].clone(),
                pop2: pops[i].clone(),
                pop3: pops[j].clone(),
                est,
                se,
                fit: fit[r],
                diff,
                z,
                p: 2.0 * normal.cdf(-z.abs()),
            }
        })
        .collect();

    let ppinv = &model.f3.ppinv;
    let score_test = match &opts.f2_blocks_test {
        Some(test) => {
            let res = precompute_f3(test, &pop_refs, opts.diag_f3)?.est - &fit;
            Some((res.transpose() * ppinv * &res)[(0, 0)])
        }
        None => None,
    };
    Ok(QpGraphResult {
        pops,
        edges: out_edges,
        score,
        f3,
        opt,
        ppinv: Array2::from_shape_fn((ppinv.nrows(), ppinv.ncols()), |(i, j)| ppinv[(i, j)]),
        score_test,
    })
}

#[cfg(test)]
mod tests;
