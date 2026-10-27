//! Admixture graphs for Python: building, fitting (qpGraph) and
//! searching (find_graphs).

use crate::compute;
use crate::f2::PyF2Blocks;
use popmix_core::python::input_error;
use pyo3::prelude::*;

/// An admixture graph: populations as leaves, joined by drift edges and
/// admixture nodes. Made by ``build_tree``.
#[pyclass(name = "AdmixGraph", module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAdmixGraph {
    inner: popmix_core::AdmixGraph,
}

#[pymethods]
impl PyAdmixGraph {
    /// Get the number of nodes.
    #[getter]
    fn node_count(&self) -> usize {
        self.inner.node_count()
    }

    /// Get the number of edges.
    #[getter]
    fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }

    /// Get leaf population names.
    #[getter]
    fn leaves(&self) -> Vec<String> {
        self.inner.leaf_names()
    }

    /// Get number of drift edges.
    #[getter]
    fn num_drift_edges(&self) -> usize {
        self.inner.num_drift_edges()
    }

    /// Get the weights of every edge.
    ///
    /// Returns
    /// -------
    /// list of float
    ///     One weight per edge (drift length, or admixture weight), in the
    ///     order the edges were added.
    fn get_edge_weights(&self) -> Vec<f64> {
        self.inner.get_edge_weights()
    }

    /// Compute the f2 between two populations that the graph predicts.
    ///
    /// The drift on the edges not shared by their paths from the root,
    /// summed over every pair of paths weighted by the admixture weights.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     Names of two nodes of the graph. An unknown name raises
    ///     ``PopulationNotFoundError``.
    ///
    /// Returns
    /// -------
    /// float
    ///     The expected f2.
    fn expected_f2(&self, pop_a: &str, pop_b: &str) -> PyResult<f64> {
        self.inner.expected_f2(pop_a, pop_b).map_err(PyErr::from)
    }
}

/// Build a tree (a graph with no admixture) from its nodes.
///
/// The tree is not validated.
///
/// Parameters
/// ----------
/// root_name
///     Name of the root node.
/// leaves
///     ``(name, parent, drift)`` tuples, one per leaf: the leaf's name, the
///     name of its parent node and the drift length of the edge to it.
/// internals
///     ``(name, parent, drift)`` tuples for the internal nodes, added before
///     the leaves, or None (the default) for none.
///
/// Returns
/// -------
/// AdmixGraph
///     The tree.
#[pyfunction]
#[pyo3(signature = (root_name, leaves, *, internals=None))]
pub(crate) fn build_tree(
    py: Python<'_>,
    root_name: &str,
    leaves: Vec<(String, String, f64)>,
    internals: Option<Vec<(String, String, f64)>>,
) -> PyResult<PyAdmixGraph> {
    compute(py, move || {
        let mut builder = popmix_core::GraphBuilder::new().root(root_name);

        // Add internal nodes first
        if let Some(ints) = internals {
            for (name, parent, drift) in ints {
                builder = builder.internal(&name, &parent, drift);
            }
        }

        // Add leaves
        for (name, parent, drift) in leaves {
            builder = builder.leaf(&name, &parent, drift);
        }

        let graph = builder.build_unchecked();

        Ok(PyAdmixGraph { inner: graph })
    })
}

/// One graph found by find_graphs.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyScoredGraph {
    /// The graph as (from, to) pairs; the path from the root to the first
    /// population comes first, so it is the f3 base.
    #[pyo3(get)]
    edges: Vec<(String, String)>,
    /// qpgraph score.
    #[pyo3(get)]
    score: f64,
}

#[pymethods]
impl PyScoredGraph {
    fn __repr__(&self) -> String {
        format!(
            "ScoredGraph(score={:.6}, {} edges)",
            self.score,
            self.edges.len()
        )
    }
}

/// Search for the admixture graphs that best fit the data.
///
/// Every candidate graph has the populations as leaves and ``nadmix``
/// admixture nodes and is fitted with qpgraph(); ``populations[0]`` is the f3
/// base population of every fit.
///
/// Parameters
/// ----------
/// f2_blocks
///     F2 blocks containing every population.
/// populations
///     The graph's leaves.
/// nadmix
///     Number of admixture nodes (default 0).
/// outgroup
///     Only graphs whose root has this population as a direct child, or
///     None (the default) for any.
/// method
///     "stepwise" (default): add the populations, then the admixture edges,
///     one at a time, keeping the ``beam_width`` best graphs after each
///     step. "exhaustive": every graph admixtools2 generate_all_graphs()
///     returns (feasible for about 5 populations with 1 admixture, or 4
///     with 2). Anything else raises ``InputError``.
/// beam_width
///     Graphs kept after each stepwise step (default 20).
/// numstart, seed, diag, diag_f3, constrained, compat
///     As for qpgraph().
///
/// Returns
/// -------
/// tuple
///     ``(graphs, n_evaluated)``: the graphs of the last step, best first,
///     as a list of ``PyScoredGraph``, and the number of qpgraph fits made.
#[pyfunction]
#[pyo3(signature = (f2_blocks, populations, nadmix=0, outgroup=None, method="stepwise", beam_width=20, numstart=10, seed=1, diag=1e-4, diag_f3=1e-5, constrained=true, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_graphs(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    populations: Vec<String>,
    nadmix: usize,
    outgroup: Option<String>,
    method: &str,
    beam_width: usize,
    numstart: usize,
    seed: i32,
    diag: f64,
    diag_f3: f64,
    constrained: bool,
    compat: Option<&str>,
) -> PyResult<(Vec<PyScoredGraph>, usize)> {
    let compat = crate::qpadm::admixtools2_compat(compat)?;
    compute(py, move || {
        let method = match method {
            "stepwise" => popmix_core::GraphSearchMethod::Stepwise { beam_width },
            "exhaustive" => popmix_core::GraphSearchMethod::Exhaustive,
            other => {
                return Err(input_error(format!(
                    "method must be \"stepwise\" or \"exhaustive\", got {other:?}"
                )))
            }
        };
        let opts = {
            let mut o = popmix_core::GraphSearchOptions::default();
            o.nadmix = nadmix;
            o.outgroup = outgroup;
            o.method = method;
            o.qpgraph = {
                let mut o = popmix_core::QpGraphOptions::default();
                o.diag = diag;
                o.diag_f3 = diag_f3;
                o.constrained = constrained;
                o.numstart = numstart;
                o.seed = seed;
                o.admix_weights = None;
                o.f2_blocks_test = None;
                o.compat = compat;
                o
            };
            o
        };
        let pop_refs: Vec<&str> = populations
            .iter()
            .map(std::string::String::as_str)
            .collect();
        let result = popmix_core::find_graphs(&f2_blocks.inner, &pop_refs, &opts)?;
        let graphs = result
            .graphs
            .into_iter()
            .map(|g| PyScoredGraph {
                edges: g.edges,
                score: g.score,
            })
            .collect();
        Ok((graphs, result.n_evaluated))
    })
}

/// Every admixture graph with these leaves and ``nadmix`` admixture nodes.
///
/// The same set as admixtools2 generate_all_graphs() returns, in its
/// order.
///
/// Parameters
/// ----------
/// leaves
///     The leaf populations.
/// nadmix
///     Number of admixture nodes (default 0).
///
/// Returns
/// -------
/// list of list of tuple
///     Each graph as a list of ``(from, to)`` edges whose path from the root
///     to ``leaves[0]`` comes first.
#[pyfunction]
#[pyo3(signature = (leaves, *, nadmix=0))]
pub(crate) fn all_graphs(
    py: Python<'_>,
    leaves: Vec<String>,
    nadmix: usize,
) -> PyResult<Vec<Vec<(String, String)>>> {
    compute(py, move || {
        let refs: Vec<&str> = leaves.iter().map(std::string::String::as_str).collect();
        popmix_core::all_graphs(&refs, nadmix).map_err(PyErr::from)
    })
}

/// The value admixtools2 graph_hash() digests for a simplified graph.
///
/// Two graphs are the same graph to admixtools2 exactly when their keys are
/// equal.
///
/// Parameters
/// ----------
/// edges
///     The graph as ``(from, to)`` pairs, already simplified (no node with
///     one parent and one child), as ``all_graphs`` and ``find_graphs``
///     return them. An invalid or unsimplified graph raises ``InputError``.
///
/// Returns
/// -------
/// list of str
///     The graph's key.
#[pyfunction]
pub(crate) fn graph_key(edges: Vec<(String, String)>) -> PyResult<Vec<String>> {
    popmix_core::graph_key(&edges).map_err(PyErr::from)
}

/// One fitted qpgraph edge.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpGraphEdge {
    #[pyo3(get)]
    from_node: String,
    #[pyo3(get)]
    to_node: String,
    /// "edge" (drift) or "admix", as in admixtools2.
    #[pyo3(get)]
    kind: String,
    /// Drift length, or admixture weight.
    #[pyo3(get)]
    weight: f64,
    /// Admixture edges: lowest and highest weight over the random starts.
    #[pyo3(get)]
    low: Option<f64>,
    #[pyo3(get)]
    high: Option<f64>,
}

/// One fitted f3(pop1; pop2, pop3) of a qpgraph fit.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpGraphF3 {
    #[pyo3(get)]
    pop1: String,
    #[pyo3(get)]
    pop2: String,
    #[pyo3(get)]
    pop3: String,
    #[pyo3(get)]
    est: f64,
    #[pyo3(get)]
    se: f64,
    #[pyo3(get)]
    fit: f64,
    #[pyo3(get)]
    diff: f64,
    #[pyo3(get)]
    z: f64,
    #[pyo3(get)]
    p: f64,
}

/// One random start of the qpgraph admixture weight optimization.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpGraphStart {
    #[pyo3(get)]
    start: Vec<f64>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    score: f64,
    /// Score evaluations (admixtools2 ``fevals``)
    #[pyo3(get)]
    evaluations: usize,
    /// R ``optim`` convergence code: 0 converged, 1 iteration limit
    #[pyo3(get)]
    convergence: i32,
}

#[pymethods]
impl PyQpGraphStart {
    /// Whether this start converged (``convergence == 0``)
    #[getter]
    fn converged(&self) -> bool {
        self.convergence == 0
    }
}

/// Result of qpgraph().
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpGraphResult {
    /// Leaves in graph order; the first is the f3 base population.
    #[pyo3(get)]
    pops: Vec<String>,
    #[pyo3(get)]
    score: f64,
    #[pyo3(get)]
    edges: Vec<PyQpGraphEdge>,
    #[pyo3(get)]
    f3: Vec<PyQpGraphF3>,
    #[pyo3(get)]
    opt: Vec<PyQpGraphStart>,
    /// Inverse of the regularized f3 covariance, rows in ``f3`` order.
    #[pyo3(get)]
    ppinv: Vec<Vec<f64>>,
    /// Score against ``f2_blocks_test``, when given.
    #[pyo3(get)]
    score_test: Option<f64>,
}

#[pymethods]
impl PyQpGraphResult {
    fn __repr__(&self) -> String {
        format!(
            "QpGraphResult(score={:.4}, pops={}, edges={})",
            self.score,
            self.pops.len(),
            self.edges.len()
        )
    }
}

/// Fit an admixture graph, as admixtools2 ``qpgraph`` does.
///
/// Drift lengths minimise the covariance-weighted residual of f3(base; i, j)
/// for all pairs of leaves, with lengths >= 0; admixture weights are optimized
/// over [0, 1] from ``numstart`` random starts. The populations are the
/// graph's leaves in order of first appearance in ``edges``; the first is the
/// f3 base population.
///
/// Parameters
/// ----------
/// f2_blocks
///     F2 blocks containing every leaf.
/// edges
///     The graph as ``(from, to)`` pairs. A node with two incoming edges is
///     an admixture node.
/// admix_weights
///     Fit at these weights instead of optimizing them, one per admixture
///     node in order of first appearance: the weight of the incoming edge
///     whose source node appears first. None (the default) to optimize.
/// numstart
///     Random starts for the admixture weights (default 10).
/// seed
///     Seed for the random starts, as R's ``set.seed`` (default 1); with
///     the same seed the starts are admixtools2's.
/// diag
///     Regularization of the edge-length fit (default 1e-4).
/// diag_f3
///     Regularization of the f3 covariance (default 1e-5).
/// constrained
///     Lengths >= 0 and weights in [0, 1] (default True).
/// f2_blocks_test
///     Also score the fit against these blocks' f3 statistics
///     (``score_test``), for an out-of-sample score; None (the default)
///     for none.
/// compat
///     None to polish the best start's weights past a zero-bound kink
///     where L-BFGS-B stalls; "admixtools2-2.0.10" to stop there, as
///     admixtools2 (ADMIXTOOLS2-0009).
///
/// Returns
/// -------
/// PyQpGraphResult
///     The fitted edges, score, f3 fits and the random starts.
#[pyfunction]
#[pyo3(signature = (f2_blocks, edges, *, admix_weights=None, numstart=10, seed=1, diag=1e-4, diag_f3=1e-5, constrained=true, f2_blocks_test=None, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpgraph(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    edges: Vec<(String, String)>,
    admix_weights: Option<Vec<f64>>,
    numstart: usize,
    seed: i32,
    diag: f64,
    diag_f3: f64,
    constrained: bool,
    f2_blocks_test: Option<&PyF2Blocks>,
    compat: Option<&str>,
) -> PyResult<PyQpGraphResult> {
    let compat = crate::qpadm::admixtools2_compat(compat)?;
    compute(py, move || {
        let opts = {
            let mut o = popmix_core::QpGraphOptions::default();
            o.diag = diag;
            o.diag_f3 = diag_f3;
            o.constrained = constrained;
            o.numstart = numstart;
            o.seed = seed;
            o.admix_weights = admix_weights;
            o.f2_blocks_test = f2_blocks_test.map(|b| b.inner.clone());
            o.compat = compat;
            o
        };
        let r = popmix_core::qpgraph(&f2_blocks.inner, &edges, &opts)?;
        Ok(PyQpGraphResult {
            pops: r.pops,
            score: r.score,
            edges: r
                .edges
                .into_iter()
                .map(|e| PyQpGraphEdge {
                    from_node: e.from,
                    to_node: e.to,
                    kind: if e.admix { "admix" } else { "edge" }.to_string(),
                    weight: e.weight,
                    low: e.low,
                    high: e.high,
                })
                .collect(),
            f3: r
                .f3
                .into_iter()
                .map(|f| PyQpGraphF3 {
                    pop1: f.pop1,
                    pop2: f.pop2,
                    pop3: f.pop3,
                    est: f.est,
                    se: f.se,
                    fit: f.fit,
                    diff: f.diff,
                    z: f.z,
                    p: f.p,
                })
                .collect(),
            opt: r
                .opt
                .into_iter()
                .map(|o| PyQpGraphStart {
                    start: o.start,
                    weights: o.weights,
                    score: o.score,
                    evaluations: o.evaluations,
                    convergence: o.convergence,
                })
                .collect(),
            ppinv: r.ppinv.outer_iter().map(|row| row.to_vec()).collect(),
            score_test: r.score_test,
        })
    })
}
