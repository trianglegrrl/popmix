//! Admixture graph search, scored by [`qpgraph`].
//!
//! Sources: `docs/cleanroom/small-helpers.md` section 5; admixtools2's
//! documentation of `generate_all_graphs()`, `graph_hash()` and
//! `find_graphs()`; and admixtools2 2.0.10 `generate_all_graphs()` run as a
//! black box: its output, graphs and their order, for the goldens of
//! `generate_all_graphs.R` and for 6 leaves with one admixture.
//!
//! Graphs are rooted DAGs: a root of out-degree 2, leaves of in-degree 1,
//! split vertices (in-degree 1, out-degree 2) and admixture vertices
//! (in-degree 2, out-degree 1). The exhaustive enumeration:
//! - trees: the first leaf below a stem vertex; each further leaf hung, in
//!   turn, on every edge of every tree so far; the stem dropped. Each rooted
//!   binary tree comes once.
//! - admixtures, one round per admixture vertex: for every graph of the
//!   previous round and every ordered pair of distinct edges from = (a, b),
//!   to = (c, d), a new vertex u on `from`, a new admixture vertex w on `to`
//!   and the edge u -> w, unless
//!   - d reaches a, so that the edge would close a cycle;
//!   - the two edges meet at a vertex (a = c, b = c or b = d): w's parents
//!     would be joined by an edge, a vertex expanded to a loop; or
//!   - either edge joins the root to a leaf.
//! - graphs with equal [`graph_key`] are kept once, the first time they
//!   appear, after each round.
//!
//! The documentation names only the loop ("a node is expanded to a loop");
//! the rule as stated, the root-to-leaf edges included, was established by
//! comparing every acyclic insertion with admixtools2's output, and with it
//! the enumeration returns the same graphs in the same order in every case
//! compared (4 leaves with 1 and 2 admixtures, 5 with 0 and 1, 6 with 1).
//! The edge layout follows admixtools2's too: a split edge is removed and
//! the new edges appended (for a leaf: to the leaf, into and out of the new
//! vertex; for an admixture: (a, u), (u, b), (c, w), (w, d), (u, w)). The
//! layout matters to fits: qpgraph's random starts follow the vertex order.
//!
//! The stepwise search adds the populations one at a time to the best trees
//! so far, then the admixture edges one at a time, keeping the
//! `beam_width` best graphs after each step.

use std::collections::HashSet;

use rayon::prelude::*;

use crate::fstats::F2Blocks;
use crate::{PopmixError, Result};

use super::qpgraph_fit::{qpgraph, QpGraphOptions};

/// How [`find_graphs`] chooses the graphs it fits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphSearchMethod {
    /// Every graph admixtools2 `generate_all_graphs()` returns. Their number
    /// grows very fast: 105 trees for 5 populations, 1485 graphs with one
    /// admixture, and 945 trees for 6.
    Exhaustive,
    /// Add populations, then admixture edges, one at a time, keeping the
    /// best `beam_width` graphs after each step.
    Stepwise { beam_width: usize },
}

/// Options for [`find_graphs`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct GraphSearchOptions {
    /// Number of admixture nodes in each graph.
    pub nadmix: usize,
    /// Only graphs whose root has this population as a direct child, as
    /// admixtools2 `find_graphs(outpop = )`.
    pub outgroup: Option<String>,
    pub method: GraphSearchMethod,
    /// How each graph is fitted.
    pub qpgraph: QpGraphOptions,
}

impl Default for GraphSearchOptions {
    fn default() -> Self {
        Self {
            nadmix: 0,
            outgroup: None,
            method: GraphSearchMethod::Stepwise { beam_width: 20 },
            qpgraph: QpGraphOptions::default(),
        }
    }
}

/// One fitted graph.
#[derive(Debug, Clone)]
pub struct ScoredGraph {
    /// The graph as qpgraph edges; the path from the root to the first
    /// population comes first, so that population is the f3 base.
    pub edges: Vec<(String, String)>,
    /// qpgraph score.
    pub score: f64,
}

/// Result of [`find_graphs`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct GraphSearchResult {
    /// The graphs fitted in the last step, best first (ties keep the order
    /// the graphs were generated in).
    pub graphs: Vec<ScoredGraph>,
    /// Number of qpgraph fits made.
    pub n_evaluated: usize,
}

/// A graph with vertices and edges in the order they were added (the layout
/// of the module doc). Internal vertices have no name.
#[derive(Debug, Clone)]
struct Topology {
    names: Vec<Option<String>>,
    edges: Vec<(usize, usize)>,
}

impl Topology {
    fn n(&self) -> usize {
        self.names.len()
    }

    fn degrees(&self) -> (Vec<usize>, Vec<usize>) {
        let mut indeg = vec![0; self.n()];
        let mut outdeg = vec![0; self.n()];
        for &(a, b) in &self.edges {
            outdeg[a] += 1;
            indeg[b] += 1;
        }
        (indeg, outdeg)
    }

    fn children(&self) -> Vec<Vec<usize>> {
        let mut out = vec![Vec::new(); self.n()];
        for &(a, b) in &self.edges {
            out[a].push(b);
        }
        out
    }

    fn root(&self) -> usize {
        let (indeg, _) = self.degrees();
        indeg
            .iter()
            .position(|&d| d == 0)
            .expect("a DAG has a root")
    }

    fn is_leaf(&self, outdeg: &[usize], v: usize) -> bool {
        outdeg[v] == 0
    }

    /// Shortest directed path lengths from `v`, `None` where unreachable.
    fn distances_from(&self, children: &[Vec<usize>], v: usize) -> Vec<Option<usize>> {
        let mut dist = vec![None; self.n()];
        dist[v] = Some(0);
        let mut queue = std::collections::VecDeque::from([v]);
        while let Some(u) = queue.pop_front() {
            let d = dist[u].unwrap() + 1;
            for &c in &children[u] {
                if dist[c].is_none() {
                    dist[c] = Some(d);
                    queue.push_back(c);
                }
            }
        }
        dist
    }

    /// Hang `leaf` from a new vertex on edge `e`. Edge `e` is removed and
    /// the edges to the leaf, into and out of the new vertex appended.
    fn insert_leaf(&self, e: usize, leaf: &str) -> Topology {
        let (a, b) = self.edges[e];
        let (m, l) = (self.n(), self.n() + 1);
        let mut names = self.names.clone();
        names.extend([None, Some(leaf.to_string())]);
        let mut edges = self.edges.clone();
        edges.remove(e);
        edges.extend([(m, l), (a, m), (m, b)]);
        Topology { names, edges }
    }

    /// One more admixture: a new vertex `u` on edge `from` = (a, b), a new
    /// admixture vertex `w` on edge `to` = (c, d) and the edge u -> w. Both
    /// edges are removed and (a, u), (u, b), (c, w), (w, d), (u, w)
    /// appended.
    fn add_admixture(&self, from: usize, to: usize) -> Topology {
        let ((a, b), (c, d)) = (self.edges[from], self.edges[to]);
        let (u, w) = (self.n(), self.n() + 1);
        let mut names = self.names.clone();
        names.extend([None, None]);
        let mut edges: Vec<(usize, usize)> = (0..self.edges.len())
            .filter(|&e| e != from && e != to)
            .map(|e| self.edges[e])
            .collect();
        edges.extend([(a, u), (u, b), (c, w), (w, d), (u, w)]);
        Topology { names, edges }
    }

    /// Every graph with one more admixture edge, for the ordered pairs of
    /// distinct edges (from, to) in order, except those the module doc
    /// excludes: `to` below `from` (a cycle), the two edges meeting at a
    /// vertex, or either joining the root to a leaf.
    fn plus_one(&self) -> Vec<Topology> {
        let children = self.children();
        let reach: Vec<Vec<Option<usize>>> = (0..self.n())
            .map(|v| self.distances_from(&children, v))
            .collect();
        let root = self.root();
        let root_to_leaf = |(x, y): (usize, usize)| x == root && children[y].is_empty();
        let mut out = Vec::new();
        for from in 0..self.edges.len() {
            let (a, b) = self.edges[from];
            for to in 0..self.edges.len() {
                let (c, d) = self.edges[to];
                let cycle = reach[d][a].is_some();
                let meet = a == c || b == c || b == d;
                if to == from
                    || cycle
                    || meet
                    || root_to_leaf(self.edges[from])
                    || root_to_leaf(self.edges[to])
                {
                    continue;
                }
                out.push(self.add_admixture(from, to));
            }
        }
        out
    }

    /// An admixture graph as the module doc describes: one root, of
    /// out-degree 2; leaves of in-degree 1; other vertices of in-degree 1
    /// and out-degree 2, or in-degree 2 and out-degree 1; no repeated edge
    /// and no cycle. Vertices with one parent and one child are allowed
    /// here ([`graph_key`] rejects them with its own message).
    fn is_valid(&self) -> bool {
        let (indeg, outdeg) = self.degrees();
        let roots: Vec<usize> = (0..self.n()).filter(|&v| indeg[v] == 0).collect();
        if roots.len() != 1 || outdeg[roots[0]] != 2 {
            return false;
        }
        let degrees_ok = (0..self.n())
            .filter(|&v| v != roots[0])
            .all(|v| matches!((indeg[v], outdeg[v]), (1, 0) | (1, 1) | (1, 2) | (2, 1)));
        let distinct: HashSet<&(usize, usize)> = self.edges.iter().collect();
        if !degrees_ok || distinct.len() != self.edges.len() {
            return false;
        }
        // acyclic: Kahn's algorithm removes every vertex
        let children = self.children();
        let mut remaining = indeg;
        let mut stack = vec![roots[0]];
        let mut removed = 0;
        while let Some(v) = stack.pop() {
            removed += 1;
            for &c in &children[v] {
                remaining[c] -= 1;
                if remaining[c] == 0 {
                    stack.push(c);
                }
            }
        }
        removed == self.n()
    }

    /// What admixtools2 `graph_hash()` digests: the sorted rows of distances
    /// from each internal vertex to each leaf ("." if unreachable), then the
    /// sorted leaf names. Rows are compared as bytes.
    fn key(&self) -> Vec<String> {
        let (_, outdeg) = self.degrees();
        let children = self.children();
        let mut leaves: Vec<(String, usize)> = (0..self.n())
            .filter(|&v| self.is_leaf(&outdeg, v))
            .map(|v| (self.names[v].clone().expect("leaves are named"), v))
            .collect();
        leaves.sort();
        let mut rows: Vec<String> = (0..self.n())
            .filter(|&v| !self.is_leaf(&outdeg, v))
            .map(|v| {
                let dist = self.distances_from(&children, v);
                leaves
                    .iter()
                    .map(|&(_, l)| dist[l].map_or(".".to_string(), |d| d.to_string()))
                    .collect()
            })
            .collect();
        rows.sort();
        rows.extend(leaves.into_iter().map(|(name, _)| name));
        rows
    }

    fn has_outgroup(&self, outgroup: &str) -> bool {
        let root = self.root();
        self.edges
            .iter()
            .any(|&(a, b)| a == root && self.names[b].as_deref() == Some(outgroup))
    }

    /// Named qpgraph edges, the path from the root to `base` first.
    /// Internal vertices are "R" (the root) and "n1", "n2", ... in vertex
    /// order, avoiding the leaf names.
    fn edge_list(&self, base: &str) -> Vec<(String, String)> {
        let leaf_names: HashSet<&str> = self.names.iter().flatten().map(|s| s.as_str()).collect();
        let root = self.root();
        let mut labels = Vec::with_capacity(self.n());
        let mut k = 0;
        for v in 0..self.n() {
            labels.push(match &self.names[v] {
                Some(name) => name.clone(),
                None if v == root => unique_label("R".to_string(), &leaf_names),
                None => {
                    k += 1;
                    unique_label(format!("n{k}"), &leaf_names)
                }
            });
        }
        // breadth-first path from the root to the base
        let children = self.children();
        let mut parent = vec![None; self.n()];
        let mut queue = std::collections::VecDeque::from([root]);
        let mut seen = vec![false; self.n()];
        seen[root] = true;
        while let Some(u) = queue.pop_front() {
            for &c in &children[u] {
                if !seen[c] {
                    seen[c] = true;
                    parent[c] = Some(u);
                    queue.push_back(c);
                }
            }
        }
        let base_v = (0..self.n())
            .find(|&v| self.names[v].as_deref() == Some(base))
            .expect("the base is a leaf of the graph");
        let mut path = Vec::new();
        let mut v = base_v;
        while let Some(p) = parent[v] {
            path.push((p, v));
            v = p;
        }
        path.reverse();
        let rest = self.edges.iter().filter(|e| !path.contains(e));
        path.iter()
            .chain(rest)
            .map(|&(a, b)| (labels[a].clone(), labels[b].clone()))
            .collect()
    }
}

fn unique_label(mut label: String, taken: &HashSet<&str>) -> String {
    while taken.contains(label.as_str()) {
        label.insert(0, '_');
    }
    label
}

fn check_leaves(leaves: &[&str]) -> Result<()> {
    if leaves.len() < 2 {
        return Err(PopmixError::too_few(
            "a graph",
            "populations",
            2,
            leaves.len(),
        ));
    }
    let mut seen = HashSet::new();
    if let Some(dup) = leaves.iter().find(|l| !seen.insert(**l)) {
        return Err(PopmixError::DuplicatePopulation(dup.to_string()));
    }
    Ok(())
}

/// Every rooted binary tree on `leaves` (at least two), each once:
/// (2n - 3)!! of them, in the module doc's order.
fn all_trees(leaves: &[&str]) -> Vec<Topology> {
    let mut trees = vec![Topology {
        names: vec![None, Some(leaves[0].to_string())],
        edges: vec![(0, 1)],
    }];
    for leaf in &leaves[1..] {
        trees = trees
            .iter()
            .flat_map(|t| (0..t.edges.len()).map(move |e| t.insert_leaf(e, leaf)))
            .collect();
    }
    for t in &mut trees {
        // the stem is vertex 0 and its only edge
        t.names.remove(0);
        t.edges = t
            .edges
            .iter()
            .filter(|&&(a, _)| a != 0)
            .map(|&(a, b)| (a - 1, b - 1))
            .collect();
    }
    trees
}

/// The graphs with distinct keys, each the first time it appears.
fn dedup(graphs: Vec<Topology>) -> Vec<Topology> {
    let mut seen = HashSet::new();
    graphs
        .into_iter()
        .filter(|g| seen.insert(g.key()))
        .collect()
}

/// The trees on `leaves`, then `nadmix` rounds of [`Topology::plus_one`]
/// over the graphs of the previous round, each round deduplicated.
fn enumerate(leaves: &[&str], nadmix: usize) -> Result<Vec<Topology>> {
    check_leaves(leaves)?;
    let mut graphs = all_trees(leaves);
    for _ in 0..nadmix {
        graphs = dedup(graphs.iter().flat_map(Topology::plus_one).collect());
    }
    Ok(graphs)
}

/// Every admixture graph with these leaves and `nadmix` admixture nodes, up
/// to [`graph_key`]: the same set as admixtools2 `generate_all_graphs()`,
/// in its order (see the module doc). Each graph is a
/// qpgraph edge list whose path from the root to `leaves[0]` comes first.
pub fn all_graphs(leaves: &[&str], nadmix: usize) -> Result<Vec<Vec<(String, String)>>> {
    Ok(enumerate(leaves, nadmix)?
        .iter()
        .map(|g| g.edge_list(leaves[0]))
        .collect())
}

/// The value admixtools2 `graph_hash()` digests for this graph: two graphs
/// are the same graph to admixtools2 exactly when their keys are equal.
///
/// admixtools2 simplifies a graph before hashing it; this takes graphs that
/// are already simplified (no node with one parent and one child), as
/// [`all_graphs`] and [`find_graphs`] return them.
pub fn graph_key(edges: &[(String, String)]) -> Result<Vec<String>> {
    let mut index = std::collections::HashMap::new();
    let mut names = Vec::new();
    let mut vertex = |name: &String| {
        *index.entry(name.clone()).or_insert_with(|| {
            names.push(name.clone());
            names.len() - 1
        })
    };
    let edge_ids: Vec<(usize, usize)> = edges.iter().map(|(a, b)| (vertex(a), vertex(b))).collect();
    let mut t = Topology {
        names: names.into_iter().map(Some).collect(),
        edges: edge_ids,
    };
    if !t.is_valid() {
        return Err(PopmixError::InvalidGraph {
            details: "not a valid admixture graph (admixtools2 is_valid)".to_string(),
        });
    }
    let (indeg, outdeg) = t.degrees();
    if (0..t.n()).any(|v| indeg[v] == 1 && outdeg[v] == 1) {
        return Err(PopmixError::InvalidGraph {
            details: "graph_key needs a simplified graph: no node with one parent and one child"
                .to_string(),
        });
    }
    for v in 0..t.n() {
        if !t.is_leaf(&outdeg, v) {
            t.names[v] = None;
        }
    }
    Ok(t.key())
}

/// Fit each graph with qpgraph; best first, ties in generation order.
type Ranked = Vec<(f64, usize)>;

/// A graph's score, or `None` for a graph that cannot be fit (its
/// edge-length system is singular or the fit does not converge); any other
/// error is in the data and stops the search.
fn candidate_score(fit: Result<f64>) -> Result<Option<f64>> {
    match fit {
        Ok(score) => Ok(Some(score)),
        Err(PopmixError::SingularMatrix { .. } | PopmixError::ConvergenceFailure { .. }) => {
            Ok(None)
        }
        Err(e) => Err(e),
    }
}

/// The graphs that can be fit, ranked; an error when none can.
fn score_all(
    f2_blocks: &F2Blocks,
    graphs: &[Topology],
    base: &str,
    options: &QpGraphOptions,
) -> Result<Ranked> {
    let scores: Vec<Option<f64>> = graphs
        .par_iter()
        .map(|g| candidate_score(qpgraph(f2_blocks, &g.edge_list(base), options).map(|r| r.score)))
        .collect::<Result<_>>()?;
    let mut ranked: Vec<(f64, usize)> = scores
        .into_iter()
        .zip(0..)
        .filter_map(|(score, i)| score.map(|s| (s, i)))
        .collect();
    if ranked.is_empty() {
        return Err(PopmixError::numerical(
            "find_graphs",
            format!("none of the {} candidate graphs could be fit", graphs.len()),
        ));
    }
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    Ok(ranked)
}

/// Search for the admixture graphs that best fit the f2 blocks.
///
/// Each graph has the populations as leaves and `nadmix` admixture nodes and
/// is fitted with [`qpgraph`]; `populations[0]` is the f3 base of every fit.
pub fn find_graphs(
    f2_blocks: &F2Blocks,
    populations: &[&str],
    options: &GraphSearchOptions,
) -> Result<GraphSearchResult> {
    let opts = options;
    check_leaves(populations)?;
    if populations.len() < 3 {
        return Err(PopmixError::too_few(
            "find_graphs",
            "populations",
            3,
            populations.len(),
        ));
    }
    for pop in populations {
        if !f2_blocks.populations.iter().any(|p| p.as_str() == *pop) {
            return Err(PopmixError::PopulationNotFound(pop.to_string()));
        }
    }
    if let Some(og) = &opts.outgroup {
        if !populations.contains(&og.as_str()) {
            return Err(PopmixError::parameter(
                "outgroup",
                format!("{og} is not one of the populations"),
            ));
        }
    }
    let base = populations[0];
    let keep = |graphs: Vec<Topology>| -> Vec<Topology> {
        match &opts.outgroup {
            Some(og) => graphs
                .into_iter()
                .filter(|g| !g.names.iter().flatten().any(|n| n == og) || g.has_outgroup(og))
                .collect(),
            None => graphs,
        }
    };
    let finish =
        |graphs: &[Topology], ranked: Vec<(f64, usize)>, n_evaluated: usize| GraphSearchResult {
            graphs: ranked
                .into_iter()
                .map(|(score, i)| ScoredGraph {
                    edges: graphs[i].edge_list(base),
                    score,
                })
                .collect(),
            n_evaluated,
        };

    match opts.method {
        GraphSearchMethod::Exhaustive => {
            let graphs = keep(enumerate(populations, opts.nadmix)?);
            let ranked = score_all(f2_blocks, &graphs, base, &opts.qpgraph)?;
            let n = graphs.len();
            Ok(finish(&graphs, ranked, n))
        }
        GraphSearchMethod::Stepwise { beam_width } => {
            if beam_width == 0 {
                return Err(PopmixError::parameter("beam_width", "must be at least 1"));
            }
            let mut n_evaluated = 0;
            let mut step = |candidates: Vec<Topology>| -> Result<(Vec<Topology>, Ranked)> {
                let candidates = keep(dedup(candidates));
                if candidates.is_empty() {
                    return Err(PopmixError::data(
                        "no graph satisfies the outgroup constraint",
                    ));
                }
                let ranked = score_all(f2_blocks, &candidates, base, &opts.qpgraph)?;
                n_evaluated += candidates.len();
                Ok((candidates, ranked))
            };
            let best = |graphs: &[Topology], ranked: &[(f64, usize)]| -> Vec<Topology> {
                ranked
                    .iter()
                    .take(beam_width)
                    .map(|&(_, i)| graphs[i].clone())
                    .collect()
            };

            let (mut graphs, mut ranked) = step(all_trees(&populations[..3]))?;
            for leaf in &populations[3..] {
                let candidates = best(&graphs, &ranked)
                    .iter()
                    .flat_map(|g| (0..g.edges.len()).map(move |e| g.insert_leaf(e, leaf)))
                    .collect();
                (graphs, ranked) = step(candidates)?;
            }
            for _ in 0..opts.nadmix {
                let mut candidates = Vec::new();
                for g in best(&graphs, &ranked) {
                    candidates.extend(g.plus_one());
                }
                (graphs, ranked) = step(candidates)?;
            }
            Ok(finish(&graphs, ranked, n_evaluated))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_graph_that_cannot_be_fit_is_passed_over() {
        let singular = PopmixError::SingularMatrix {
            context: "x".into(),
        };
        let stuck = PopmixError::ConvergenceFailure {
            iterations: 1,
            context: "x".into(),
        };
        assert_eq!(candidate_score(Ok(0.5)).unwrap(), Some(0.5));
        assert_eq!(candidate_score(Err(singular)).unwrap(), None);
        assert_eq!(candidate_score(Err(stuck)).unwrap(), None);
        // an error in the data stops the search
        assert!(candidate_score(Err(PopmixError::PopulationNotFound("A".into()))).is_err());
    }

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("P{i}")).collect()
    }

    #[test]
    fn tree_counts_are_the_rooted_binary_tree_counts() {
        for (n, count) in [(2, 1), (3, 3), (4, 15), (5, 105), (6, 945)] {
            let leaves = names(n);
            let refs: Vec<&str> = leaves.iter().map(|s| s.as_str()).collect();
            assert_eq!(all_graphs(&refs, 0).unwrap().len(), count, "{n} leaves");
        }
    }

    #[test]
    fn graph_counts_match_admixtools2_generate_all_graphs() {
        // admixtools2 2.0.10: length(generate_all_graphs(LETTERS[1:n], nadmix))
        for (n, nadmix, count) in [(4, 1, 84), (5, 1, 1485), (4, 2, 1302), (6, 1, 24615)] {
            let leaves = names(n);
            let refs: Vec<&str> = leaves.iter().map(|s| s.as_str()).collect();
            let graphs = all_graphs(&refs, nadmix).unwrap();
            assert_eq!(graphs.len(), count, "{n} leaves, {nadmix} admixture");
            for g in &graphs {
                let targets: HashSet<&String> = g.iter().map(|(_, b)| b).collect();
                let admix = targets
                    .iter()
                    .filter(|t| g.iter().filter(|(_, b)| b == **t).count() == 2)
                    .count();
                assert_eq!(admix, nadmix);
            }
        }
    }

    #[test]
    fn edge_list_puts_the_base_first() {
        let refs = ["A", "B", "C", "D"];
        for g in all_graphs(&refs, 1).unwrap() {
            let leaves: HashSet<&str> = refs.into_iter().collect();
            let first_leaf = g
                .iter()
                .flat_map(|(a, b)| [a.as_str(), b.as_str()])
                .find(|v| leaves.contains(v));
            assert_eq!(first_leaf, Some("A"));
            assert_eq!(g[0].0, "R");
        }
    }

    #[test]
    fn graph_key_is_invariant_to_internal_names_and_edge_order() {
        let e = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
            pairs
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect()
        };
        let one = e(&[("R", "A"), ("R", "x"), ("x", "B"), ("x", "C")]);
        let two = e(&[("x", "C"), ("root", "y"), ("y", "B"), ("root", "A")]);
        let two = two
            .into_iter()
            .map(|(a, b)| (a.replace('y', "x"), b.replace('y', "x")))
            .collect::<Vec<_>>();
        assert_eq!(graph_key(&one).unwrap(), graph_key(&two).unwrap());
        let other = e(&[("R", "B"), ("R", "x"), ("x", "A"), ("x", "C")]);
        assert_ne!(graph_key(&one).unwrap(), graph_key(&other).unwrap());
    }

    #[test]
    fn graph_key_rejects_an_invalid_graph() {
        let edges = vec![("R".to_string(), "A".to_string())];
        assert!(graph_key(&edges).is_err());
    }

    #[test]
    fn enumeration_rejects_duplicate_leaves() {
        assert!(all_graphs(&["A", "B", "A"], 0).is_err());
    }

    #[test]
    fn internal_labels_avoid_leaf_names() {
        let graphs = all_graphs(&["R", "n1", "n2"], 0).unwrap();
        for g in graphs {
            let mut seen = HashSet::new();
            let leaves = g
                .iter()
                .filter(|(_, b)| ["R", "n1", "n2"].contains(&b.as_str()));
            for (_, b) in leaves {
                assert!(seen.insert(b.clone()), "leaf {b} appears twice");
            }
            assert!(g
                .iter()
                .all(|(a, _)| !["R", "n1", "n2"].contains(&a.as_str())));
        }
    }

    /// A topology from named edges; vertices without children are leaves.
    fn topology<'a>(edges: &[(&'a str, &'a str)]) -> Topology {
        let mut names: Vec<&'a str> = Vec::new();
        let mut id = |v: &'a str| match names.iter().position(|n| *n == v) {
            Some(i) => i,
            None => {
                names.push(v);
                names.len() - 1
            }
        };
        let edges: Vec<(usize, usize)> = edges.iter().map(|(a, b)| (id(a), id(b))).collect();
        let parents: HashSet<usize> = edges.iter().map(|e| e.0).collect();
        Topology {
            names: (0..names.len())
                .map(|v| (!parents.contains(&v)).then(|| names[v].to_string()))
                .collect(),
            edges,
        }
    }

    #[test]
    fn admixture_edges_join_edges_that_do_not_meet_or_touch_a_root_leaf() {
        // ((A, B), (C, D)): a leaf edge of one side with a leaf edge of the
        // other (8), or a root edge with the other side's leaf edges (8)
        let balanced = topology(&[
            ("r", "x"),
            ("r", "y"),
            ("x", "A"),
            ("x", "B"),
            ("y", "C"),
            ("y", "D"),
        ]);
        assert_eq!(balanced.plus_one().len(), 16);
        // (A, (B, (C, D))): r -> A takes no part
        let caterpillar = topology(&[
            ("r", "A"),
            ("r", "x"),
            ("x", "B"),
            ("x", "y"),
            ("y", "C"),
            ("y", "D"),
        ]);
        let graphs = caterpillar.plus_one();
        assert_eq!(graphs.len(), 6);
        for g in &graphs {
            assert!(g.is_valid());
            // the new edges come last: (a, u), (u, b), (c, w), (w, d), (u, w)
            let (u, w) = (g.n() - 2, g.n() - 1);
            assert_eq!(g.edges[g.edges.len() - 1], (u, w));
            assert_eq!(g.edges[g.edges.len() - 3].1, w);
        }
    }

    #[test]
    fn every_insertion_keeps_the_graph_valid() {
        for t in all_trees(&["A", "B", "C", "D"]) {
            assert!(t.is_valid());
            for g in t.plus_one() {
                assert!(g.is_valid());
                assert!(g.plus_one().iter().all(|h| h.is_valid()));
            }
        }
    }
}
