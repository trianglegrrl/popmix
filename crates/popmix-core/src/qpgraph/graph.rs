//! Admixture graph data structure.
//!
//! An admixture graph is a directed acyclic graph where:
//! - Leaf nodes represent sampled populations
//! - Internal nodes represent ancestral populations
//! - Drift edges represent genetic drift (with branch length)
//! - Admixture nodes have exactly two incoming edges with mixing proportions
//!
//! The graph flows from root (ancestor) to leaves (sampled populations).

use std::collections::{HashMap, HashSet};

use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;

use crate::{PopmixError, Result};

/// Unique identifier for a node in the graph.
pub type NodeId = NodeIndex;

/// Type of edge in the admixture graph.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeType {
    /// Normal drift edge with branch length.
    Drift,
    /// Admixture edge with mixing proportion (0 to 1).
    Admixture,
}

/// Edge data containing type and weight.
#[derive(Debug, Clone, Copy)]
pub struct EdgeData {
    /// Type of edge (drift or admixture).
    pub edge_type: EdgeType,
    /// Weight: branch length for drift, mixing proportion for admixture.
    pub weight: f64,
}

impl EdgeData {
    /// Create a new drift edge.
    pub fn drift(length: f64) -> Self {
        Self {
            edge_type: EdgeType::Drift,
            weight: length,
        }
    }

    /// Create a new admixture edge.
    pub fn admixture(proportion: f64) -> Self {
        Self {
            edge_type: EdgeType::Admixture,
            weight: proportion,
        }
    }
}

/// Node in the admixture graph.
#[derive(Debug, Clone)]
pub struct AdmixNode {
    /// Node name/label.
    pub name: String,
    /// Whether this is a leaf (sampled population).
    pub is_leaf: bool,
    /// Whether this is the root.
    pub is_root: bool,
    /// Whether this is an admixture node.
    pub is_admixed: bool,
}

impl AdmixNode {
    /// Create a leaf node (sampled population).
    pub fn leaf(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            is_leaf: true,
            is_root: false,
            is_admixed: false,
        }
    }

    /// Create an internal node.
    pub fn internal(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            is_leaf: false,
            is_root: false,
            is_admixed: false,
        }
    }

    /// Create the root node.
    pub fn root(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            is_leaf: false,
            is_root: true,
            is_admixed: false,
        }
    }

    /// Create an admixture node.
    pub fn admixed(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            is_leaf: false,
            is_root: false,
            is_admixed: true,
        }
    }
}

/// An admixture graph representing population relationships.
///
/// The graph is a directed acyclic graph (DAG) flowing from root to leaves.
/// Edges can be either drift edges (branch length) or admixture edges (mixing proportion).
#[derive(Debug, Clone)]
pub struct AdmixGraph {
    /// The underlying directed graph.
    graph: DiGraph<AdmixNode, EdgeData>,
    /// Map from population name to node index.
    name_to_node: HashMap<String, NodeId>,
    /// Root node index.
    root: Option<NodeId>,
}

impl AdmixGraph {
    /// Create a new empty admixture graph.
    pub fn new() -> Self {
        Self {
            graph: DiGraph::new(),
            name_to_node: HashMap::new(),
            root: None,
        }
    }

    /// Add a node to the graph.
    pub fn add_node(&mut self, node: AdmixNode) -> NodeId {
        let name = node.name.clone();
        let is_root = node.is_root;
        let idx = self.graph.add_node(node);
        self.name_to_node.insert(name, idx);
        if is_root {
            self.root = Some(idx);
        }
        idx
    }

    /// Add a drift edge from parent to child.
    pub fn add_drift_edge(&mut self, parent: NodeId, child: NodeId, length: f64) {
        self.graph.add_edge(parent, child, EdgeData::drift(length));
    }

    /// Add an admixture edge from source to admixed node.
    pub fn add_admix_edge(&mut self, source: NodeId, target: NodeId, proportion: f64) {
        self.graph
            .add_edge(source, target, EdgeData::admixture(proportion));
        // Mark target as admixed
        if let Some(node) = self.graph.node_weight_mut(target) {
            node.is_admixed = true;
        }
    }

    /// Get a node by name.
    pub fn get_node(&self, name: &str) -> Option<NodeId> {
        self.name_to_node.get(name).copied()
    }

    /// Get node data.
    pub fn node(&self, id: NodeId) -> Option<&AdmixNode> {
        self.graph.node_weight(id)
    }

    /// Get the root node.
    pub fn root(&self) -> Option<NodeId> {
        self.root
    }

    /// Get all leaf nodes.
    pub fn leaves(&self) -> Vec<NodeId> {
        self.graph
            .node_indices()
            .filter(|&idx| {
                self.graph
                    .node_weight(idx)
                    .map(|n| n.is_leaf)
                    .unwrap_or(false)
            })
            .collect()
    }

    /// Get leaf population names.
    pub fn leaf_names(&self) -> Vec<String> {
        self.leaves()
            .iter()
            .filter_map(|&idx| self.graph.node_weight(idx).map(|n| n.name.clone()))
            .collect()
    }

    /// Get number of nodes.
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get number of edges.
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Get incoming edges for a node.
    pub fn incoming_edges(&self, node: NodeId) -> Vec<(NodeId, EdgeData)> {
        self.graph
            .edges_directed(node, Direction::Incoming)
            .map(|e| (e.source(), *e.weight()))
            .collect()
    }

    /// Get outgoing edges for a node.
    pub fn outgoing_edges(&self, node: NodeId) -> Vec<(NodeId, EdgeData)> {
        self.graph
            .edges_directed(node, Direction::Outgoing)
            .map(|e| (e.target(), *e.weight()))
            .collect()
    }

    /// Check if a node is an admixture node (has 2 incoming edges).
    pub fn is_admixed(&self, node: NodeId) -> bool {
        self.graph.edges_directed(node, Direction::Incoming).count() == 2
    }

    /// Validate the graph structure.
    pub fn validate(&self) -> Result<()> {
        // Check for root
        if self.root.is_none() {
            return Err(PopmixError::InvalidGraph {
                details: "Graph has no root node".to_string(),
            });
        }

        let root = self.root.unwrap();

        // Root should have no incoming edges
        if self.graph.edges_directed(root, Direction::Incoming).count() > 0 {
            return Err(PopmixError::InvalidGraph {
                details: "Root node has incoming edges".to_string(),
            });
        }

        // Check for cycles using DFS
        let mut visited = HashSet::new();
        let mut rec_stack = HashSet::new();
        if self.has_cycle(root, &mut visited, &mut rec_stack) {
            return Err(PopmixError::InvalidGraph {
                details: "Graph contains a cycle".to_string(),
            });
        }

        // Check admixture nodes have exactly 2 incoming edges
        for idx in self.graph.node_indices() {
            let node = self.graph.node_weight(idx).unwrap();
            let in_count = self.graph.edges_directed(idx, Direction::Incoming).count();

            if node.is_admixed && in_count != 2 {
                return Err(PopmixError::InvalidAdmixtureNode {
                    node: node.name.clone(),
                });
            }
        }

        // Check admixture proportions sum to 1 for each admixed node
        for idx in self.graph.node_indices() {
            if self.is_admixed(idx) {
                let props: Vec<f64> = self
                    .graph
                    .edges_directed(idx, Direction::Incoming)
                    .filter_map(|e| {
                        let data = e.weight();
                        if data.edge_type == EdgeType::Admixture {
                            Some(data.weight)
                        } else {
                            None
                        }
                    })
                    .collect();

                if props.len() == 2 {
                    let sum = props[0] + props[1];
                    if (sum - 1.0).abs() > 1e-6 {
                        let name = self.graph.node_weight(idx).map(|n| &n.name).unwrap();
                        return Err(PopmixError::InvalidGraph {
                            details: format!(
                                "Admixture proportions for {} sum to {} (expected 1.0)",
                                name, sum
                            ),
                        });
                    }
                }
            }
        }

        Ok(())
    }

    /// Check for cycles using DFS.
    fn has_cycle(
        &self,
        node: NodeId,
        visited: &mut HashSet<NodeId>,
        rec_stack: &mut HashSet<NodeId>,
    ) -> bool {
        visited.insert(node);
        rec_stack.insert(node);

        for edge in self.graph.edges_directed(node, Direction::Outgoing) {
            let child = edge.target();
            if !visited.contains(&child) {
                if self.has_cycle(child, visited, rec_stack) {
                    return true;
                }
            } else if rec_stack.contains(&child) {
                return true;
            }
        }

        rec_stack.remove(&node);
        false
    }

    /// Find the path(s) from root to a leaf node.
    ///
    /// For admixed nodes, there are multiple paths weighted by admixture proportions.
    /// Returns a list of (path, weight) tuples where weight is the product of
    /// admixture proportions along the path.
    pub fn paths_to_leaf(&self, leaf: NodeId) -> Vec<(Vec<NodeId>, f64)> {
        let mut results = Vec::new();
        self.find_paths_up(leaf, vec![leaf], 1.0, &mut results);
        results
    }

    /// Recursive helper for finding paths (going backwards from leaf to root).
    fn find_paths_up(
        &self,
        node: NodeId,
        current_path: Vec<NodeId>,
        current_weight: f64,
        results: &mut Vec<(Vec<NodeId>, f64)>,
    ) {
        let incoming: Vec<_> = self.incoming_edges(node);

        if incoming.is_empty() {
            // Reached root
            let mut path = current_path;
            path.reverse();
            results.push((path, current_weight));
            return;
        }

        for (parent, edge_data) in incoming {
            let mut new_path = current_path.clone();
            new_path.push(parent);

            let new_weight = if edge_data.edge_type == EdgeType::Admixture {
                current_weight * edge_data.weight
            } else {
                current_weight
            };

            self.find_paths_up(parent, new_path, new_weight, results);
        }
    }

    /// Compute the expected f2 between two populations based on graph structure.
    ///
    /// f2(A,B) = sum of drift along unique path segments
    pub fn expected_f2(&self, pop_a: &str, pop_b: &str) -> Result<f64> {
        let node_a = self
            .get_node(pop_a)
            .ok_or_else(|| PopmixError::PopulationNotFound(pop_a.to_string()))?;
        let node_b = self
            .get_node(pop_b)
            .ok_or_else(|| PopmixError::PopulationNotFound(pop_b.to_string()))?;

        // Get all paths for both populations
        let paths_a = self.paths_to_leaf(node_a);
        let paths_b = self.paths_to_leaf(node_b);

        // For each pair of paths, compute the drift on non-shared segments
        let mut total_f2 = 0.0;

        for (path_a, weight_a) in &paths_a {
            for (path_b, weight_b) in &paths_b {
                let combined_weight = weight_a * weight_b;
                let path_f2 = self.compute_path_drift(path_a, path_b);
                total_f2 += combined_weight * path_f2;
            }
        }

        Ok(total_f2)
    }

    /// Compute drift along non-shared path segments.
    ///
    /// For f2(A,B), we need to sum drift on the path segments that are unique
    /// to each population (i.e., from LCA to A and from LCA to B).
    fn compute_path_drift(&self, path_a: &[NodeId], path_b: &[NodeId]) -> f64 {
        // Paths go from root to leaf: [root, ..., internal, ..., leaf]
        // Find the LCA as the last common node in both paths
        let path_b_set: HashSet<_> = path_b.iter().copied().collect();

        // Find LCA: the deepest node that appears in both paths
        // Since paths go root→leaf, LCA is the last common node
        let mut lca_idx_a = 0;
        for (i, node) in path_a.iter().enumerate() {
            if path_b_set.contains(node) {
                lca_idx_a = i;
            }
        }

        let lca = path_a[lca_idx_a];
        let lca_idx_b = path_b.iter().position(|&n| n == lca).unwrap_or(0);

        // Sum drift on path from LCA to A (exclusive of LCA)
        let mut drift = 0.0;
        for i in lca_idx_a..path_a.len().saturating_sub(1) {
            let from = path_a[i];
            let to = path_a[i + 1];
            drift += self.get_edge_drift(from, to);
        }

        // Sum drift on path from LCA to B (exclusive of LCA)
        for i in lca_idx_b..path_b.len().saturating_sub(1) {
            let from = path_b[i];
            let to = path_b[i + 1];
            drift += self.get_edge_drift(from, to);
        }

        drift
    }

    /// Get drift weight for an edge between two adjacent nodes.
    fn get_edge_drift(&self, from: NodeId, to: NodeId) -> f64 {
        for edge in self.graph.edges_directed(from, Direction::Outgoing) {
            if edge.target() == to && edge.weight().edge_type == EdgeType::Drift {
                return edge.weight().weight;
            }
        }
        0.0
    }

    /// Get all edge weights as a vector (for optimization).
    pub fn get_edge_weights(&self) -> Vec<f64> {
        self.graph.edge_weights().map(|e| e.weight).collect()
    }

    /// Set edge weights from a vector (for optimization).
    pub fn set_edge_weights(&mut self, weights: &[f64]) {
        for (idx, edge_idx) in self.graph.edge_indices().enumerate() {
            if idx < weights.len() {
                if let Some(edge) = self.graph.edge_weight_mut(edge_idx) {
                    edge.weight = weights[idx];
                }
            }
        }
    }

    /// Get number of drift edges (parameters to optimize).
    pub fn num_drift_edges(&self) -> usize {
        self.graph
            .edge_weights()
            .filter(|e| e.edge_type == EdgeType::Drift)
            .count()
    }
}

impl Default for AdmixGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// Builder for creating admixture graphs more conveniently.
pub struct GraphBuilder {
    graph: AdmixGraph,
}

impl GraphBuilder {
    /// Create a new graph builder.
    pub fn new() -> Self {
        Self {
            graph: AdmixGraph::new(),
        }
    }

    /// Add the root node.
    pub fn root(mut self, name: &str) -> Self {
        self.graph.add_node(AdmixNode::root(name));
        self
    }

    /// Add a leaf node with drift edge from parent.
    pub fn leaf(mut self, name: &str, parent: &str, drift: f64) -> Self {
        let leaf_id = self.graph.add_node(AdmixNode::leaf(name));
        if let Some(parent_id) = self.graph.get_node(parent) {
            self.graph.add_drift_edge(parent_id, leaf_id, drift);
        }
        self
    }

    /// Add an internal node with drift edge from parent.
    pub fn internal(mut self, name: &str, parent: &str, drift: f64) -> Self {
        let node_id = self.graph.add_node(AdmixNode::internal(name));
        if let Some(parent_id) = self.graph.get_node(parent) {
            self.graph.add_drift_edge(parent_id, node_id, drift);
        }
        self
    }

    /// Add an admixture node with two source edges.
    pub fn admixed(
        mut self,
        name: &str,
        source1: &str,
        prop1: f64,
        source2: &str,
        prop2: f64,
    ) -> Self {
        let node_id = self.graph.add_node(AdmixNode::admixed(name));
        if let Some(src1_id) = self.graph.get_node(source1) {
            self.graph.add_admix_edge(src1_id, node_id, prop1);
        }
        if let Some(src2_id) = self.graph.get_node(source2) {
            self.graph.add_admix_edge(src2_id, node_id, prop2);
        }
        self
    }

    /// Build the graph.
    pub fn build(self) -> Result<AdmixGraph> {
        self.graph.validate()?;
        Ok(self.graph)
    }

    /// Build without validation (for testing).
    pub fn build_unchecked(self) -> AdmixGraph {
        self.graph
    }
}

impl Default for GraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Create a simple tree: R -> (A, B)
    fn simple_tree() -> AdmixGraph {
        GraphBuilder::new()
            .root("R")
            .leaf("A", "R", 0.1)
            .leaf("B", "R", 0.2)
            .build_unchecked()
    }

    /// Create a tree with internal node: R -> I -> (A, B)
    fn tree_with_internal() -> AdmixGraph {
        GraphBuilder::new()
            .root("R")
            .internal("I", "R", 0.05)
            .leaf("A", "I", 0.1)
            .leaf("B", "I", 0.2)
            .build_unchecked()
    }

    /// Create a graph with admixture:
    /// R -> I1 -> A
    /// R -> I2 -> B
    /// I1 + I2 (50/50) -> Admix -> C
    fn admixed_graph() -> AdmixGraph {
        let mut graph = AdmixGraph::new();

        let r = graph.add_node(AdmixNode::root("R"));
        let i1 = graph.add_node(AdmixNode::internal("I1"));
        let i2 = graph.add_node(AdmixNode::internal("I2"));
        let a = graph.add_node(AdmixNode::leaf("A"));
        let b = graph.add_node(AdmixNode::leaf("B"));
        let admix = graph.add_node(AdmixNode::admixed("Admix"));
        let c = graph.add_node(AdmixNode::leaf("C"));

        graph.add_drift_edge(r, i1, 0.05);
        graph.add_drift_edge(r, i2, 0.05);
        graph.add_drift_edge(i1, a, 0.1);
        graph.add_drift_edge(i2, b, 0.1);
        graph.add_admix_edge(i1, admix, 0.5);
        graph.add_admix_edge(i2, admix, 0.5);
        graph.add_drift_edge(admix, c, 0.02);

        graph
    }

    #[test]
    fn test_simple_tree_structure() {
        let graph = simple_tree();

        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 2);
        assert!(graph.root().is_some());
        assert_eq!(graph.leaves().len(), 2);
    }

    #[test]
    fn test_tree_with_internal_structure() {
        let graph = tree_with_internal();

        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.edge_count(), 3);
        assert_eq!(graph.leaves().len(), 2);
        assert_eq!(graph.leaf_names(), vec!["A", "B"]);
    }

    #[test]
    fn test_admixed_graph_structure() {
        let graph = admixed_graph();

        assert_eq!(graph.node_count(), 7);
        assert!(graph.is_admixed(graph.get_node("Admix").unwrap()));
        assert_eq!(graph.leaves().len(), 3);
    }

    #[test]
    fn test_graph_validation_valid() {
        let graph = simple_tree();
        assert!(graph.validate().is_ok());

        let graph = admixed_graph();
        assert!(graph.validate().is_ok());
    }

    #[test]
    fn test_paths_to_leaf_simple() {
        let graph = simple_tree();

        let leaf_a = graph.get_node("A").unwrap();
        let paths = graph.paths_to_leaf(leaf_a);

        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].1, 1.0); // Weight is 1.0 for non-admixed path
    }

    #[test]
    fn test_paths_to_leaf_admixed() {
        let graph = admixed_graph();

        let leaf_c = graph.get_node("C").unwrap();
        let paths = graph.paths_to_leaf(leaf_c);

        // C has two paths through the admixture node
        assert_eq!(paths.len(), 2);

        // Both paths should have weight 0.5
        for (_, weight) in &paths {
            assert!((weight - 0.5).abs() < 1e-10);
        }
    }

    #[test]
    fn test_expected_f2_simple_tree() {
        let graph = simple_tree();

        // f2(A, B) = drift(R->A) + drift(R->B) = 0.1 + 0.2 = 0.3
        let f2 = graph.expected_f2("A", "B").unwrap();
        assert!((f2 - 0.3).abs() < 1e-10, "Expected 0.3, got {}", f2);
    }

    #[test]
    fn test_expected_f2_self() {
        let graph = simple_tree();

        // f2(A, A) = 0
        let f2 = graph.expected_f2("A", "A").unwrap();
        assert!(f2.abs() < 1e-10, "Expected 0, got {}", f2);
    }

    #[test]
    fn test_expected_f2_tree_with_internal() {
        let graph = tree_with_internal();

        // f2(A, B) = drift(I->A) + drift(I->B) = 0.1 + 0.2 = 0.3
        // (shared drift R->I doesn't contribute)
        let f2 = graph.expected_f2("A", "B").unwrap();
        assert!((f2 - 0.3).abs() < 1e-10, "Expected 0.3, got {}", f2);
    }

    #[test]
    fn test_get_set_edge_weights() {
        let mut graph = simple_tree();

        let weights = graph.get_edge_weights();
        assert_eq!(weights.len(), 2);

        // Modify weights
        let new_weights = vec![0.5, 0.6];
        graph.set_edge_weights(&new_weights);

        let updated = graph.get_edge_weights();
        assert_eq!(updated, new_weights);
    }

    #[test]
    fn test_num_drift_edges() {
        let graph = simple_tree();
        assert_eq!(graph.num_drift_edges(), 2);

        let graph = admixed_graph();
        assert_eq!(graph.num_drift_edges(), 5); // 5 drift edges, 2 admixture edges
    }

    #[test]
    fn test_builder_pattern() {
        let graph = GraphBuilder::new()
            .root("Root")
            .internal("Anc", "Root", 0.1)
            .leaf("Pop1", "Anc", 0.05)
            .leaf("Pop2", "Anc", 0.08)
            .build();

        assert!(graph.is_ok());
        let graph = graph.unwrap();
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.leaves().len(), 2);
    }
}
