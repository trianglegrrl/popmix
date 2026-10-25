//! qpGraph admixture graph fitting module.
//!
//! This module implements admixture graph construction and fitting:
//! - AdmixGraph: Directed acyclic graph with drift and admixture edges
//! - Path computation for expected f-statistics
//! - Edge length optimization to fit observed data
//! - qpgraph: admixtools2's graph fit
//! - Graph search over topologies, scored by qpgraph

mod fitting;
mod graph;
mod graph_search;
mod qpgraph_fit;

pub use fitting::{fit_graph, worst_residual, FitOptions, GraphFitResult};
pub use graph::{AdmixGraph, AdmixNode, EdgeType, GraphBuilder, NodeId};
pub use graph_search::{
    all_graphs, find_graphs, graph_key, GraphSearchMethod, GraphSearchOptions, GraphSearchResult,
    ScoredGraph,
};
pub use qpgraph_fit::{
    qpgraph, QpGraphEdge, QpGraphF3, QpGraphOptions, QpGraphResult, QpGraphStart,
};
