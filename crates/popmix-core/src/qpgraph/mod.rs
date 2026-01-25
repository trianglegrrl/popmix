//! qpGraph admixture graph fitting module.
//!
//! This module implements admixture graph construction and fitting:
//! - AdmixGraph: Directed acyclic graph with drift and admixture edges
//! - Path computation for expected f-statistics
//! - Edge length optimization to fit observed data
//! - Graph search algorithms for finding optimal topologies

mod graph;
mod fitting;
mod search;

pub use graph::{AdmixGraph, AdmixNode, EdgeType, GraphBuilder, NodeId};
pub use fitting::{fit_graph, GraphFitResult, FitOptions};
pub use search::{find_graphs, GraphSearchOptions, GraphSearchResult};
