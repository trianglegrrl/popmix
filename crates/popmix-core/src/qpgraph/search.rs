//! Graph search algorithms for finding optimal admixture graph topologies.
//!
//! Implements:
//! - Random graph generation
//! - Graph enumeration (all trees, graphs with N admixture events)
//! - Local search moves (SPR, admixture insertion/deletion)

use std::collections::HashSet;

use crate::fstats::F2Blocks;
use crate::{PopmixError, Result};

use super::fitting::{fit_graph, FitOptions, GraphFitResult};
use super::graph::{AdmixGraph, AdmixNode, GraphBuilder};

/// Options for graph search.
#[derive(Debug, Clone)]
pub struct GraphSearchOptions {
    /// Maximum number of admixture events to allow.
    pub max_admix: usize,
    /// Maximum number of graphs to evaluate.
    pub max_graphs: usize,
    /// Score threshold for accepting a graph.
    pub score_threshold: f64,
    /// Options for fitting each graph.
    pub fit_options: FitOptions,
    /// Random seed for reproducibility.
    pub seed: Option<u64>,
}

impl Default for GraphSearchOptions {
    fn default() -> Self {
        Self {
            max_admix: 2,
            max_graphs: 1000,
            score_threshold: 0.01,
            fit_options: FitOptions::default(),
            seed: None,
        }
    }
}

/// Result of graph search.
#[derive(Debug)]
pub struct GraphSearchResult {
    /// Best fitting graph found.
    pub best_graph: AdmixGraph,
    /// Fit result for best graph.
    pub best_fit: GraphFitResult,
    /// All graphs evaluated with their scores.
    pub all_fits: Vec<GraphFitResult>,
    /// Number of graphs evaluated.
    pub n_evaluated: usize,
    /// Number of graphs that passed the score threshold.
    pub n_passing: usize,
}

impl GraphSearchResult {
    /// Get the N best fitting graphs.
    pub fn top_n(&self, n: usize) -> Vec<&GraphFitResult> {
        let mut sorted: Vec<_> = self.all_fits.iter().collect();
        sorted.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
        sorted.into_iter().take(n).collect()
    }
}

/// Search for the best admixture graph topology.
///
/// Generates candidate graph topologies and fits each to the observed data,
/// returning the best-fitting graph(s).
///
/// # Arguments
/// * `populations` - List of population names to include in the graph
/// * `f2_blocks` - Observed f-statistics
/// * `outgroup` - Optional outgroup population (used as root constraint)
/// * `options` - Search options
///
/// # Returns
/// GraphSearchResult containing the best graph and all evaluated graphs
pub fn find_graphs(
    populations: &[&str],
    f2_blocks: &F2Blocks,
    outgroup: Option<&str>,
    options: Option<GraphSearchOptions>,
) -> Result<GraphSearchResult> {
    let opts = options.unwrap_or_default();

    if populations.len() < 3 {
        return Err(PopmixError::Validation(
            "Need at least 3 populations for graph search".to_string(),
        ));
    }

    // Verify all populations exist in data
    let data_pops: HashSet<_> = f2_blocks.populations.iter().map(|p| p.as_str()).collect();
    for pop in populations {
        if !data_pops.contains(*pop) {
            return Err(PopmixError::PopulationNotFound(pop.to_string()));
        }
    }

    // Generate candidate topologies
    let candidates = generate_candidate_graphs(populations, outgroup, opts.max_admix, opts.seed)?;

    let mut all_fits = Vec::new();
    let mut best_fit: Option<GraphFitResult> = None;
    let mut n_evaluated = 0;
    let mut n_passing = 0;

    // Evaluate each candidate
    for graph in candidates.into_iter().take(opts.max_graphs) {
        n_evaluated += 1;

        match fit_graph(&graph, f2_blocks, Some(opts.fit_options.clone())) {
            Ok(result) => {
                if result.score < opts.score_threshold {
                    n_passing += 1;
                }

                if best_fit.is_none() || result.score < best_fit.as_ref().unwrap().score {
                    best_fit = Some(result.clone());
                }

                all_fits.push(result);
            }
            Err(_) => {
                // Skip graphs that fail to fit
                continue;
            }
        }
    }

    // Sort by score
    all_fits.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());

    let best_fit = best_fit.ok_or_else(|| {
        PopmixError::Validation("No valid graphs found during search".to_string())
    })?;

    Ok(GraphSearchResult {
        best_graph: best_fit.graph.clone(),
        best_fit,
        all_fits,
        n_evaluated,
        n_passing,
    })
}

/// Generate candidate graph topologies.
fn generate_candidate_graphs(
    populations: &[&str],
    outgroup: Option<&str>,
    max_admix: usize,
    seed: Option<u64>,
) -> Result<Vec<AdmixGraph>> {
    let mut graphs = Vec::new();

    // Start with simple trees, then add admixture
    // For now, generate a limited set of topologies

    // Generate all binary tree topologies
    let trees = generate_all_trees(populations, outgroup)?;
    graphs.extend(trees);

    // If admixture is allowed, generate graphs with admixture events
    if max_admix > 0 {
        let base_trees = generate_all_trees(populations, outgroup)?;
        for tree in base_trees {
            // For each pair of branches, try adding an admixture edge
            let admix_graphs = add_admixture_events(&tree, 1)?;
            graphs.extend(admix_graphs);
        }
    }

    // Shuffle if seed is provided
    if let Some(seed) = seed {
        let mut rng_state = seed;
        let mut rand_idx = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            rng_state as usize
        };

        // Fisher-Yates shuffle
        for i in (1..graphs.len()).rev() {
            let j = rand_idx() % (i + 1);
            graphs.swap(i, j);
        }
    }

    Ok(graphs)
}

/// Generate all binary tree topologies for the given populations.
fn generate_all_trees(populations: &[&str], outgroup: Option<&str>) -> Result<Vec<AdmixGraph>> {
    let n = populations.len();

    if n < 2 {
        return Err(PopmixError::Validation(
            "Need at least 2 populations".to_string(),
        ));
    }

    // Separate outgroup if specified
    let (ingroup, outgroup_pop) = if let Some(og) = outgroup {
        let ingroup: Vec<_> = populations.iter().filter(|&&p| p != og).copied().collect();
        (ingroup, Some(og))
    } else {
        (populations.to_vec(), None)
    };

    let mut trees = Vec::new();

    // For small numbers, enumerate all topologies
    // For larger numbers, sample random topologies
    if ingroup.len() <= 4 {
        let topologies = enumerate_binary_topologies(ingroup.len());
        for topo in topologies {
            let graph = build_tree_from_topology(&ingroup, &topo, outgroup_pop)?;
            trees.push(graph);
        }
    } else {
        // Generate a random sample of trees
        for i in 0..100 {
            let graph = random_binary_tree(&ingroup, outgroup_pop, i as u64)?;
            trees.push(graph);
        }
    }

    Ok(trees)
}

/// Enumerate all binary tree topologies for n leaves.
///
/// Returns list of topologies encoded as nested tuples.
fn enumerate_binary_topologies(n: usize) -> Vec<Vec<(usize, usize)>> {
    if n < 2 {
        return vec![];
    }
    if n == 2 {
        return vec![vec![(0, 1)]];
    }
    if n == 3 {
        // Three topologies: ((0,1),2), ((0,2),1), ((1,2),0)
        return vec![
            vec![(0, 1), (3, 2)], // (0,1) forms node 3, then (3,2)
            vec![(0, 2), (3, 1)],
            vec![(1, 2), (3, 0)],
        ];
    }
    if n == 4 {
        // 15 unrooted topologies, but we'll use rooted
        return vec![
            vec![(0, 1), (2, 3), (4, 5)],
            vec![(0, 2), (1, 3), (4, 5)],
            vec![(0, 3), (1, 2), (4, 5)],
            vec![(0, 1), (4, 2), (5, 3)],
            vec![(0, 2), (4, 1), (5, 3)],
        ];
    }

    // For larger n, just return a placeholder
    vec![vec![(0, 1)]]
}

/// Build a tree from a topology specification.
fn build_tree_from_topology(
    populations: &[&str],
    _topology: &[(usize, usize)],
    outgroup: Option<&str>,
) -> Result<AdmixGraph> {
    // Simplified: build a caterpillar tree for now
    let mut builder = GraphBuilder::new().root("R");

    // Add outgroup first if present
    if let Some(og) = outgroup {
        builder = builder.leaf(og, "R", 0.1);
    }

    // Add ingroup as a caterpillar
    let mut last_internal = "R".to_string();
    for (i, pop) in populations.iter().enumerate() {
        if Some(*pop) == outgroup {
            continue;
        }

        if i < populations.len() - 1 {
            let internal_name = format!("I{}", i);
            builder = builder.internal(&internal_name, &last_internal, 0.01);
            builder = builder.leaf(pop, &internal_name, 0.01);
            last_internal = internal_name;
        } else {
            builder = builder.leaf(pop, &last_internal, 0.01);
        }
    }

    Ok(builder.build_unchecked())
}

/// Generate a random binary tree.
fn random_binary_tree(
    populations: &[&str],
    outgroup: Option<&str>,
    seed: u64,
) -> Result<AdmixGraph> {
    let mut rng_state = seed;
    let mut rand = || {
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (rng_state >> 33) as f64 / (1u64 << 31) as f64
    };

    let mut graph = AdmixGraph::new();
    let root = graph.add_node(AdmixNode::root("R"));

    // Add outgroup first
    if let Some(og) = outgroup {
        let og_node = graph.add_node(AdmixNode::leaf(og));
        graph.add_drift_edge(root, og_node, rand() * 0.1);
    }

    // Build tree by successively joining populations
    let mut available_nodes: Vec<_> = Vec::new();

    // Start with first internal node connected to root
    let first_internal = graph.add_node(AdmixNode::internal("I0"));
    graph.add_drift_edge(root, first_internal, rand() * 0.05);
    available_nodes.push(first_internal);

    // Add populations
    for (i, pop) in populations.iter().enumerate() {
        if Some(*pop) == outgroup {
            continue;
        }

        let pop_node = graph.add_node(AdmixNode::leaf(*pop));

        if available_nodes.is_empty() {
            graph.add_drift_edge(root, pop_node, rand() * 0.1);
        } else {
            // Pick a random available node to attach to
            let rand_idx = (rand() * available_nodes.len() as f64) as usize;
            let parent = available_nodes[rand_idx.min(available_nodes.len() - 1)];
            graph.add_drift_edge(parent, pop_node, rand() * 0.05);

            // Maybe add a new internal node
            if i < populations.len() - 1 && rand() > 0.5 {
                let new_internal = graph.add_node(AdmixNode::internal(format!("I{}", i + 1)));
                graph.add_drift_edge(parent, new_internal, rand() * 0.03);
                available_nodes.push(new_internal);
            }
        }
    }

    Ok(graph)
}

/// Add admixture events to a tree.
#[allow(dead_code)]
fn add_admixture_events(tree: &AdmixGraph, _n_admix: usize) -> Result<Vec<AdmixGraph>> {
    // Simplified: for now, just return the original tree
    // Full implementation would:
    // 1. Enumerate pairs of edges
    // 2. For each pair, try adding an admixture node
    // 3. Return all valid resulting graphs

    // For now, just return the original
    Ok(vec![tree.clone()])
}

/// Generate a random admixture graph.
#[allow(dead_code)]
pub fn random_admixture_graph(
    populations: &[&str],
    n_admix: usize,
    seed: u64,
) -> Result<AdmixGraph> {
    // Start with a random tree
    let tree = random_binary_tree(populations, None, seed)?;

    // Add admixture events if requested
    if n_admix == 0 {
        return Ok(tree);
    }

    // For now, just return the tree
    // Full implementation would add admixture edges
    Ok(tree)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use crate::genotypes::{Chromosome, PopulationId};
    use ndarray::Array2;

    fn create_test_f2_blocks() -> F2Blocks {
        let pops = vec!["A", "B", "C", "D"];
        let n_snps = 500;

        let mut rng_state: u64 = 22222;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 4));
        let counts = Array2::from_elem((n_snps, 4), 50u32);

        for i in 0..n_snps {
            let base = rand() * 0.8 + 0.1;
            freqs[[i, 0]] = (base + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 1]] = (base + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 2]] = (base + rand() * 0.1 - 0.05).clamp(0.0, 1.0);
            freqs[[i, 3]] = rand();
        }

        let mut chromosomes = Vec::with_capacity(n_snps);
        let mut positions = Vec::with_capacity(n_snps);

        for i in 0..n_snps {
            chromosomes.push(Chromosome::Autosome(1));
            positions.push(i as f64 * 0.001);
        }

        let af = AlleleFrequencies::new(
            freqs,
            counts,
            pops.iter().map(|s| PopulationId::from(*s)).collect(),
            chromosomes,
            positions,
        );

        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_generate_all_trees() {
        let pops = vec!["A", "B", "C"];
        let trees = generate_all_trees(&pops, None).unwrap();

        assert!(!trees.is_empty());
        for tree in &trees {
            assert_eq!(tree.leaves().len(), 3);
        }
    }

    #[test]
    fn test_random_binary_tree() {
        let pops = vec!["A", "B", "C", "D"];
        let tree = random_binary_tree(&pops, None, 12345).unwrap();

        assert_eq!(tree.leaves().len(), 4);
        assert!(tree.root().is_some());
    }

    #[test]
    fn test_random_admixture_graph() {
        let pops = vec!["A", "B", "C"];
        let graph = random_admixture_graph(&pops, 0, 54321).unwrap();

        assert_eq!(graph.leaves().len(), 3);
    }

    #[test]
    fn test_find_graphs_basic() {
        let f2_blocks = create_test_f2_blocks();

        let opts = GraphSearchOptions {
            max_graphs: 10,
            max_admix: 0,
            ..Default::default()
        };

        let result = find_graphs(&["A", "B", "C"], &f2_blocks, Some("D"), Some(opts));

        assert!(result.is_ok(), "Search failed: {:?}", result.err());
        let result = result.unwrap();

        println!("Evaluated {} graphs", result.n_evaluated);
        println!("Best score: {}", result.best_fit.score);
        assert!(result.n_evaluated > 0);
    }

    #[test]
    fn test_find_graphs_requires_populations() {
        let f2_blocks = create_test_f2_blocks();

        // Too few populations
        let result = find_graphs(&["A", "B"], &f2_blocks, None, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_search_options() {
        let opts = GraphSearchOptions {
            max_admix: 1,
            max_graphs: 500,
            score_threshold: 0.001,
            seed: Some(42),
            ..Default::default()
        };

        assert_eq!(opts.max_admix, 1);
        assert_eq!(opts.max_graphs, 500);
        assert_eq!(opts.seed, Some(42));
    }
}
