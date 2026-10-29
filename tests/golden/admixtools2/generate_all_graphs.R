#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
# Calls admixtools2's internals (admixtools:::), so GPL-3 like admixtools2
# (COPYING here); test tooling, in no popmix package (PROVENANCE.md).
#
# Generate the admixtools2 reference for popmix's exhaustive graph search.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_all_graphs.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/all_graphs.json. For each
# (leaves, nadmix) case it holds every graph generate_all_graphs() returns, in
# its order:
#   - key: what graph_hash() digests (the sorted rows of shortest-path
#     distances from each internal node to each leaf, then the sorted leaves),
#     with the rows sorted in C order so the key is portable; two graphs are
#     the same graph in admixtools2 exactly when their keys are equal
#   - edges: the graph's edge list, reordered so the path from the root to the
#     first leaf comes first; qpgraph() takes its leaves in vertex order and
#     fits f3 with the first one as base, so this fixes the base
#   - score: qpgraph() score of that edge list (seed 1), and the spread of the
#     scores over its random starts
#   - weights: its best admixture weights, one per admixture node in vertex
#     order: the weight of the node's incoming edge whose source vertex comes
#     first (popmix's admix_weights), so the fit at those weights can be
#     checked without an optimizer
#
# The f2 blocks are those of the dg_6pop config in generate_reference.R.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
})

# As in generate_reference.R: admixtools2 2.0.10 qpgraph() with igraph >= 2.
local({
  imports <- parent.env(asNamespace("admixtools"))
  shim <- function(graph, vp, ...) igraph::get_edge_ids(graph, as.vector(vp), ...)
  unlockBinding("get.edge.ids", imports)
  assign("get.edge.ids", shim, envir = imports)
  lockBinding("get.edge.ids", imports)
})

# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache;
# its manifest.json records the release files' checksums.
AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
AADR_PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
AADR_DATASET <- jsonlite::read_json(file.path(AADR_DIR, "manifest.json"))
POPS <- c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF")
CASES <- list(
  list(leaves = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 0),
  list(leaves = c("Chimp.REF", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 1),
  list(leaves = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 1),
  list(leaves = c("Chimp.REF", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 2)
)

graph_key <- function(graph) {
  leaves <- sort(admixtools:::get_leafnames(graph), method = "radix")
  internal <- setdiff(names(V(graph)), leaves)
  rows <- apply(distances(graph, internal, leaves, mode = "out"), 1,
                function(x) paste(sub("Inf", ".", as.character(x)), collapse = ""))
  c(sort(unname(rows), method = "radix"), leaves)
}

base_first <- function(graph, base) {
  el <- as_edgelist(graph)
  root <- names(which(degree(graph, mode = "in") == 0))
  p <- names(shortest_paths(graph, root, base, mode = "out")$vpath[[1]])
  on_path <- match(paste(head(p, -1), p[-1]), paste(el[, 1], el[, 2]))
  el[c(on_path, setdiff(seq_len(nrow(el)), on_path)), , drop = FALSE]
}

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
dest <- file.path(script_dir, paste0("v", version), "dg_6pop")

f2_blocks <- f2_from_geno(AADR_PREFIX, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                          verbose = FALSE)

cases <- lapply(CASES, function(case) {
  graphs <- generate_all_graphs(case$leaves, case$nadmix, verbose = FALSE)
  keys <- lapply(graphs, graph_key)
  stopifnot(!anyDuplicated(vapply(keys, paste, "", collapse = "\n")))
  t0 <- Sys.time()
  rows <- lapply(seq_along(graphs), function(i) {
    el <- base_first(graphs[[i]], case$leaves[1])
    fit <- qpgraph(f2_blocks, el, seed = 1)
    stopifnot(fit$edges$from == el[, 1], fit$edges$to == el[, 2])
    opt_scores <- if (is.null(fit$opt)) fit$score else fit$opt$value
    g <- graph_from_edgelist(el)
    admixnodes <- which(degree(g, mode = "in") == 2)
    # row of el for that edge: the incoming edge with the earliest source vertex
    first_in <- vapply(admixnodes, function(v) {
      rows <- which(el[, 2] == names(V(g))[v])
      rows[which.min(match(el[rows, 1], names(V(g))))]
    }, 1L)
    list(key = keys[[i]], edges = el, score = fit$score,
         score_spread = max(opt_scores) - min(opt_scores), weights = I(fit$edges$weight[first_in]))
  })
  cat(sprintf("%d leaves, %d admixture: %d graphs, %.0f s\n", length(case$leaves), case$nadmix,
              length(graphs), as.numeric(Sys.time() - t0, units = "secs")))
  list(leaves = case$leaves, nadmix = case$nadmix, graphs = rows)
})

write_json(list(admixtools_version = version, dataset = AADR_DATASET, pops = POPS, cases = cases),
           file.path(dest, "all_graphs.json"), digits = I(17), auto_unbox = TRUE, matrix = "rowmajor")
