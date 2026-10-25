#!/usr/bin/env Rscript
#
# Reference for popmix's evolutionary graph search: admixtools2 find_graphs()
# on the dg_6pop f2 blocks, for two cases of all_graphs.json (5 populations
# with 1 admixture, 4 with 2), Chimp.REF as outgroup, default settings
# (numgraphs 10, stop_gen2 15, numstart 1), after set.seed(seed) for seeds 1
# to 8: the best score and the number of graphs fitted.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_find_graphs.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/find_graphs.json.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
})
future::plan("sequential")

local({
  imports <- parent.env(asNamespace("admixtools"))
  shim <- function(graph, vp, ...) igraph::get_edge_ids(graph, as.vector(vp), ...)
  unlockBinding("get.edge.ids", imports)
  assign("get.edge.ids", shim, envir = imports)
  lockBinding("get.edge.ids", imports)
})

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"
POPS <- c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF")
f2_blocks <- f2_from_geno(AADR_PREFIX, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                          verbose = FALSE)
CASES <- list(
  list(leaves = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 1),
  list(leaves = c("Chimp.REF", "Yoruba.DG", "French.DG", "Han.DG"), nadmix = 2)
)

cases <- lapply(CASES, function(case) {
  runs <- lapply(1:8, function(seed) {
    set.seed(seed)
    res <- find_graphs(f2_blocks[case$leaves, case$leaves, ], numadmix = case$nadmix, outpop = case$leaves[1],
                       verbose = FALSE)
    cat(length(case$leaves), case$nadmix, "seed", seed, min(res$score, na.rm = TRUE), nrow(res), "\n")
    list(seed = seed, best = sprintf("%.17g", min(res$score, na.rm = TRUE)), n_fitted = nrow(res),
         generations = max(res$generation))
  })
  c(case, list(runs = runs))
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
writeLines(toJSON(list(admixtools_version = as.character(packageVersion("admixtools")), cases = cases),
                  auto_unbox = TRUE, pretty = FALSE),
           file.path(dest, "find_graphs.json"))
cat("wrote", file.path(dest, "find_graphs.json"), "\n")
