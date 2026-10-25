#!/usr/bin/env Rscript
#
# Reference for popmix's graph constraints: summarize_numadmix,
# summarize_zerof4 and satisfies_oneevent (behind satisfies_eventorder).
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_graph_constraints.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/graph_constraints.json:
# for every 4th graph of each case in all_graphs.json, the number of
# admixture events on each population's paths, the zero f4-statistics
# (summarize_zerof4 draws random edge lengths: set.seed(1) first; the set is
# the same for almost every draw), and the status (TRUE, FALSE or NA) of
# each event-order constraint, of both types.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
ref <- fromJSON(file.path(dest, "all_graphs.json"), simplifyVector = FALSE)

events <- function(leaves) {
  # (earlier1, earlier2, later1, later2); NA: the parent of the first
  ev <- list(c("Chimp.REF", "Yoruba.DG", "French.DG", "Han.DG"),
             c("Yoruba.DG", "French.DG", "French.DG", "Han.DG"),
             c("French.DG", "Han.DG", "Han.DG", NA),
             c("Yoruba.DG", NA, "French.DG", "Han.DG"),
             c("Chimp.REF", NA, "Yoruba.DG", NA),
             c("Han.DG", "French.DG", "Chimp.REF", "Yoruba.DG"))
  if ("Mbuti.DG" %in% leaves) ev <- c(ev, list(c("Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG"),
                                               c("Mbuti.DG", NA, "Han.DG", "French.DG")))
  ev
}

cases <- lapply(ref$cases, function(case) {
  idx <- seq(1, length(case$graphs), by = 4)
  graphs <- lapply(idx, function(i) {
    el <- do.call(rbind, lapply(case$graphs[[i]]$edges, unlist))
    g <- graph_from_edgelist(el)
    na <- summarize_numadmix(g)
    set.seed(1)
    z <- summarize_zerof4(g)
    ev <- events(case$leaves)
    status <- lapply(ev, function(e) lapply(1:2, function(type) {
      s <- admixtools:::satisfies_oneevent(g, e[1], e[2], e[3], e[4], type = type)
      if (is.na(s)) NULL else s
    }))
    list(index = i - 1, edges = case$graphs[[i]]$edges,
         numadmix = setNames(as.list(na$nadmix), na$pop),
         zerof4 = lapply(seq_len(nrow(z)), function(k) unlist(z[k, ])),
         events = lapply(seq_along(ev), function(k) list(event = ev[[k]], type1 = status[[k]][[1]],
                                                         type2 = status[[k]][[2]])))
  })
  cat(length(case$leaves), "leaves,", case$nadmix, "admixture:", length(graphs), "graphs\n")
  list(leaves = case$leaves, nadmix = case$nadmix, graphs = graphs)
})

writeLines(toJSON(list(admixtools_version = as.character(packageVersion("admixtools")), cases = cases),
                  auto_unbox = TRUE, null = "null", na = "null", pretty = FALSE),
           file.path(dest, "graph_constraints.json"))
cat("wrote", file.path(dest, "graph_constraints.json"), "\n")
