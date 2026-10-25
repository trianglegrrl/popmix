#!/usr/bin/env Rscript
#
# Reference for popmix's DOT output: admixtools2 write_dot() of the dg_6pop
# qpgraph fits in qpgraph.json (with their fitted edges) and of a bare
# graph, and parse_dot() of each.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_dot.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/dot.json.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
  library(dplyr)
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
graphs <- fromJSON(file.path(dest, "qpgraph.json"), simplifyVector = FALSE)$graphs

dot_text <- function(x, ...) {
  f <- tempfile(fileext = ".dot")
  write_dot(x, outfile = f, ...)
  txt <- paste(readLines(f), collapse = "\n")
  list(text = txt, parsed = lapply(seq_len(nrow(as_edgelist(parse_dot(f)))), function(i) as_edgelist(parse_dot(f))[i, ]))
}

out <- lapply(graphs, function(g) {
  fe <- g$fit_edges
  fit <- tibble(from = vapply(fe, `[[`, "", "from"), to = vapply(fe, `[[`, "", "to"),
                type = vapply(fe, `[[`, "", "type"), weight = vapply(fe, `[[`, 0, "weight"))
  el <- do.call(rbind, lapply(g$edges, unlist))
  list(name = g$name, fit_edges = g$fit_edges,
       fitted = dot_text(fit, title = g$name),
       bare = dot_text(graph_from_edgelist(el), size1 = 5, size2 = 4.25))
})
writeLines(toJSON(list(admixtools_version = as.character(packageVersion("admixtools")), graphs = out),
                  auto_unbox = TRUE, digits = NA, pretty = FALSE), file.path(dest, "dot.json"))
cat("wrote", file.path(dest, "dot.json"), "\n")
