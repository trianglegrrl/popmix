#!/usr/bin/env Rscript
#
# Reference for popmix's msprime script writers: admixtools2 msprime_sim()
# and msprime_genome() for three dg_6pop graphs (the qpgraph.json tree and
# one-admixture graph, and a two-admixture graph of all_graphs.json), with
# their defaults and with other options.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_msprime.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/msprime.json: each
# script's text, the output prefix replaced by {outpref}.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
qg <- fromJSON(file.path(dest, "qpgraph.json"), simplifyVector = FALSE)$graphs
ag <- fromJSON(file.path(dest, "all_graphs.json"), simplifyVector = FALSE)$cases[[4]]$graphs
el <- function(edges) do.call(rbind, lapply(edges, unlist))
graphs <- list(tree = el(qg[[1]]$edges), admix1 = el(qg[[2]]$edges), admix2 = el(ag[[7]]$edges))

outdir <- tempfile("msprime")
dir.create(outdir)
outpref <- file.path(outdir, "sim")
text <- function(path) gsub(outpref, "{outpref}", paste(readLines(path), collapse = "\n"), fixed = TRUE)

OPTIONS <- list(
  list(fun = "msprime_sim", args = list()),
  list(fun = "msprime_sim", args = list(nsnps = 200, neff = 500, ind_per_pop = 2, mutation_rate = 1e-4,
                                        time = 500, fix_leaf = TRUE, admix_default = 0.3, numcores = 4)),
  list(fun = "msprime_genome", args = list()),
  list(fun = "msprime_genome", args = list(nchr = 3, seq_length = 1e6, recomb_rate_chr = 1e-8, ind_per_pop = 3)))

out <- list()
for (name in names(graphs)) {
  for (o in OPTIONS) {
    f <- get(o$fun, asNamespace("admixtools"))
    path <- do.call(f, c(list(graphs[[name]], outpref = outpref), o$args))
    out[[length(out) + 1]] <- list(graph = name, edges = lapply(seq_len(nrow(graphs[[name]])), function(i) graphs[[name]][i, ]),
                                   fun = o$fun, args = o$args, text = text(path))
  }
}
writeLines(toJSON(list(admixtools_version = as.character(packageVersion("admixtools")), scripts = out),
                  auto_unbox = TRUE, digits = NA, pretty = FALSE), file.path(dest, "msprime.json"))
cat("wrote", file.path(dest, "msprime.json"), "\n")
