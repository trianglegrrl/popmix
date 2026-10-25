#!/usr/bin/env Rscript
#
# Reference for popmix's _resample_inds functions: qpAdm, qpWave, qpdstat,
# qp3pop and qpGraph with each individual left out in turn.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_resample_inds.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/resample_inds.json.
# The samples are the first three of each dg_6pop population in the .ind
# file, and Chimp.REF. For each one in a population of more than one,
# admixtools2's f2_from_geno(inds, pops) without it (the dg_6pop parameters),
# then each statistic on those f2 blocks. This is popmix's definition of
# leaving an individual out; admixtools2's own *_resample_inds read
# extract_counts() data, a different f2 estimator (DEVIATIONS.md). Numbers
# are "%.17g" strings.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
  library(dplyr)
})

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
ind <- read.table(paste0(AADR_PREFIX, ".ind"), col.names = c("id", "sex", "pop"), colClasses = "character")
ind <- ind[ind$pop %in% POPS, ]
poplist <- do.call(rbind, lapply(POPS, function(p) head(ind[ind$pop == p, c("id", "pop")], 3)))

num17 <- function(x) sprintf("%.17g", x)
el <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
GRAPH <- el("R", "Chimp.REF", "R", "n1", "n1", "Mbuti.DG", "n1", "n2", "n2", "Yoruba.DG", "n2", "n3",
            "n3", "eW", "n3", "Han.DG", "eW", "adm", "n2", "adm", "adm", "French.DG")
QPADM <- list(target = "French.DG", left = c("Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Papuan.DG", "Chimp.REF"))
QPWAVE <- list(left = c("French.DG", "Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Papuan.DG", "Chimp.REF"))
QPDSTAT <- c("Chimp.REF", "Mbuti.DG", "French.DG", "Han.DG")
QP3POP <- c("French.DG", "Han.DG", "Yoruba.DG")

counts <- table(poplist$pop)
left_out <- poplist$id[counts[poplist$pop] > 1]
rows <- lapply(left_out, function(id) {
  keep <- poplist$id != id
  f2 <- f2_from_geno(AADR_PREFIX, inds = poplist$id[keep], pops = poplist$pop[keep], blgsize = 0.05, maxmiss = 0,
                     minmaf = 0, maxmaf = 0.5, auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2",
                     apply_corr = TRUE, verbose = FALSE)
  qa <- qpadm(f2, target = QPADM$target, left = QPADM$left, right = QPADM$right, verbose = FALSE)
  qw <- qpwave(f2, left = QPWAVE$left, right = QPWAVE$right, verbose = FALSE)
  qd <- qpdstat(f2, QPDSTAT[1], QPDSTAT[2], QPDSTAT[3], QPDSTAT[4])
  q3 <- qp3pop(f2, QP3POP[1], QP3POP[2], QP3POP[3])
  qg <- qpgraph(f2, GRAPH, seed = 1)
  adm <- qg$edges %>% filter(type == "admix")
  cat(id, "\n")
  list(ind = id, pop = poplist$pop[poplist$id == id], n_blocks = dim(f2)[3],
       qpadm = list(weights = num17(qa$weights$weight), se = num17(qa$weights$se),
                    chisq = num17(qa$rankdrop$chisq[1]), p = num17(qa$rankdrop$p[1])),
       qpwave = list(f4rank = qw$rankdrop$f4rank, chisq = num17(qw$rankdrop$chisq)),
       qpdstat = list(est = num17(qd$est), se = num17(qd$se)),
       qp3pop = list(est = num17(q3$est), se = num17(q3$se)),
       qpgraph = list(score = num17(qg$score), weight = num17(adm$weight[1]),
                      low = num17(adm$low[1]), high = num17(adm$high[1])))
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
out <- list(admixtools_version = as.character(packageVersion("admixtools")), prefix = AADR_PREFIX, dataset = AADR_DATASET,
            inds = poplist$id, pops = poplist$pop,
            params = list(blgsize = 0.05, maxmiss = 0, auto_only = TRUE, poly_only = "f2"),
            graph = lapply(seq_len(nrow(GRAPH)), function(i) GRAPH[i, ]), qpadm = QPADM, qpwave = QPWAVE,
            qpdstat = QPDSTAT, qp3pop = QP3POP, rows = rows)
writeLines(toJSON(out, auto_unbox = TRUE, null = "null", pretty = FALSE), file.path(dest, "resample_inds.json"))
cat("wrote", file.path(dest, "resample_inds.json"), "\n")
