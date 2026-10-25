#!/usr/bin/env Rscript
#
# Reference for popmix's resampling functions: boo_list,
# qpgraph_resample_snps, qpgraph_resample_multi, compare_fits and the
# qpadm / qpwave / qpdstat / qp3pop _resample_snps functions.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_resample.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/resample.json, on the
# dg_6pop f2 blocks of generate_all_graphs.R. Each bootstrap starts from
# set.seed(<seed>) just before the call, as recorded; qpgraph gets seed = 1,
# so its random starts are reproducible inside the resampling (which
# otherwise draws them from furrr's parallel streams). Numbers are "%.17g"
# strings, so they round-trip exactly.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(igraph)
  library(dplyr)
  library(tidyr)
})
future::plan("sequential")

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
f2_blocks <- f2_from_geno(AADR_PREFIX, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                          verbose = FALSE)

num17 <- function(x) sprintf("%.17g", x)
el <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
GRAPHS <- list(
  tree = el("R", "Chimp.REF", "R", "n1", "n1", "Mbuti.DG", "n1", "n2", "n2", "Yoruba.DG", "n2", "n3",
            "n3", "French.DG", "n3", "Han.DG"),
  admix = el("R", "Chimp.REF", "R", "n1", "n1", "Mbuti.DG", "n1", "n2", "n2", "Yoruba.DG", "n2", "n3",
             "n3", "eW", "n3", "Han.DG", "eW", "adm", "n2", "adm", "adm", "French.DG")
)
admix_weight <- function(edges, field = "weight") {
  # the weight of the admixture node's incoming edge from the vertex listed
  # first, or its low/high over the random starts
  e <- edges %>% filter(type == "admix")
  if (nrow(e) == 0) return(NULL)
  num17(e[[field]][1])
}
qpgraph_rows <- function(res) lapply(seq_len(nrow(res)), function(i) list(
  id = res$id[i], sel = if ("sel" %in% names(res)) I(res$sel[[i]]) else NULL,
  score = num17(res$score[i]),
  score_test = if ("score_test" %in% names(res)) num17(res$score_test[i]) else NULL,
  weight = admix_weight(res$edges[[i]]),
  low = admix_weight(res$edges[[i]], "low"), high = admix_weight(res$edges[[i]], "high"),
  error = !is.null(res$error[[i]])))

out <- list(admixtools_version = as.character(packageVersion("admixtools")), dataset = AADR_DATASET, pops = POPS,
            graphs = lapply(GRAPHS, function(g) lapply(seq_len(nrow(g)), function(i) g[i, ])))

set.seed(1)
boo <- boo_list(f2_blocks, nboot = 4)
out$boo_list <- list(seed = 1, nboot = 4, sel = lapply(boo$sel, I))

set.seed(2)
fits <- qpgraph_resample_multi(f2_blocks, GRAPHS, nboot = 6, seed = 1)
out$qpgraph_resample_multi <- list(seed = 2, nboot = 6, qpgraph_seed = 1,
                                   fits = setNames(lapply(fits, qpgraph_rows), names(GRAPHS)))
cf <- compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
out$compare_fits <- lapply(cf, num17)
# ties and a zero difference, which compare_fits drops
cf2 <- compare_fits(c(1, 2, 3, 4, 5, 6), c(0.5, 2, 4, 3.5, 5.5, 1))
out$compare_fits_small <- list(scores1 = c(1, 2, 3, 4, 5, 6), scores2 = c(0.5, 2, 4, 3.5, 5.5, 1),
                               result = lapply(cf2, num17))

out$qpgraph_resample_snps_jack <- qpgraph_rows(qpgraph_resample_snps(f2_blocks, graph = GRAPHS$admix, seed = 1))
set.seed(3)
out$qpgraph_resample_snps_boot <- list(seed = 3, boot = 5, rows = qpgraph_rows(
  qpgraph_resample_snps(f2_blocks, boot = 5, graph = GRAPHS$admix, seed = 1)))

QPADM <- list(target = "French.DG", left = c("Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Papuan.DG", "Chimp.REF"))
qa <- qpadm_resample_snps(f2_blocks, target = QPADM$target, left = QPADM$left, right = QPADM$right)
out$qpadm_resample_snps <- c(QPADM, list(rows = lapply(seq_len(nrow(qa)), function(i) list(
  id = qa$id[i], weights = num17(qa$weights[[i]]$weight),
  chisq = num17(qa$rankdrop[[i]]$chisq[1]), p = num17(qa$rankdrop[[i]]$p[1])))))

QPWAVE <- list(left = c("French.DG", "Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Papuan.DG", "Chimp.REF"))
set.seed(4)
qw <- qpwave_resample_snps(f2_blocks, boot = 5, left = QPWAVE$left, right = QPWAVE$right)
out$qpwave_resample_snps <- c(QPWAVE, list(seed = 4, boot = 5, rows = lapply(seq_len(nrow(qw)), function(i) list(
  id = qw$id[i], sel = I(qw$sel[[i]]), f4rank = qw$rankdrop[[i]]$f4rank,
  chisq = num17(qw$rankdrop[[i]]$chisq), p = num17(qw$rankdrop[[i]]$p)))))

qd <- qpdstat_resample_snps(f2_blocks, pop1 = "Chimp.REF", pop2 = "Mbuti.DG", pop3 = "French.DG", pop4 = "Han.DG")
out$qpdstat_resample_snps <- list(pops = c("Chimp.REF", "Mbuti.DG", "French.DG", "Han.DG"),
                                  rows = lapply(seq_len(nrow(qd)), function(i) list(
                                    id = qd$id[i], est = num17(qd$est[i]), se = num17(qd$se[i]))))
q3 <- qp3pop_resample_snps(f2_blocks, pop1 = "French.DG", pop2 = "Han.DG", pop3 = "Yoruba.DG")
out$qp3pop_resample_snps <- list(pops = c("French.DG", "Han.DG", "Yoruba.DG"),
                                 rows = lapply(seq_len(nrow(q3)), function(i) list(
                                   id = q3$id[i], est = num17(q3$est[i]), se = num17(q3$se[i]))))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
writeLines(toJSON(out, auto_unbox = TRUE, null = "null", pretty = FALSE), file.path(dest, "resample.json"))
cat("wrote", file.path(dest, "resample.json"), "\n")
