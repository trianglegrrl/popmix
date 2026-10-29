#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
# Calls admixtools2's internals (admixtools:::), so GPL-3 like admixtools2
# (COPYING here); test tooling, in no popmix package (PROVENANCE.md).
#
# Generate admixtools2 golden master reference data.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_reference.R [config...]
#
# Each config writes to tests/golden/admixtools2/v<version>/<config>/:
#
# f2.json holds every stage of admixtools2's f2_from_geno pipeline, so a
# mismatch can be traced to the stage that causes it:
#   - ploidy assigned to each sample (adjust_pseudohaploid)
#   - SNP counts read, kept after discard_from_aftable, and polymorphic
#   - block lengths (polymorphic SNPs per block)
#   - the f2 block array, as returned by f2_from_geno
#   - f2() estimate and SE for every pair
#
# f3_f4.json holds f3() and f4() for every population combination, computed
# from the same f2 blocks.
#
# qpwave.json holds, for each qpWave model, the stages of qpwave() on the same
# f2 blocks: the f4 matrix and its jackknife covariance
# (f2blocks_to_f4stats), then the rank tests.
#
# qpadm.json holds qpadm() models (weights and SEs, rankdrop, popdrop) and
# qpf4ratio() estimates on the same f2 blocks.
#
# qpgraph.json holds qpgraph() fits of fixed graphs on the same f2 blocks:
# the f3-statistics it fits and the inverse of their covariance, the fitted
# edges and score, and the multistart table (seed 1). Because the admixture
# weights come from R's optim, each graph is also evaluated at fixed weights
# (qpgraph()'s own best weights and all 0.3) through the internals qpgraph()
# uses, which gives exact, optimizer-free reference values.
#
# fst.json holds fst() for every pair. admixtools2 computes FST from the
# genotype files on its own SNP set (all SNPs kept by maxmiss, monomorphic
# included) and blocks; fst() is given all of the config's FST populations so
# the SNP set is that of f2_from_geno on them. Its blocks
# (f2_from_geno(fst = TRUE)) are recorded too. Chimp.REF is a single
# pseudo-haploid sample: its p(1-p)/(n-1) is 0/0 at every SNP, so admixtools2
# discards every FST block and refuses; the FST populations leave it out.
#
# Numbers are written with 17 significant digits, enough to round-trip a
# double exactly.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

# admixtools2 2.0.10 passes igraph vertex sequences to get.edge.ids(), which
# igraph >= 2 rejects ("Only two-column data.frames and matrices, and vectors
# are allowed"). Passing the same vertex ids as a plain vector keeps igraph
# 1.x's meaning (consecutive pairs), so qpgraph() runs unchanged.
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

# A qpgraph graph: (from, to) pairs, one edge per row.
edgelist <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)

# Per-config params override PARAMS; constrained qpadm models with three or
# more sources are run through admixtools2's internals (see qpadm_internal).
CONFIGS <- list(
  dg_6pop = list(
    prefix = AADR_PREFIX,
    pops = c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF"),
    fst_pops = c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG"),
    qpwave = list(
      list(left = c("Han.DG", "French.DG"), right = c("Yoruba.DG", "Mbuti.DG", "Chimp.REF")),
      list(left = c("French.DG", "Han.DG", "Papuan.DG"), right = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG"))
    ),
    qpadm = list(
      list(target = "French.DG", left = c("Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Chimp.REF", "Papuan.DG")),
      list(target = "Papuan.DG", left = c("Han.DG", "French.DG"), right = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG")),
      list(target = "Papuan.DG", left = c("Han.DG", "French.DG"), right = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG"),
           constrained = TRUE)
    ),
    f4ratio = list(
      c("Han.DG", "Mbuti.DG", "French.DG", "Papuan.DG", "Yoruba.DG"),
      c("Chimp.REF", "Mbuti.DG", "French.DG", "Yoruba.DG", "Han.DG")
    ),
    qpgraph = list(
      tree = edgelist("R","Chimp.REF", "R","n1", "n1","Mbuti.DG", "n1","n2", "n2","Yoruba.DG", "n2","n3", "n3","French.DG", "n3","n4", "n4","Han.DG", "n4","Papuan.DG"),
      # French as a mix of a West Eurasian and a basal Eurasian lineage
      admix1 = edgelist("R","Chimp.REF", "R","n1", "n1","Mbuti.DG", "n1","n2", "n2","Yoruba.DG", "n2","n2b", "n2b","bE", "n2b","n3", "n3","n4", "n3","eW", "n4","Han.DG", "n4","Papuan.DG", "eW","adm", "bE","adm", "adm","French.DG")
    )
  ),
  dg_11pop = list(
    prefix = AADR_PREFIX,
    pops = c("Druze.DG", "Basque.DG", "Sardinian.DG", "Yoruba.DG", "Han.DG", "Dinka.DG",
             "Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF"),
    fst_pops = c("Druze.DG", "Sardinian.DG", "Dinka.DG", "Han.DG"),
    qpwave = list(
      list(left = c("Druze.DG", "Sardinian.DG", "Dinka.DG", "Han.DG"),
           right = c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF"))
    ),
    qpadm = list(
      list(target = "Druze.DG", left = c("Sardinian.DG", "Dinka.DG", "Han.DG"),
           right = c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF")),
      # Unconstrained Han weight is negative: constrained pins it at 0 and
      # leaves the other two inside (0, 1).
      list(target = "Druze.DG", left = c("Sardinian.DG", "Dinka.DG", "Han.DG"),
           right = c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF"),
           constrained = TRUE),
      # Unconstrained weights already inside [0, 1]: constrained must agree.
      list(target = "Druze.DG", left = c("Sardinian.DG", "Dinka.DG"),
           right = c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF"),
           constrained = TRUE)
    ),
    f4ratio = list(c("Chimp.REF", "Mbuti.DG", "Druze.DG", "Dinka.DG", "Sardinian.DG")),
    qpgraph = list(
      # Karitiana: East + West Eurasian; Sardinian: West + basal Eurasian
      admix2 = edgelist("R","Chimp.REF", "R","a", "a","Mbuti.DG", "a","b", "b","Yoruba.DG", "b","c", "c","bE", "c","ooa", "ooa","east", "ooa","west", "east","Han.DG", "east","e2", "e2","Papuan.DG", "e2","kar", "west","w1", "west","w2", "w1","kar", "w1","Basque.DG", "w2","sar", "w2","Druze.DG", "bE","sar", "kar","Karitiana.DG", "sar","Sardinian.DG")
    )
  ),
  # DATES' simulated data: 54 samples, so each PACKEDANCESTRYMAP record is
  # padded to 48 bytes, not ceil(54 / 4) = 14. Paths are relative to the
  # repository root.
  family_4pop = list(
    prefix = "tests/golden/dates/family_chr20_22",
    pops = c("French", "Yoruba", "Simulation", "Simulation.1"),
    fst_pops = c("French", "Yoruba", "Simulation.1"),
    qpwave = list(
      list(left = c("Simulation", "French"), right = c("Yoruba", "Simulation.1"))
    ),
    qpadm = list(
      list(target = "Simulation", left = c("French"), right = c("Yoruba", "Simulation.1"))
    ),
    qpgraph = list(
      tree = edgelist("R","Yoruba", "R","a", "a","French", "a","b", "b","Simulation", "b","Simulation.1")
    )
  ),
  # As ancient_5pop, but keeping SNPs missing in one population (of five for
  # f2, of four for FST).
  ancient_maxmiss = list(
    prefix = AADR_PREFIX,
    params = list(maxmiss = 0.25),
    pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG", "Chimp.REF"),
    fst_pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG"),
    qpwave = list(
      list(left = c("French.DG", "Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG"),
           right = c("Mbuti.DG", "Chimp.REF"))
    ),
    qpadm = list(
      list(target = "French.DG", left = c("Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG"),
           right = c("Chimp.REF", "Mbuti.DG"))
    ),
    f4ratio = list(
      c("Chimp.REF", "Mbuti.DG", "French.DG", "Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG")
    )
  ),
  ancient_5pop = list(
    prefix = AADR_PREFIX,
    pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG", "Chimp.REF"),
    fst_pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG"),
    qpwave = list(
      list(left = c("French.DG", "Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG"),
           right = c("Mbuti.DG", "Chimp.REF"))
    ),
    qpadm = list(
      list(target = "French.DG", left = c("Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG"),
           right = c("Chimp.REF", "Mbuti.DG"))
    ),
    f4ratio = list(
      c("Chimp.REF", "Mbuti.DG", "French.DG", "Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG")
    ),
    qpgraph = list(
      # French as a mix of Yamnaya- and Sardinian-related lineages
      admix1 = edgelist("R","Chimp.REF", "R","a", "a","Mbuti.DG", "a","b", "b","y", "b","s", "y","Russia_Samara_EBA_Yamnaya.AG", "y","fy", "s","Sardinian.DG", "s","fs", "fy","fr", "fs","fr", "fr","French.DG")
    )
  )
)

# qpgraph() fit at fixed admixture weights `w`, through the internals
# qpgraph() uses after its weight optimization (default diag, diag_f3).
qpgraph_fixed <- function(f2_blocks, edges, w) {
  g <- igraph::graph_from_edgelist(edges)
  pops <- admixtools:::get_leafnames(g)
  pre <- admixtools:::qpgraph_precompute_f3(f2_blocks, pops)
  nedges <- length(igraph::E(g))
  admixnodes <- which(igraph::degree(g, mode = "in") == 2)
  admixedgesfull <- sapply(seq_along(admixnodes),
                           function(i) igraph::incident_edges(g, admixnodes, mode = "in")[[i]][1:2])
  normedges <- setdiff(1:nedges, admixedgesfull)
  cmb <- combn(0:(length(pops) - 1), 2) + (1:0)
  pwts <- admixtools:::graph_to_pwts(g, pops)
  if (length(admixnodes) > 0) {
    wi <- admixtools:::graph_to_weightind(g)
    pwts <- admixtools:::cpp_fill_pwts(pwts, w, wi[[1]], wi[[2]], wi[[3]])
  }
  pwts <- pwts[, -1] - pwts[, 1]
  ppwts_2d <- t(pwts[, cmb[1, ]] * pwts[, cmb[2, ]])
  bl <- admixtools:::cpp_opt_edge_lengths(ppwts_2d, pre$ppinv, pre$f3_est, admixtools:::qpsolve,
                                         rep(0, length(normedges)),
                                         rep(.Machine$integer.max, length(normedges)), 1e-4, TRUE)
  f3_fit <- c(ppwts_2d %*% bl)
  weight <- rep(NA, nedges)
  weight[normedges] <- bl
  if (length(admixnodes) > 0) weight[c(t(admixedgesfull))] <- c(w, 1 - w)
  list(admix_weights = w, weight = weight, f3_fit = f3_fit,
       score = admixtools:::qpgraph_score(f3_fit, pre$f3_est, pre$ppinv))
}

# admixtools2 defaults, stated explicitly so the reference records them.
PARAMS <- list(blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
               auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2",
               apply_corr = TRUE)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
out_root <- file.path(script_dir, paste0("v", version))

run_config <- function(name, cfg) {
  cat("==", name, "\n")
  pref <- cfg$prefix
  pops <- cfg$pops
  PARAMS <- modifyList(PARAMS, if (is.null(cfg$params)) list() else cfg$params)

  # Ploidy, exactly as packedancestrymap_to_afs assigns it
  ind <- read.table(paste0(pref, ".ind"), col.names = c("id", "sex", "pop"), colClasses = "character")
  nsnp <- as.numeric(system(paste("wc -l <", shQuote(paste0(pref, ".snp"))), intern = TRUE))
  ip <- admixtools:::match_samples(ind$id, ind$pop, NULL, pops)
  indvec <- ip$indvec - 1
  ploidy <- admixtools:::cpp_packedancestrymap_ploidy(paste0(pref, ".geno"), nsnp, nrow(ind), indvec, 1000)
  sel <- which(indvec != -1)
  ploidy_by_sample <- setNames(as.list(ploidy[sel]), ind$id[sel])

  # SNP filtering stages, as in extract_f2
  afdat <- admixtools:::anygeno_to_aftable(pref, pops = pops, adjust_pseudohaploid = PARAMS$adjust_pseudohaploid,
                                           verbose = FALSE)
  n_read <- nrow(afdat$afs)
  afdat <- admixtools:::discard_from_aftable(afdat, maxmiss = PARAMS$maxmiss, minmaf = PARAMS$minmaf,
                                             maxmaf = PARAMS$maxmaf, auto_only = PARAMS$auto_only,
                                             poly_only = FALSE)
  n_kept <- nrow(afdat$afs)
  poly <- as.logical(admixtools:::cpp_is_polymorphic(afdat$afs))
  n_poly <- sum(poly)
  # Blocks before f2_from_geno(remove_na = TRUE) drops those where any pair
  # has no SNP with data in both populations (possible only with maxmiss > 0)
  block_lengths_all <- admixtools:::get_block_lengths(afdat$snpfile[poly, ], blgsize = PARAMS$blgsize,
                                                      verbose = FALSE)

  # The public pipeline
  f2_blocks <- f2_from_geno(pref, pops = pops, blgsize = PARAMS$blgsize, maxmiss = PARAMS$maxmiss,
                            minmaf = PARAMS$minmaf, maxmaf = PARAMS$maxmaf, auto_only = PARAMS$auto_only,
                            adjust_pseudohaploid = PARAMS$adjust_pseudohaploid,
                            poly_only = PARAMS$poly_only, apply_corr = PARAMS$apply_corr, verbose = FALSE)
  block_lengths <- readr::parse_number(dimnames(f2_blocks)[[3]])
  stopifnot(sum(block_lengths_all) == n_poly)

  blocks <- list()
  for (i in seq_along(pops)) for (j in seq_along(pops)) if (i < j) {
    blocks[[paste(pops[i], pops[j], sep = "|")]] <- unname(f2_blocks[pops[i], pops[j], ])
  }

  f2_stats <- f2(f2_blocks)

  out <- list(
    admixtools_version = version,
    prefix = pref,
    pops = pops,
    params = PARAMS,
    ploidy = ploidy_by_sample,
    n_snps_read = n_read,
    n_snps_kept = n_kept,
    n_snps_polymorphic = n_poly,
    block_lengths_before_remove_na = block_lengths_all,
    block_lengths = block_lengths,
    f2_blocks = blocks,
    f2 = f2_stats
  )
  if (identical(pref, AADR_PREFIX)) out$dataset <- AADR_DATASET
  dest <- file.path(out_root, name)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  write_json(out, file.path(dest, "f2.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE)

  fstats <- list(
    admixtools_version = version,
    f3 = select(f3(f2_blocks), pop1, pop2, pop3, est, se),
    f4 = select(f4(f2_blocks), pop1, pop2, pop3, pop4, est, se)
  )
  write_json(fstats, file.path(dest, "f3_f4.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("  %d f3, %d f4\n", nrow(fstats$f3), nrow(fstats$f4)))

  models <- lapply(cfg$qpwave, function(m) {
    st <- admixtools:::f2blocks_to_f4stats(f2_blocks, m$left, m$right)
    qw <- qpwave(f2_blocks, m$left, m$right, verbose = FALSE)
    # f4_est is (left[-1] x right[-1]); f4_var is ordered as c(t(f4_est)),
    # i.e. by left population first, then right.
    list(left = m$left, right = m$right,
         f4_est = unname(st$est), f4_var = unname(st$var),
         rankdrop = select(qw$rankdrop, f4rank, dof, chisq, p))
  })
  write_json(list(admixtools_version = version, models = models), file.path(dest, "qpwave.json"),
             digits = I(17), auto_unbox = TRUE, pretty = TRUE, matrix = "rowmajor")
  cat(sprintf("  %d qpwave models\n", length(models)))

  qpadm_models <- lapply(cfg$qpadm, function(m) {
    constrained <- isTRUE(m$constrained)
    if (constrained && length(m$left) >= 3) {
      # qpadm(constrained = TRUE) fails with three or more sources: drop_ranks
      # refits the lower ranks with the constrained solver, where the
      # quadratic program is singular and quadprog stops ("matrix D in
      # quadratic function is not positive definite"). The model's weights,
      # SEs and popdrop are computed with the same internals qpadm() calls
      # for them. Its rankdrop is drop_ranks() run unconstrained: rankdrop
      # (f4rank, dof, chisq, p) comes from the rank fit alone, and the
      # weights drop_ranks() computes along the way are discarded
      # (addweights = FALSE), so this is what qpadm() would report.
      err <- tryCatch({qpadm(f2_blocks, m$left, m$right, m$target, constrained = TRUE, verbose = FALSE); NA},
                      error = function(e) conditionMessage(e))
      st <- admixtools:::f2blocks_to_f4stats(f2_blocks, c(m$target, m$left), m$right)
      v <- st$var
      diag(v) <- diag(v) + 1e-4 * sum(diag(v))
      qinv <- solve(v)
      w <- admixtools:::cpp_qpadm_weights(st$est, qinv, length(m$left) - 1, 1e-4, constrained = TRUE,
                                          qpsolve = admixtools:::qpsolve)$weights
      se <- sqrt(diag(admixtools:::cpp_get_weights_covariance(st$f4_lo, qinv, st$block_lengths, 1e-4,
                                                              constrained = TRUE,
                                                              qpsolve = admixtools:::qpsolve)))
      rankdrop <- admixtools:::drop_ranks(st$est, qinv, 1e-4, constrained = FALSE, cpp = TRUE)
      popdrop <- admixtools:::drop_pops(st$est, qinv, 1e-4, constrained = TRUE, cpp = TRUE, m$left)
      return(list(target = m$target, left = m$left, right = m$right, constrained = TRUE,
                  qpadm_error = err, weights = tibble(left = m$left, weight = c(w), se = se),
                  rankdrop = select(rankdrop, f4rank, dof, chisq, p), popdrop = popdrop))
    }
    qa <- qpadm(f2_blocks, m$left, m$right, m$target, constrained = constrained, verbose = FALSE)
    # I(): a single source stays a JSON array
    list(target = m$target, left = I(m$left), right = I(m$right), constrained = constrained,
         weights = select(qa$weights, left, weight, se),
         rankdrop = select(qa$rankdrop, f4rank, dof, chisq, p),
         popdrop = qa$popdrop)
  })
  f4ratios <- bind_rows(lapply(cfg$f4ratio, function(p) qpf4ratio(f2_blocks, p)))
  write_json(list(admixtools_version = version, models = qpadm_models, f4ratio = f4ratios),
             file.path(dest, "qpadm.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("  %d qpadm models, %d f4ratio\n", length(qpadm_models), nrow(f4ratios)))

  qpgraph_fits <- lapply(names(cfg$qpgraph), function(gname) {
    el <- cfg$qpgraph[[gname]]
    fit <- qpgraph(f2_blocks, el, seed = 1)
    # the first incoming edge of each admixture node carries w
    g <- igraph::graph_from_edgelist(el)
    admixnodes <- which(igraph::degree(g, mode = "in") == 2)
    first_in <- vapply(seq_along(admixnodes),
                       function(i) as.integer(igraph::incident_edges(g, admixnodes, mode = "in")[[i]][1]), 1L)
    w_best <- fit$edges$weight[first_in]
    at_best <- qpgraph_fixed(f2_blocks, el, w_best)
    stopifnot(isTRUE(all.equal(at_best$score, fit$score, tolerance = 1e-12)),
              isTRUE(all.equal(at_best$weight, fit$edges$weight, tolerance = 1e-12)))
    list(name = gname, edges = el, pops = admixtools:::get_leafnames(g),
         score = fit$score, fit_edges = fit$edges, f3 = fit$f3, opt = fit$opt,
         ppinv = unname(fit$ppinv),
         fixed = list(at_best, qpgraph_fixed(f2_blocks, el, rep(0.3, length(admixnodes)))))
  })
  write_json(list(admixtools_version = version, graphs = qpgraph_fits), file.path(dest, "qpgraph.json"),
             digits = I(17), auto_unbox = TRUE, pretty = TRUE, matrix = "rowmajor", na = "null")
  cat(sprintf("  %d qpgraph fits\n", length(qpgraph_fits)))

  fpops <- cfg$fst_pops
  fst_blocks <- f2_from_geno(pref, pops = fpops, blgsize = PARAMS$blgsize, maxmiss = PARAMS$maxmiss,
                             auto_only = PARAMS$auto_only, adjust_pseudohaploid = PARAMS$adjust_pseudohaploid,
                             fst = TRUE, verbose = FALSE)
  fst_bl <- readr::parse_number(dimnames(fst_blocks)[[3]])
  fblocks <- list()
  for (i in seq_along(fpops)) for (j in seq_along(fpops)) if (i < j) {
    fblocks[[paste(fpops[i], fpops[j], sep = "|")]] <- unname(fst_blocks[fpops[i], fpops[j], ])
  }
  fst_stats <- fst(pref, pop1 = fpops, verbose = FALSE, blgsize = PARAMS$blgsize, maxmiss = PARAMS$maxmiss,
                   auto_only = PARAMS$auto_only, adjust_pseudohaploid = PARAMS$adjust_pseudohaploid)
  write_json(list(admixtools_version = version, pops = fpops, block_lengths = fst_bl, fst_blocks = fblocks,
                  fst = select(fst_stats, pop1, pop2, est, se)),
             file.path(dest, "fst.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("  %d fst pairs, %d fst blocks (%d SNPs)\n", nrow(fst_stats), length(fst_bl), sum(fst_bl)))
  cat(sprintf("  %d read, %d kept, %d polymorphic, %d blocks\n", n_read, n_kept, n_poly, length(block_lengths)))
  print(f2_stats)
}

selected <- if (length(args) == 0) names(CONFIGS) else args
for (name in selected) run_config(name, CONFIGS[[name]])
