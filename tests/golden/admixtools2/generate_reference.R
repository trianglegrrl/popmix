#!/usr/bin/env Rscript
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

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"

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
    f4ratio = list(c("Chimp.REF", "Mbuti.DG", "Druze.DG", "Dinka.DG", "Sardinian.DG"))
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
    )
  )
)

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
      # quadratic function is not positive definite"). The model's weights
      # and SEs are computed with the same internals qpadm() calls for them.
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
      return(list(target = m$target, left = m$left, right = m$right, constrained = TRUE,
                  qpadm_error = err, weights = tibble(left = m$left, weight = c(w), se = se)))
    }
    qa <- qpadm(f2_blocks, m$left, m$right, m$target, constrained = constrained, verbose = FALSE)
    list(target = m$target, left = m$left, right = m$right, constrained = constrained,
         weights = select(qa$weights, left, weight, se),
         rankdrop = select(qa$rankdrop, f4rank, dof, chisq, p),
         popdrop = qa$popdrop)
  })
  f4ratios <- bind_rows(lapply(cfg$f4ratio, function(p) qpf4ratio(f2_blocks, p)))
  write_json(list(admixtools_version = version, models = qpadm_models, f4ratio = f4ratios),
             file.path(dest, "qpadm.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("  %d qpadm models, %d f4ratio\n", length(qpadm_models), nrow(f4ratios)))

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
