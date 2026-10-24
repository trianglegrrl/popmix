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
# Numbers are written with 17 significant digits, enough to round-trip a
# double exactly.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"

CONFIGS <- list(
  dg_6pop = list(
    prefix = AADR_PREFIX,
    pops = c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF"),
    qpwave = list(
      list(left = c("Han.DG", "French.DG"), right = c("Yoruba.DG", "Mbuti.DG", "Chimp.REF")),
      list(left = c("French.DG", "Han.DG", "Papuan.DG"), right = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG"))
    )
  ),
  ancient_5pop = list(
    prefix = AADR_PREFIX,
    pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG", "Chimp.REF"),
    qpwave = list(
      list(left = c("French.DG", "Sardinian.DG", "Russia_Samara_EBA_Yamnaya.AG"),
           right = c("Mbuti.DG", "Chimp.REF"))
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
  n_poly <- sum(admixtools:::cpp_is_polymorphic(afdat$afs))

  # The public pipeline
  f2_blocks <- f2_from_geno(pref, pops = pops, blgsize = PARAMS$blgsize, maxmiss = PARAMS$maxmiss,
                            minmaf = PARAMS$minmaf, maxmaf = PARAMS$maxmaf, auto_only = PARAMS$auto_only,
                            adjust_pseudohaploid = PARAMS$adjust_pseudohaploid,
                            poly_only = PARAMS$poly_only, apply_corr = PARAMS$apply_corr, verbose = FALSE)
  block_lengths <- readr::parse_number(dimnames(f2_blocks)[[3]])
  stopifnot(sum(block_lengths) == n_poly)

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
  cat(sprintf("  %d read, %d kept, %d polymorphic, %d blocks\n", n_read, n_kept, n_poly, length(block_lengths)))
  print(f2_stats)
}

selected <- if (length(args) == 0) names(CONFIGS) else args
for (name in selected) run_config(name, CONFIGS[[name]])
