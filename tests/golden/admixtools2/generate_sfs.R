#!/usr/bin/env Rscript
#
# Generate the admixtools2 reference for popmix's joint site-frequency
# spectrum.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_sfs.R [config...]
#
# admixtools2's joint_sfs() takes allele frequencies and counts
# (packedancestrymap_to_afs: every SNP, pseudo-haploid samples counted once)
# and counts, for each combination of per-population allele counts from 0 to
# the largest count seen, the SNPs with that combination of afs * counts. It
# documents joint_sfs(pref = prefix), but the geno_to_afs() that path calls
# does not exist in 2.0.10, so the reference calls packedancestrymap_to_afs()
# and passes its result.
#
# Writes tests/golden/admixtools2/v<version>/<config>/sfs.json with:
#   - pops, max_counts: the populations and the largest allele count of each
#   - sfs: one row per combination, in joint_sfs() order (the first
#     population varies slowest), with its count n
#   - n_snps: SNPs read; n_in_sfs: SNPs counted in the spectrum (SNPs where a
#     population has no data are left out, as are SNPs whose afs * counts is
#     not an exact integer)
#   - not_integer: SNPs (and the populations) whose afs * counts is not an
#     exact integer, which joint_sfs() drops

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"

CONFIGS <- list(
  family_3pop = list(prefix = "tests/golden/dates/family_chr20_22",
                     pops = c("Simulation", "Simulation.1", "Yoruba")),
  dg_3pop = list(prefix = AADR_PREFIX, pops = c("Mbuti.DG", "Han.DG", "Chimp.REF")),
  # pseudo-haploid samples with missing data
  ancient_2pop = list(prefix = AADR_PREFIX, pops = c("Russia_Samara_EBA_Yamnaya.AG", "French.DG"))
)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
out_root <- file.path(script_dir, paste0("v", version))

run_config <- function(name, cfg) {
  cat("==", name, "\n")
  afs <- admixtools:::packedancestrymap_to_afs(cfg$prefix, pops = cfg$pops, verbose = FALSE)
  stopifnot(identical(colnames(afs$afs), cfg$pops))
  sfs <- joint_sfs(afs)
  prod <- afs$afs * afs$counts
  off <- which(is.finite(prod) & prod != round(prod), arr.ind = TRUE)
  not_integer <- list(snp = unname(off[, 1]), pop = cfg$pops[off[, 2]],
                      value = prod[off], count = afs$counts[off])
  out <- list(
    admixtools_version = version,
    prefix = cfg$prefix,
    pops = I(cfg$pops),
    max_counts = unname(apply(afs$counts, 2, max)),
    n_snps = nrow(afs$afs),
    n_in_sfs = sum(sfs$n),
    sfs = lapply(as.list(sfs), unname),
    not_integer = not_integer
  )
  dest <- file.path(out_root, name)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  write_json(out, file.path(dest, "sfs.json"), digits = I(17), auto_unbox = TRUE)
  cat(sprintf("  %d SNPs, %d in the spectrum, %d cells, %d values not integer\n",
              out$n_snps, out$n_in_sfs, nrow(sfs), length(not_integer$snp)))
}

selected <- if (length(args) == 0) names(CONFIGS) else args
for (name in selected) run_config(name, CONFIGS[[name]])
