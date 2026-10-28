#!/usr/bin/env Rscript
#
# admixtools2 extract_f2's SNP filters (TASKS.md H1): transitions,
# transversions, minmaf/maxmaf, keepsnps, each with maxmiss.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_snp_filters.R
#
# Writes tests/golden/admixtools2/v<version>/snp_filters/<case>.json: the
# block lengths (polymorphic SNPs per block), every pair's f2 blocks and
# f2() estimate and SE, and fst() for every pair, from f2_from_geno and
# fst() given the case's filters; and keepsnps.txt, the SNP IDs the keepsnps
# cases keep (every 20th SNP of chromosomes 1 and 2 and the first 1000 of
# chromosome 23, which keepsnps keeps despite auto_only).

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
POPS <- c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG")

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
dest <- file.path(script_dir, paste0("v", version), "snp_filters")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

snp <- read.table(paste0(PREFIX, ".snp"), colClasses = c("character", "character", "numeric", "numeric",
                                                         "character", "character"))
keep <- c(snp$V1[snp$V2 %in% c("1", "2")][c(TRUE, rep(FALSE, 19))], snp$V1[snp$V2 == "23"][1:1000])
writeLines(keep, file.path(dest, "keepsnps.txt"))

CASES <- list(
  transversions = list(transitions = FALSE),
  transitions = list(transversions = FALSE),
  maf = list(minmaf = 0.05, maxmaf = 0.4),
  transversions_maf_maxmiss = list(transitions = FALSE, minmaf = 0.05, maxmiss = 0.25),
  keepsnps = list(keepsnps = keep),
  keepsnps_maf = list(keepsnps = keep, minmaf = 0.2, transitions = FALSE)
)

for (name in names(CASES)) {
  args <- CASES[[name]]
  f2_blocks <- do.call(f2_from_geno, c(list(PREFIX, pops = POPS, verbose = FALSE), args))
  block_lengths <- readr::parse_number(dimnames(f2_blocks)[[3]])
  blocks <- list()
  for (i in seq_along(POPS)) for (j in seq_along(POPS)) if (i < j) {
    blocks[[paste(POPS[i], POPS[j], sep = "|")]] <- unname(f2_blocks[POPS[i], POPS[j], ])
  }
  fst_stats <- do.call(fst, c(list(PREFIX, pop1 = POPS, verbose = FALSE), args))
  params <- args
  if (!is.null(params$keepsnps)) params$keepsnps <- "keepsnps.txt"
  out <- list(admixtools_version = version, pops = POPS, params = params, block_lengths = block_lengths,
              f2_blocks = blocks, f2 = f2(f2_blocks), fst = select(fst_stats, pop1, pop2, est, se))
  write_json(out, file.path(dest, paste0(name, ".json")), digits = I(17), auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("%s: %d SNPs in %d blocks\n", name, sum(block_lengths), length(block_lengths)))
}
