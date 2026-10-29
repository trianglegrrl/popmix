#!/usr/bin/env Rscript
#
# admixtools2 f2_from_geno's SNP filters on the small family data
# (tests/golden/dates/family_chr20_22, tests/golden/dates/README.md): the
# polymorphic SNPs in its blocks with each filter, which
# tests/python/test_snp_filters.py checks popmix's extract_f2 against.
# Public API only.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_family_snp_filters.R
#
# Writes tests/golden/admixtools2/v<version>/family_4pop/snp_filters.json.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

PREFIX <- "tests/golden/dates/family_chr20_22"
POPS <- c("French", "Yoruba", "Simulation", "Simulation.1")
CASES <- list(
  all = list(),
  transversions = list(transitions = FALSE),
  transitions = list(transversions = FALSE),
  minmaf_0.1 = list(minmaf = 0.1),
  maxmaf_0.2 = list(maxmaf = 0.2),
  maf_0.05_0.3 = list(minmaf = 0.05, maxmaf = 0.3)
)

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
out <- list(admixtools_version = version, prefix = PREFIX, pops = POPS, cases = list())
for (name in names(CASES)) {
  f2_blocks <- do.call(f2_from_geno, c(list(PREFIX, pops = POPS, verbose = FALSE), CASES[[name]]))
  n <- sum(readr::parse_number(dimnames(f2_blocks)[[3]]))
  out$cases[[name]] <- list(params = CASES[[name]], snps = n)
  cat(sprintf("%s: %d SNPs\n", name, n))
}
dest <- file.path(script_dir, paste0("v", version), "family_4pop", "snp_filters.json")
write_json(out, dest, auto_unbox = TRUE, pretty = TRUE)
