#!/usr/bin/env Rscript
#
# Reference for popmix's lazadm, qpf4diff and qpfstats.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_more_stats.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/more_stats.json, on the
# dg_6pop f2 blocks of generate_reference.R (and its genotypes, for
# qpfstats). Numbers are "%.17g" strings.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"
POPS <- c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF")
f2_blocks <- f2_from_geno(AADR_PREFIX, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                          verbose = FALSE)
num17 <- function(x) sprintf("%.17g", x)
tab <- function(d) {
  d <- as.data.frame(d)
  for (k in names(d)) if (is.numeric(d[[k]])) d[[k]] <- num17(d[[k]])
  unname(lapply(split(d, seq_len(nrow(d))), as.list))
}

LAZ <- list(target = "French.DG", left = c("Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Papuan.DG", "Chimp.REF"))
lazadm_ref <- lapply(c(TRUE, FALSE), function(con) c(LAZ, list(constrained = con, result = tab(
  lazadm(f2_blocks, LAZ$left, LAZ$right, LAZ$target, constrained = con)))))

F4DIFF <- rbind(c("Chimp.REF", "Mbuti.DG", "French.DG", "Han.DG", "Yoruba.DG"),
                c("Mbuti.DG", "Yoruba.DG", "Papuan.DG", "Chimp.REF", "Han.DG"))
f4diff_ref <- tab(qpf4diff(f2_blocks, F4DIFF))

QPF <- c("Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG")
qf <- qpfstats(AADR_PREFIX, QPF, verbose = FALSE)
qpfstats_ref <- list(pops = QPF, sorted_pops = dimnames(qf)[[1]], block_lengths = as.numeric(gsub("l", "", dimnames(qf)[[3]])),
                     f2 = lapply(setNames(nm = combn(dimnames(qf)[[1]], 2, paste, collapse = "|")), function(p) {
                       ab <- strsplit(p, "\\|")[[1]]; num17(qf[ab[1], ab[2], ])
                     }))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop")
writeLines(toJSON(list(admixtools_version = as.character(packageVersion("admixtools")), lazadm = lazadm_ref,
                       qpf4diff = list(pops = lapply(seq_len(nrow(F4DIFF)), function(i) F4DIFF[i, ]), result = f4diff_ref),
                       qpfstats = qpfstats_ref), auto_unbox = TRUE, pretty = FALSE),
           file.path(dest, "more_stats.json"))
cat("wrote", file.path(dest, "more_stats.json"), "\n")
