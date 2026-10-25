#!/usr/bin/env Rscript
#
# Reference for popmix's admixtools2 f2 directories: admixtools2
# extract_f2() of the dg_6pop populations and parameters (generate_reference.R)
# into tests/golden/admixtools2/v<version>/dg_6pop/f2_dir.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_f2_dir.R

suppressPackageStartupMessages(library(admixtools))

AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"
POPS <- c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF")
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop", "f2_dir")
unlink(dest, recursive = TRUE)
extract_f2(AADR_PREFIX, dest, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
           auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE, verbose = FALSE)
cat("wrote", dest, "\n")
