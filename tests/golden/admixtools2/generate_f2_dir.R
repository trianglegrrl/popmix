#!/usr/bin/env Rscript
#
# Reference for popmix's admixtools2 f2 directories: admixtools2
# extract_f2() of the dg_6pop populations and parameters (generate_reference.R)
# into tests/golden/admixtools2/v<version>/dg_6pop/f2_dir.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_f2_dir.R

suppressPackageStartupMessages(library(admixtools))

# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache;
AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
AADR_PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
POPS <- c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF")
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_6pop", "f2_dir")
unlink(dest, recursive = TRUE)
extract_f2(AADR_PREFIX, dest, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
           auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE, verbose = FALSE)
cat("wrote", dest, "\n")
