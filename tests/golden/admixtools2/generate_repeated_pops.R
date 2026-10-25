#!/usr/bin/env Rscript
#
# Generate the admixtools2 reference for f3 and f4 from f2 blocks when a
# population appears more than once (f3(A; A, B), f3(A; B, B),
# f4(A, B; A, C), f4(A, A; B, C), f4(A, B; A, B)): admixtools2 sets a
# population's f2 with itself to 0, so f3(A; A, B) = 0, f3(A; B, B) = f2(A, B),
# f4(A, B; A, C) = f3(A; B, C) and f4(A, B; A, B) = f2(A, B).
#
# The small family data (in the repository), f2 blocks with the settings of
# generate_reference.R. One f3() or f4() call per combination, so none is
# dropped as a duplicate.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_repeated_pops.R
#
# Writes tests/golden/admixtools2/v<version>/family_4pop/repeated_pops.json.
# Numbers are written with 17 significant digits.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
pops <- c("French", "Yoruba", "Simulation", "Simulation.1")
f2_blocks <- f2_from_geno("tests/golden/dates/family_chr20_22", pops = pops, blgsize = 0.05, maxmiss = 0,
                          minmaf = 0, maxmaf = 0.5, auto_only = TRUE, adjust_pseudohaploid = TRUE,
                          poly_only = "f2", apply_corr = TRUE, verbose = FALSE)

num <- function(x) sprintf("%.17g", x)
f3s <- list(c("French", "French", "Yoruba"), c("French", "Yoruba", "French"), c("Yoruba", "French", "French"),
            c("Simulation.1", "Simulation.1", "French"), c("French", "Yoruba", "Yoruba"),
            c("French", "Yoruba", "Simulation"))
f4s <- list(c("French", "Yoruba", "French", "Simulation"), c("French", "Yoruba", "Simulation", "French"),
            c("French", "French", "Yoruba", "Simulation"), c("French", "Yoruba", "French", "Yoruba"),
            c("Simulation", "French", "Yoruba", "Simulation.1"))
row <- function(r, n) c(as.list(setNames(unlist(r[1, paste0("pop", seq_len(n))]), paste0("pop", seq_len(n)))),
                        list(est = num(r$est), se = num(r$se)))
out <- list(
  admixtools_version = version,
  f3 = lapply(f3s, function(p) row(f3(f2_blocks, p[1], p[2], p[3]), 3)),
  f4 = lapply(f4s, function(p) row(f4(f2_blocks, p[1], p[2], p[3], p[4]), 4))
)
dest <- file.path(script_dir, paste0("v", version), "family_4pop", "repeated_pops.json")
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE), dest)
cat("wrote", dest, "\n")
