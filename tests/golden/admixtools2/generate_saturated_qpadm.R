#!/usr/bin/env Rscript
#
# Generate the admixtools2 reference for a saturated qpAdm model: as many
# right populations as left, so the full-rank fit has 0 degrees of freedom
# and nothing is tested. admixtools2 reports its p as NA (popdrop and
# rankdrop), not 1.
#
# French.DG = Han.DG + Yoruba.DG with right Mbuti.DG, Chimp.REF, on the
# dg_6pop admixtools2 f2 directory (generate_f2_dir.R).
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_saturated_qpadm.R
#
# Writes tests/golden/admixtools2/v<version>/dg_6pop/qpadm_saturated.json.
# Numbers are written with 17 significant digits; NA as "NA".

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
dir <- file.path(script_dir, paste0("v", version), "dg_6pop")
f2 <- read_f2(file.path(dir, "f2_dir"))
model <- list(target = "French.DG", left = c("Han.DG", "Yoruba.DG"), right = c("Mbuti.DG", "Chimp.REF"))
r <- qpadm(f2, left = model$left, right = model$right, target = model$target, verbose = FALSE)

num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
out <- c(model, list(
  admixtools_version = version,
  weights = lapply(seq_len(nrow(r$weights)), function(i) list(left = r$weights$left[i],
    weight = num(r$weights$weight[i]), se = num(r$weights$se[i]))),
  rankdrop = lapply(seq_len(nrow(r$rankdrop)), function(i) list(f4rank = r$rankdrop$f4rank[i],
    dof = r$rankdrop$dof[i], chisq = num(r$rankdrop$chisq[i]), p = num(r$rankdrop$p[i]))),
  popdrop = lapply(seq_len(nrow(r$popdrop)), function(i) list(pat = r$popdrop$pat[i],
    dof = r$popdrop$dof[i], chisq = num(r$popdrop$chisq[i]), p = num(r$popdrop$p[i]),
    feasible = r$popdrop$feasible[i]))
))
dest <- file.path(dir, "qpadm_saturated.json")
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE), dest)
cat("wrote", dest, "\n")
