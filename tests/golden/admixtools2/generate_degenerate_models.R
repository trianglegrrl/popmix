#!/usr/bin/env Rscript
#
# admixtools2 qpadm() and qpwave() at the edges of what a model can be
# (TASKS.md K4): no degrees of freedom left, a single source, as many
# sources as right populations, and a population in both left and right.
#
# f2 blocks of AADR v66.p1 (f2_from_geno defaults), so the numbers are
# admixtools2's qpadm(f2_from_geno(prefix, pops)). A model admixtools2
# refuses is recorded with its error message.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_degenerate_models.R
#
# Writes tests/golden/admixtools2/v<version>/degenerate_models/models.json.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
dest <- file.path(script_dir, paste0("v", version), "degenerate_models")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

YAM <- "Russia_Samara_EBA_Yamnaya.AG"
LBK <- "Germany_HalberstadtSonntagsfeld_EN_LBK.AG"
POPS <- c("French.DG", YAM, LBK, "Mbuti.DG", "Han.DG", "Papuan.DG", "Karitiana.DG", "CEU.DG")
f2 <- f2_from_geno(PREFIX, pops = POPS, verbose = FALSE)

QPADM <- list(
  dof_zero = list(target = "French.DG", left = c(YAM, LBK), right = c("Mbuti.DG", "Han.DG")),
  one_source = list(target = "French.DG", left = "CEU.DG", right = c("Mbuti.DG", "Han.DG", "Papuan.DG")),
  sources_as_many_as_right = list(target = "French.DG", left = c(YAM, LBK, "CEU.DG"),
                                  right = c("Mbuti.DG", "Han.DG", "Papuan.DG")),
  source_in_right = list(target = "French.DG", left = c(YAM, LBK), right = c("Mbuti.DG", "Han.DG", "Papuan.DG", LBK)),
  target_in_right = list(target = "French.DG", left = c(YAM, LBK),
                         right = c("Mbuti.DG", "Han.DG", "Papuan.DG", "French.DG"))
)
QPWAVE <- list(
  dof_zero = list(left = c("French.DG", "CEU.DG"), right = c("Mbuti.DG", "Han.DG")),
  one_left = list(left = "French.DG", right = c("Mbuti.DG", "Han.DG", "Papuan.DG")),
  left_in_right = list(left = c("French.DG", "CEU.DG", "Han.DG"),
                       right = c("Mbuti.DG", "Han.DG", "Papuan.DG", "Karitiana.DG"))
)

attempt <- function(rec, f) {
  tryCatch(c(rec, f()), error = function(e) c(rec, list(error = conditionMessage(e))))
}
qpadm_out <- lapply(QPADM, function(m) attempt(m, function() {
  r <- qpadm(f2, m$left, m$right, m$target, verbose = FALSE)
  list(weights = r$weights$weight, se = r$weights$se, rankdrop = r$rankdrop)
}))
qpwave_out <- lapply(QPWAVE, function(m) attempt(m, function() {
  list(rankdrop = qpwave(f2, m$left, m$right, verbose = FALSE)$rankdrop)
}))
write_json(list(admixtools_version = version, prefix = "aadr_v66p1", populations = POPS,
                qpadm = qpadm_out, qpwave = qpwave_out),
           file.path(dest, "models.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("wrote", file.path(dest, "models.json"), "\n")
