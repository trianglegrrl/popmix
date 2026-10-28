#!/usr/bin/env Rscript
#
# admixtools2 qpadm() and qpwave() given a genotype prefix (TASKS.md H11):
# each f4 computed from the genotypes, on the SNPs every population of the
# model has data for (allsnps = FALSE, the default) or on those of its four
# populations (allsnps = TRUE), rather than from f2 blocks.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_qpadm_geno.R
#
# Writes tests/golden/admixtools2/v<version>/qpadm_geno/models.json: for each
# model and allsnps, the weights and SEs, rankdrop and popdrop of qpadm(), and
# (allsnps = FALSE only: qpwave() has no allsnps) the rank tests of qpwave()
# on the target and sources.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
})

AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
dest <- file.path(script_dir, paste0("v", version), "qpadm_geno")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

MODELS <- list(
  list(target = "French.DG", left = c("Russia_Samara_EBA_Yamnaya.AG", "Germany_HalberstadtSonntagsfeld_EN_LBK.AG"),
       right = c("Mbuti.DG", "Han.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG")),
  list(target = "Papuan.DG", left = c("Han.DG", "French.DG"), right = c("Chimp.REF", "Mbuti.DG", "Yoruba.DG")),
  list(target = "Basque.DG", left = c("France_Yonne_N.AG", "Russia_Samara_EBA_Yamnaya.AG", "Druze.DG"),
       right = c("Mbuti.DG", "Han.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Dinka.DG"))
)

out <- list()
for (m in MODELS) for (allsnps in c(FALSE, TRUE)) {
  r <- qpadm(PREFIX, m$left, m$right, m$target, allsnps = allsnps, verbose = FALSE)
  rec <- c(m, list(allsnps = allsnps, weights = r$weights$weight, se = r$weights$se, rankdrop = r$rankdrop,
                   popdrop = r$popdrop))
  # qpwave() takes no allsnps: its f4 are on the model's SNPs
  if (!allsnps) rec$qpwave <- qpwave(PREFIX, c(m$target, m$left), m$right, verbose = FALSE)$rankdrop
  out[[length(out) + 1]] <- rec
  cat(m$target, allsnps, r$rankdrop$chisq[1], "\n")
}
write_json(list(admixtools_version = version, prefix = "aadr_v66p1", models = out),
           file.path(dest, "models.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE, na = "null")
