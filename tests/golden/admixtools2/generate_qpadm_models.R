#!/usr/bin/env Rscript
#
# Reference for popmix's admixtools2 model-set functions: qpadm_p,
# test_cladality, qpwave_pairs, qpadm_rotate, rotate_models, qpadm_multi and
# qpadm_models.
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_qpadm_models.R
#
# Writes tests/golden/admixtools2/v<version>/dg_11pop/qpadm_models.json, on
# the dg_11pop f2 blocks of generate_reference.R (same populations and
# parameters). Numbers are "%.17g" strings, so they round-trip exactly.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
  library(purrr)
})
future::plan("sequential")

# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache;
# its manifest.json records the release files' checksums.
AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
AADR_PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
AADR_DATASET <- jsonlite::read_json(file.path(AADR_DIR, "manifest.json"))
POPS <- c("Druze.DG", "Basque.DG", "Sardinian.DG", "Yoruba.DG", "Han.DG", "Dinka.DG",
          "Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF")
RIGHT <- c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG", "Chimp.REF")

f2_blocks <- f2_from_geno(AADR_PREFIX, pops = POPS, blgsize = 0.05, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                          verbose = FALSE)

num17 <- function(x) sprintf("%.17g", x)
# a qpadm_p-style row: populations as lists, numbers as strings
row_out <- function(r) {
  r <- as.list(r)
  lapply(r, function(v) if (is.list(v)) unlist(v) else if (is.numeric(v)) num17(v) else v)
}
table_out <- function(d) lapply(seq_len(nrow(d)), function(i) row_out(d[i, ]))

qpadm_p_cases <- list(
  list(left = c("Sardinian.DG", "Dinka.DG", "Han.DG"), right = RIGHT, target = "Druze.DG"),
  list(left = c("Sardinian.DG", "Dinka.DG"), right = RIGHT, target = "Druze.DG"),
  list(left = c("Sardinian.DG", "Dinka.DG"), right = RIGHT, target = "Druze.DG", constrained = TRUE),
  list(left = c("Basque.DG", "Yoruba.DG"), right = RIGHT, target = "Druze.DG", rnk = 0),
  list(left = c("Sardinian.DG", "Basque.DG", "Druze.DG"), right = RIGHT, target = NULL)
)
qpadm_p_ref <- lapply(qpadm_p_cases, function(cs) {
  r <- qpadm_p(f2_blocks, cs$left, cs$right, cs$target,
               constrained = if (is.null(cs$constrained)) FALSE else cs$constrained,
               rnk = if (is.null(cs$rnk)) length(setdiff(cs$left, cs$target)) - 1 else cs$rnk, weights = TRUE)
  w <- r %>% select(-any_of(c("f4rank", "dof", "chisq", "p", "feasible")))
  c(cs, list(f4rank = r$f4rank, dof = r$dof, chisq = num17(r$chisq), p = num17(r$p), feasible = r$feasible,
             weights = if (ncol(w)) setNames(as.list(num17(unlist(w))), names(w)) else NULL))
})

cladality <- test_cladality(f2_blocks, c("Druze.DG", "Basque.DG"), RIGHT)
cladality <- list(left = c("Druze.DG", "Basque.DG"), right = RIGHT, chisq = num17(cladality$chisq),
                  p = num17(cladality$p), dof = cladality$dof)

pairs_left <- c("Druze.DG", "Basque.DG", "Sardinian.DG", "Han.DG", "Dinka.DG")
pairs <- qpwave_pairs(f2_blocks, pairs_left, RIGHT)

rot <- list(leftright = c("Sardinian.DG", "Dinka.DG", "Han.DG", "Basque.DG", "Yoruba.DG"),
            target = "Druze.DG", rightfix = c("Mbuti.DG", "Chimp.REF"))
rotate <- qpadm_rotate(f2_blocks, rot$leftright, rot$target, rot$rightfix, verbose = FALSE)
rmodels <- rotate_models(rot$leftright, rot$target, rot$rightfix)

models <- tibble(left = list(c("Sardinian.DG", "Dinka.DG"), c("Basque.DG", "Yoruba.DG", "Han.DG"), "Sardinian.DG"),
                 right = list(RIGHT, c("Mbuti.DG", "Papuan.DG", "Karitiana.DG", "Chimp.REF"), RIGHT),
                 target = c("Druze.DG", "Druze.DG", "Basque.DG"))
multi <- qpadm_multi(f2_blocks, models, full_results = FALSE, verbose = FALSE)
multi_full <- qpadm_multi(f2_blocks, models, full_results = TRUE, verbose = FALSE)
multi_full <- lapply(multi_full, function(m) list(
  weights = table_out(m$weights), rankdrop = table_out(m$rankdrop),
  popdrop = table_out(m$popdrop %>% select(pat, wt, dof, chisq, p, f4rank, feasible))))

qm <- function(pops, allpops, more_right) {
  d <- qpadm_models(pops, allpops = allpops, more_right = more_right)
  list(pops = pops, allpops = allpops, more_right = more_right,
       models = lapply(seq_len(nrow(d)), function(i) list(left = d$l[[i]], right = d$r[[i]])))
}
models_ref <- list(qm(c("A", "B", "C", "D", "E"), TRUE, TRUE), qm(c("A", "B", "C", "D"), FALSE, TRUE),
                   qm(c("A", "B", "C", "D"), TRUE, FALSE))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("v", packageVersion("admixtools")), "dg_11pop")
out <- list(admixtools_version = as.character(packageVersion("admixtools")), dataset = AADR_DATASET, pops = POPS,
            qpadm_p = qpadm_p_ref, test_cladality = cladality,
            qpwave_pairs = list(left = pairs_left, right = RIGHT, result = table_out(pairs)),
            qpadm_rotate = c(rot, list(result = table_out(rotate))),
            rotate_models = c(rot, list(models = table_out(rmodels))),
            qpadm_multi = list(models = table_out(models), result = table_out(multi), full = multi_full),
            qpadm_models = models_ref)
writeLines(toJSON(out, auto_unbox = TRUE, null = "null", pretty = FALSE), file.path(dest, "qpadm_models.json"))
cat("wrote", file.path(dest, "qpadm_models.json"), "\n")
