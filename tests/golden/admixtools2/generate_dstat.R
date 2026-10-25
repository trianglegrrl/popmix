#!/usr/bin/env Rscript
#
# Generate the admixtools2 reference for f4-, D- and f3-statistics computed
# directly from genotype files (qpdstat() and qp3pop() given a genotype
# prefix rather than f2 blocks).
#
# Usage: /usr/bin/Rscript tests/golden/admixtools2/generate_dstat.R [config...]
#
# Writes tests/golden/admixtools2/v<version>/<config>/dstat.json with:
#   - qpdstat_f4, qpdstat_d: qpdstat(prefix, pops) with f4mode TRUE and FALSE,
#     for every unique combination of the config's populations, and
#     qpdstat_grid_d: qpdstat(prefix, pop1, pop2, pop3, pop4, f4mode = FALSE)
#     for vectors (expand_grid order, repeated populations included)
#   - qp3pop, qp3pop_outgroup: qp3pop(prefix, pops) with outgroupmode FALSE
#     (the default: f3 normalized by the target's heterozygosity) and TRUE
#   - f4_blocks, f3_blocks: the per-block data (f4blockdat_from_geno and
#     f3blockdat_from_geno, as qpdstat and qp3pop call them) for a few of the
#     combinations, so a mismatch can be traced to a block
#   - ploidy: the ploidy admixtools2 assigns each sample (used by qp3pop)
#
# These differ from the f-statistics computed from f2 blocks
# (generate_reference.R): every autosomal SNP is used (no polymorphism or
# maxmiss filter), blocks are formed over all of them, each statistic uses the
# SNPs with data in its own populations (allsnps = TRUE), and allele
# frequencies count every sample as diploid (qpdstat) or adjust for
# pseudo-haploid samples only in the f3 bias correction (qp3pop).
#
# The "mixed_ploidy" config relabels AADR samples so that one population mixes
# pseudo-haploid and diploid samples; its .ind is written next to symlinks to
# the AADR .geno and .snp, and the relabeling is recorded in the reference so
# the tests can recreate it.
#
# Numbers are written with 17 significant digits.

suppressPackageStartupMessages({
  library(admixtools)
  library(jsonlite)
  library(dplyr)
})

# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache;
# its manifest.json records the release files' checksums.
AADR_DIR <- file.path(Sys.getenv("POPMIX_REFERENCE_CACHE", path.expand("~/.cache/popmix-references")), "aadr_v66p1")
AADR_PREFIX <- file.path(AADR_DIR, "aadr_v66p1")
AADR_DATASET <- jsonlite::read_json(file.path(AADR_DIR, "manifest.json"))

CONFIGS <- list(
  family_4pop = list(
    prefix = "tests/golden/dates/family_chr20_22",
    pops = c("French", "Yoruba", "Simulation", "Simulation.1"),
    grid = list(pop1 = c("Simulation", "French"), pop2 = c("Yoruba", "French"),
                pop3 = "Simulation.1", pop4 = c("Yoruba", "Simulation"))
  ),
  dg_6pop = list(
    prefix = AADR_PREFIX,
    pops = c("Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF"),
    grid = list(pop1 = c("Han.DG", "French.DG"), pop2 = c("Papuan.DG", "Han.DG"),
                pop3 = c("Mbuti.DG", "Yoruba.DG"), pop4 = "Chimp.REF")
  ),
  ancient_5pop = list(
    prefix = AADR_PREFIX,
    pops = c("Russia_Samara_EBA_Yamnaya.AG", "Sardinian.DG", "French.DG", "Mbuti.DG", "Chimp.REF"),
    grid = list(pop1 = "Russia_Samara_EBA_Yamnaya.AG", pop2 = c("Sardinian.DG", "French.DG"),
                pop3 = "Mbuti.DG", pop4 = "Chimp.REF")
  ),
  # pseudo-haploid Samara Yamnaya and five diploid shotgun Yamnaya genomes
  # (v62.0's Russia_Kalmykia_EBA_Yamnaya.SG and Russia_Remontnoye_EBA_Yamnaya.SG,
  # relabelled in v66.p1)
  mixed_ploidy = list(
    prefix = AADR_PREFIX,
    relabel = list(Yamnaya_mixed = c("Russia_Samara_EBA_Yamnaya.AG", "Russia_EBA_Kalmykia_Yamnaya.SG",
                                     "Russia_EBA_Yamnaya.SG")),
    pops = c("Yamnaya_mixed", "French.DG", "Mbuti.DG", "Han.DG"),
    grid = list(pop1 = "Yamnaya_mixed", pop2 = "French.DG", pop3 = "Mbuti.DG", pop4 = "Han.DG")
  )
)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
version <- as.character(packageVersion("admixtools"))
out_root <- file.path(script_dir, paste0("v", version))

# A prefix whose .ind gives the samples of each population in `relabel` the
# population's name; .geno and .snp are symlinks to the original files.
relabeled_prefix <- function(prefix, relabel) {
  dir <- file.path(tempdir(), "relabeled")
  dir.create(dir, showWarnings = FALSE)
  out <- file.path(dir, basename(prefix))
  for (ext in c(".geno", ".snp")) file.symlink(normalizePath(paste0(prefix, ext)), paste0(out, ext))
  ind <- read.table(paste0(prefix, ".ind"), col.names = c("id", "sex", "pop"), colClasses = "character")
  for (name in names(relabel)) ind$pop[ind$pop %in% relabel[[name]]] <- name
  write.table(ind, paste0(out, ".ind"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  out
}

sample_ploidy <- function(pref, pops) {
  ind <- read.table(paste0(pref, ".ind"), col.names = c("id", "sex", "pop"), colClasses = "character")
  nsnp <- as.numeric(system(paste("wc -l <", shQuote(paste0(pref, ".snp"))), intern = TRUE))
  sel <- ind$pop %in% pops
  ploidy <- admixtools:::cpp_packedancestrymap_ploidy(paste0(pref, ".geno"), nsnp, nrow(ind), sel + 0, 1000)
  setNames(as.list(ploidy[sel]), ind$id[sel])
}

run_config <- function(name, cfg) {
  cat("==", name, "\n")
  pref <- cfg$prefix
  if (!is.null(cfg$relabel)) pref <- relabeled_prefix(pref, cfg$relabel)
  pops <- cfg$pops
  g <- cfg$grid

  t0 <- Sys.time()
  f4 <- qpdstat(pref, pops, verbose = FALSE)
  d <- qpdstat(pref, pops, f4mode = FALSE, verbose = FALSE)
  d_grid <- qpdstat(pref, g$pop1, g$pop2, g$pop3, g$pop4, f4mode = FALSE, verbose = FALSE)
  f3 <- qp3pop(pref, pops)
  f3_out <- qp3pop(pref, pops, outgroupmode = TRUE)
  cat(sprintf("  %d f4, %d D, %d grid D, %d f3 (%.0f s)\n", nrow(f4), nrow(d), nrow(d_grid), nrow(f3),
              as.numeric(Sys.time() - t0, units = "secs")))

  # Per-block data, exactly as qpdstat_geno and qp3pop_geno compute it
  # (allsnps = TRUE); est is the per-block mean, NA for blocks with no SNP.
  f4_combs <- d %>% select(pop1:pop4) %>% slice(1:2)
  f4_blocks <- admixtools:::f4blockdat_from_geno(pref, f4_combs, f4mode = FALSE, allsnps = TRUE, verbose = FALSE)
  f3_combs <- f3 %>% select(pop1:pop3) %>% slice(1:2)
  f3_blocks <- admixtools:::f3blockdat_from_geno(pref, f3_combs, allsnps = TRUE, verbose = FALSE)

  out <- list(
    admixtools_version = version,
    prefix = cfg$prefix,
    relabel = if (is.null(cfg$relabel)) NULL else lapply(cfg$relabel, I),
    pops = pops,
    grid = lapply(g, I),
    ploidy = sample_ploidy(pref, pops),
    qpdstat_f4 = f4,
    qpdstat_d = d,
    qpdstat_grid_d = d_grid,
    qp3pop = f3,
    qp3pop_outgroup = f3_out,
    f4_blocks = f4_blocks,
    f3_blocks = f3_blocks
  )
  if (identical(cfg$prefix, AADR_PREFIX)) out$dataset <- AADR_DATASET
  dest <- file.path(out_root, name)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  write_json(out, file.path(dest, "dstat.json"), digits = I(17), auto_unbox = TRUE, pretty = TRUE,
             na = "null")
  print(head(d))
  print(head(f3))
}

selected <- if (length(args) == 0) names(CONFIGS) else args
for (name in selected) run_config(name, CONFIGS[[name]])
