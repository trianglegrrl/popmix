#!/usr/bin/env Rscript
#
# ADMIXTOOLS 7.0.2 rexpfit.r's fit (affine = TRUE), printing the estimates
# with %.17g: DEoptim (seeded) picks the start, nls fits
# wcorr ~ C + A exp(-m dist / 100) to the column `col` of a ROLLOFF output
# file over lval <= dist <= hval (cM). Usage:
#   rexpfit_full.R input col lval hval seed
# Prints one line: C A m (or "NA NA NA" if nls fails).

args <- commandArgs(TRUE)
input <- args[1]; col <- as.numeric(args[2]); lval <- as.numeric(args[3]); hval <- as.numeric(args[4])
seed <- as.integer(args[5])
.libPaths(c("~/R/popmix-references", .libPaths()))
suppressPackageStartupMessages(library(DEoptim))

data <- read.table(input, header = FALSE)
dist <- data[, 1]
wcorr <- data[, col]
ndist <- length(dist)
data.sub <- data
if ((lval > dist[1]) || (hval < dist[ndist])) data.sub <- subset(data, ((dist <= hval) & (dist >= lval)))
dist <- data.sub[, 1]
wcorr <- data.sub[, col]

set.seed(seed)
fm1 <- function(x) x[1] + x[2] * exp(-x[3] * dist / 100)
fm2 <- function(x) sum((wcorr - fm1(x))^2)
fm3 <- DEoptim(fm2, lower = rep(0.0001, 3), upper = c(10e7, 10, 500), control = list(trace = FALSE))
par1 <- fm3$optim$bestmem
names(par1) <- c("C", "A", "m")
fit1 <- tryCatch(nls(wcorr ~ (C + A * exp(-m * dist / 100)), start = par1,
                     control = list(maxiter = 1000, warnOnly = TRUE)), error = function(e) NULL)
if (is.null(fit1)) cat("NA NA NA\n") else cat(sprintf("%.17g", coef(fit1)), "\n")
