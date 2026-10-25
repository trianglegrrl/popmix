#!/usr/bin/env Rscript
#
# Reference for popmix's R-compatible optimizer and random numbers.
#
# Usage: /usr/bin/Rscript tests/golden/rlang/generate_optim.R
#
# Writes tests/golden/rlang/r-<version>/optim.json:
#   - runif: set.seed(seed); runif(n) for several seeds
#   - sample: set.seed(seed), then sample(seq_len(n), size, replace = TRUE)
#     `reps` times
#   - lbfgsb: optim(method = "L-BFGS-B") on bounded test problems, with an
#     analytic gradient or R's finite differences (ndeps, fnscale, parscale),
#     with the result and every point the objective was evaluated at

suppressPackageStartupMessages(library(jsonlite))

rosen <- function(x) sum(100 * (x[-1] - x[-length(x)]^2)^2 + (1 - x[-length(x)])^2)
rosen_gr <- function(x) {
  n <- length(x)
  g <- numeric(n)
  g[-n] <- -400 * x[-n] * (x[-1] - x[-n]^2) - 2 * (1 - x[-n])
  g[-1] <- g[-1] + 200 * (x[-1] - x[-n]^2)
  g
}
quad <- function(x) sum((seq_along(x)) * (x - c(2, -1, 0.5, 3)[seq_along(x)])^2) + x[1] * x[2]
quad_gr <- function(x) {
  g <- 2 * seq_along(x) * (x - c(2, -1, 0.5, 3)[seq_along(x)])
  g[1] <- g[1] + x[2]
  g[2] <- g[2] + x[1]
  g
}
# a bumpy function in [0, 1]^2 with several local minima, like qpGraph's score
bumpy <- function(x) 1e-4 * (sin(7 * x[1]) * cos(5 * x[2]) + (x[1] - 0.3)^2 + 2 * (x[2] - 0.6)^2 + 1)

problems <- list(
  list(name = "rosen2_grad", fn = rosen, gr = rosen_gr, par = c(-1.2, 1), lower = c(-2, -2), upper = c(2, 2),
       control = list()),
  list(name = "rosen2_bound_grad", fn = rosen, gr = rosen_gr, par = c(-1.2, 1), lower = c(-2, 1.5), upper = c(2, 3),
       control = list()),
  list(name = "rosen5_grad", fn = rosen, gr = rosen_gr, par = c(-1.2, 1, -1.2, 1, 0.5),
       lower = rep(-0.5, 5), upper = rep(0.9, 5), control = list()),
  list(name = "quad4_halfbounded_grad", fn = quad, gr = quad_gr, par = c(0, 0, 0, 0),
       lower = c(-Inf, 0, -Inf, -Inf), upper = c(1, Inf, Inf, 2), control = list()),
  list(name = "rosen2_fd", fn = rosen, gr = NULL, par = c(-1.2, 1), lower = c(-2, -2), upper = c(2, 2),
       control = list()),
  list(name = "rosen5_fd_scaled", fn = rosen, gr = NULL, par = c(-1.2, 1, -1.2, 1, 0.5),
       lower = rep(-0.5, 5), upper = rep(0.9, 5), control = list(fnscale = 1e-6, maxit = 1e4)),
  list(name = "rosen2_fd_maxit", fn = rosen, gr = NULL, par = c(-1.2, 1), lower = c(-2, -2), upper = c(2, 2),
       control = list(maxit = 5)),
  list(name = "quad4_fd_lmm3", fn = quad, gr = NULL, par = c(0.5, 0.5, 0.5, 0.5),
       lower = c(-Inf, 0, -Inf, -Inf), upper = c(1, Inf, Inf, 2), control = list(lmm = 3, factr = 1e3)),
  list(name = "bumpy_fd_start1", fn = bumpy, gr = NULL, par = c(0.9, 0.1), lower = c(0, 0), upper = c(1, 1),
       control = list(maxit = 1e4, fnscale = 1e-6)),
  list(name = "bumpy_fd_start2", fn = bumpy, gr = NULL, par = c(0.05, 0.95), lower = c(0, 0), upper = c(1, 1),
       control = list(maxit = 1e4, fnscale = 1e-6)),
  list(name = "bumpy_fd_at_bound", fn = bumpy, gr = NULL, par = c(1, 0), lower = c(0, 0), upper = c(1, 1),
       control = list(maxit = 1e4, fnscale = 1e-6))
)

lbfgsb <- lapply(problems, function(p) {
  evals <- list()
  fn <- function(x) {
    evals[[length(evals) + 1]] <<- x
    p$fn(x)
  }
  res <- optim(p$par, fn, p$gr, method = "L-BFGS-B", lower = p$lower, upper = p$upper, control = p$control)
  list(name = p$name, gradient = if (is.null(p$gr)) "fd" else "analytic", par0 = p$par,
       lower = p$lower, upper = p$upper, control = p$control,
       par = res$par, value = res$value, counts = unname(res$counts), convergence = res$convergence,
       message = res$message, evals = do.call(rbind, evals))
})

# sample(replace = TRUE): R >= 3.6 draws each index by rejection from
# 16-bit chunks; sizes cover one chunk, two chunks (n > 2^16) and a power of 2
sample_ref <- lapply(list(list(seed = 1L, n = 710L, reps = 3L), list(seed = 42L, n = 5L, reps = 4L),
                          list(seed = 7L, n = 65536L, reps = 1L), list(seed = 3L, n = 100000L, reps = 1L),
                          list(seed = 9L, n = 1L, reps = 2L)), function(s) {
  set.seed(s$seed)
  c(s, list(values = lapply(seq_len(s$reps), function(i) I(sample(seq_len(s$n), min(s$n, 2000L), replace = TRUE)))))
})

runif_ref <- lapply(c(1L, 42L, -7L, 123456789L, 0L), function(s) {
  set.seed(s)
  list(seed = s, values = runif(1500))
})

# Numbers go to JSON as "%.17g" strings, which round-trip exactly (jsonlite's
# own 17-digit output does not), and keep Inf.
num17 <- function(x) sprintf("%.17g", x)
lbfgsb <- lapply(lbfgsb, function(p) {
  for (k in c("par0", "lower", "upper", "par", "value")) p[[k]] <- num17(p[[k]])
  p$control <- lapply(p$control, num17)
  p$evals <- lapply(seq_len(nrow(p$evals)), function(i) num17(p$evals[i, ]))
  p
})
runif_ref <- lapply(runif_ref, function(r) { r$values <- num17(r$values); r })

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
dest <- file.path(script_dir, paste0("r-", getRversion()))
dir.create(dest, showWarnings = FALSE)
writeLines(toJSON(list(r_version = R.version.string, runif = runif_ref, sample = sample_ref, lbfgsb = lbfgsb),
                  auto_unbox = TRUE, pretty = FALSE),
           file.path(dest, "optim.json"))
cat("wrote", file.path(dest, "optim.json"), "\n")
