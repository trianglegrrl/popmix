#!/usr/bin/env Rscript
#
# Generate golden master reference data from admixtools2.
#
# Usage: Rscript generate_reference.R
#
# This script runs admixtools2 on test data and saves the outputs
# as JSON files for comparison with popmix.

library(admixtools)
library(jsonlite)

# Test data path
AADR_PREFIX <- "/references/aadr/v62.0_1240k_public"

# Output directory
OUTPUT_DIR <- dirname(sys.frame(1)$ofile)
if (is.null(OUTPUT_DIR)) OUTPUT_DIR <- "."

cat("Generating golden master data from admixtools2...\n")

# Define test populations (use actual AADR population names)
# These are the ".DG" suffix populations which are the high-quality ones
TEST_POPS <- list(
  han = "Han.DG",
  yoruba = "Yoruba.DG",
  mbuti = "Mbuti.DG",
  french = "French.DG"
)

# Check if populations exist in the data
cat("Loading genotype data...\n")
geno <- read_f2(AADR_PREFIX, pops = unlist(TEST_POPS), maxmem = 4000)

cat("Available populations: ", length(unique(geno$pops)), "\n")

# ============================================================================
# Test 1: F2 statistics
# ============================================================================
cat("\n=== Generating F2 reference data ===\n")

f2_results <- list()

# F2 between Han and Yoruba
f2_han_yoruba <- f2(geno, pop1 = TEST_POPS$han, pop2 = TEST_POPS$yoruba)
f2_results$han_yoruba <- list(
  pop1 = TEST_POPS$han,
  pop2 = TEST_POPS$yoruba,
  estimate = f2_han_yoruba$est,
  std_err = f2_han_yoruba$se,
  z_score = f2_han_yoruba$z
)
cat("f2(Han, Yoruba) = ", f2_han_yoruba$est, " +/- ", f2_han_yoruba$se, "\n")

# F2 between Han and Mbuti
f2_han_mbuti <- f2(geno, pop1 = TEST_POPS$han, pop2 = TEST_POPS$mbuti)
f2_results$han_mbuti <- list(
  pop1 = TEST_POPS$han,
  pop2 = TEST_POPS$mbuti,
  estimate = f2_han_mbuti$est,
  std_err = f2_han_mbuti$se,
  z_score = f2_han_mbuti$z
)
cat("f2(Han, Mbuti) = ", f2_han_mbuti$est, " +/- ", f2_han_mbuti$se, "\n")

# Save F2 reference
write_json(f2_results, file.path(OUTPUT_DIR, "f2_reference.json"), pretty = TRUE)

# ============================================================================
# Test 1b: FST statistics
# ============================================================================
cat("\n=== Generating FST reference data ===\n")

fst_results <- list()

# FST between Han and Yoruba
fst_han_yoruba <- fst(geno, pop1 = TEST_POPS$han, pop2 = TEST_POPS$yoruba)
fst_results$han_yoruba <- list(
  pop1 = TEST_POPS$han,
  pop2 = TEST_POPS$yoruba,
  estimate = fst_han_yoruba$est,
  std_err = fst_han_yoruba$se
)
cat("FST(Han, Yoruba) = ", fst_han_yoruba$est, " +/- ", fst_han_yoruba$se, "\n")

# FST between Han and Mbuti
fst_han_mbuti <- fst(geno, pop1 = TEST_POPS$han, pop2 = TEST_POPS$mbuti)
fst_results$han_mbuti <- list(
  pop1 = TEST_POPS$han,
  pop2 = TEST_POPS$mbuti,
  estimate = fst_han_mbuti$est,
  std_err = fst_han_mbuti$se
)
cat("FST(Han, Mbuti) = ", fst_han_mbuti$est, " +/- ", fst_han_mbuti$se, "\n")

# FST between Yoruba and Mbuti (closer populations)
fst_yoruba_mbuti <- fst(geno, pop1 = TEST_POPS$yoruba, pop2 = TEST_POPS$mbuti)
fst_results$yoruba_mbuti <- list(
  pop1 = TEST_POPS$yoruba,
  pop2 = TEST_POPS$mbuti,
  estimate = fst_yoruba_mbuti$est,
  std_err = fst_yoruba_mbuti$se
)
cat("FST(Yoruba, Mbuti) = ", fst_yoruba_mbuti$est, " +/- ", fst_yoruba_mbuti$se, "\n")

# Save FST reference
write_json(fst_results, file.path(OUTPUT_DIR, "fst_reference.json"), pretty = TRUE)

# ============================================================================
# Test 2: F3 statistics
# ============================================================================
cat("\n=== Generating F3 reference data ===\n")

f3_results <- list()

# F3(French; Han, Yoruba) - admixture test
f3_french <- f3(geno, pop1 = TEST_POPS$french, pop2 = TEST_POPS$han, pop3 = TEST_POPS$yoruba)
f3_results$french_han_yoruba <- list(
  target = TEST_POPS$french,
  source1 = TEST_POPS$han,
  source2 = TEST_POPS$yoruba,
  estimate = f3_french$est,
  std_err = f3_french$se,
  z_score = f3_french$z
)
cat("f3(French; Han, Yoruba) = ", f3_french$est, " +/- ", f3_french$se, "\n")

# Save F3 reference
write_json(f3_results, file.path(OUTPUT_DIR, "f3_reference.json"), pretty = TRUE)

# ============================================================================
# Test 3: F4 statistics
# ============================================================================
cat("\n=== Generating F4 reference data ===\n")

f4_results <- list()

# F4(French, Han; Yoruba, Mbuti)
f4_result <- f4(geno,
                pop1 = TEST_POPS$french, pop2 = TEST_POPS$han,
                pop3 = TEST_POPS$yoruba, pop4 = TEST_POPS$mbuti)
f4_results$french_han_yoruba_mbuti <- list(
  pop1 = TEST_POPS$french,
  pop2 = TEST_POPS$han,
  pop3 = TEST_POPS$yoruba,
  pop4 = TEST_POPS$mbuti,
  estimate = f4_result$est,
  std_err = f4_result$se,
  z_score = f4_result$z
)
cat("f4(French, Han; Yoruba, Mbuti) = ", f4_result$est, " +/- ", f4_result$se, "\n")

# Save F4 reference
write_json(f4_results, file.path(OUTPUT_DIR, "f4_reference.json"), pretty = TRUE)

# ============================================================================
# Test 4: qpWave
# ============================================================================
cat("\n=== Generating qpWave reference data ===\n")

qpwave_results <- list()

# qpWave with French, Han as left; Yoruba, Mbuti as right
left <- c(TEST_POPS$french, TEST_POPS$han)
right <- c(TEST_POPS$yoruba, TEST_POPS$mbuti)

qpwave_result <- qpwave(geno, left = left, right = right)

qpwave_results$test1 <- list(
  left = left,
  right = right,
  rank_tests = lapply(1:nrow(qpwave_result$rankdrop), function(i) {
    list(
      rank = qpwave_result$rankdrop$rank[i],
      chi_sq = qpwave_result$rankdrop$chisq[i],
      df = qpwave_result$rankdrop$df[i],
      p_value = qpwave_result$rankdrop$p[i]
    )
  })
)

cat("qpWave rank tests:\n")
print(qpwave_result$rankdrop)

# Save qpWave reference
write_json(qpwave_results, file.path(OUTPUT_DIR, "qpwave_reference.json"), pretty = TRUE)

# ============================================================================
# Test 5: qpAdm
# ============================================================================
cat("\n=== Generating qpAdm reference data ===\n")

qpadm_results <- list()

# qpAdm: French as target, Han and Yoruba as sources, Mbuti as outgroup
target <- TEST_POPS$french
left <- c(TEST_POPS$han, TEST_POPS$yoruba)
right <- c(TEST_POPS$mbuti)

# Note: qpAdm needs more outgroups for a proper model
# This is just a simple test case
qpadm_result <- tryCatch({
  qpadm(geno, target = target, left = left, right = right)
}, error = function(e) {
  cat("qpAdm error: ", e$message, "\n")
  NULL
})

if (!is.null(qpadm_result)) {
  qpadm_results$test1 <- list(
    target = target,
    sources = left,
    outgroups = right,
    weights = qpadm_result$weights$weight,
    weight_se = qpadm_result$weights$se,
    p_value = qpadm_result$p
  )

  cat("qpAdm weights:\n")
  print(qpadm_result$weights)
  cat("p-value: ", qpadm_result$p, "\n")
}

# Save qpAdm reference
write_json(qpadm_results, file.path(OUTPUT_DIR, "qpadm_reference.json"), pretty = TRUE)

cat("\n=== Golden master data generation complete ===\n")
cat("Reference files saved to: ", OUTPUT_DIR, "\n")
