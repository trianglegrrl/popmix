# Plan to finish popmix

Every feature is done the way the existing ones were:

1. **Reference first.** A generator script runs the reference tool on real
   data (AADR v62; the 54-sample DATES data for small-file cases) and stores
   every intermediate stage, not only the final numbers.
2. **Golden tests at the CLAUDE.md tolerances** (statistics rtol 1e-10, SEs
   1e-6, weights 1e-6, p-values 1e-4), written before the code.
3. Rust core, Python API, CLI command, nf-core module with nf-test, then
   docs (TESTING.md, DEVIATIONS.md).
4. A bug in the reference is fixed by default and reproduced under a compat
   flag (DEVIATIONS.md process).
5. One commit per feature.

Where a reference tool optimizes from random starts (ADMIXTURE, qpGraph),
the parts with no optimizer in them must match exactly, and the optimum must
be at least as good; matching the optimizer itself is a separate step.

## Phase A: new analyses

Ordered so that each builds on the last and the references that are already
installed come first.

### A1. D-statistic and outgroup-f3 scans, with heatmaps

- `qpdstat` over lists or all combinations: f4 mode from the f2 blocks
  (exists), and **D mode** (f4 / its normalizing denominator, per block),
  which needs a new per-block numerator/denominator from allele frequencies,
  as admixtools2 `qpdstat(f4mode = FALSE)` computes from genotypes.
- `qp3pop` scans; outgroup-f3 matrix f3(O; A, B) for all pairs.
- Output: tables (TSV/JSON) and heatmaps (f3 matrix; D/z for a fixed
  pair), with clustering order. Plotting is optional (matplotlib extra).
- **Reference:** admixtools2 `qpdstat()` (both modes) and `qp3pop()`.
- **Tolerance:** as f4/f3; Z and p as SE and p.

### A2. PCA with projection (smartpca)

- Genotype normalization as smartpca, eigen-decomposition, eigenvalues and
  Tracy-Widom statistics (twstats), outlier removal iterations, LD
  regression (`nsnpldregress`) if wanted, and **projection of ancient
  samples** (`lsqproject: YES`, `shrinkmode`).
- Randomized or streaming eigensolver for AADR-sized data, checked against
  the exact one.
- **Reference:** EIGENSOFT `smartpca` (installed, `/usr/bin/smartpca`).
- **Tolerance:** eigenvalues rtol 1e-6; eigenvectors up to sign, to
  smartpca's printed precision (it writes a fixed number of decimals, so
  the generator asks for the most it can print).

### A3. Kinship

- KING-robust kinship and IBS0 for diploid data.
- For ancient pseudo-haploid data, KING does not apply; the usual tool is
  pairwise mismatch normalized by the population (READ/READv2).
  Both are included (see Decisions).
- **Reference:** PLINK 2 `--make-king` / `--make-king-table` (installed;
  KING-robust), and KING itself if installed.
- **Tolerance:** kinship and IBS0 rtol 1e-10 (binary output).

### A4. Site-frequency spectra

- 1-D and joint SFS (folded and unfolded), from genotypes or allele counts;
  projection down to a smaller sample size (hypergeometric, as dadi/easySFS),
  needed for missing data and pseudo-haploid samples.
- Summary statistics from the SFS (π, Watterson's θ, Tajima's D).
- **Reference:** admixtools2 `joint_sfs()`, PLINK 2 `--freq counts` (exact
  counts), and dadi for projection and summaries.
- **Tolerance:** counts exact; projected spectra and statistics rtol 1e-10.

### A5. ADMIXTURE

- Maximum-likelihood ancestry proportions (Q) and allele frequencies (P):
  block relaxation with quasi-Newton acceleration, as ADMIXTURE; supervised
  mode; projection mode (fixed P); cross-validation error; bootstrap SEs.
- **Reference:** ADMIXTURE 1.3.0 (installed, `/usr/local/bin/admixture`).
- **Tolerance:** log-likelihood exact for given P and Q; projection mode
  (fixed P, a convex problem) Q to 1e-6; unsupervised runs: log-likelihood
  at convergence at least ADMIXTURE's, Q matched up to cluster order;
  cross-validation error to 1e-6 for the same folds.

### A6. LD-based admixture tools: ALDER and Rolloff

- ALDER: weighted LD curves (one or two references), the admixture test
  with its pre-test for shared LD, amplitude and date. Rolloff: the older
  correlation-based curve. Both share most of DATES' machinery.
- **Reference:** neither is installed. ALDER 1.03 builds from source (needs
  GSL and FFTW); Rolloff ships with older ADMIXTOOLS releases. Both are
  built locally as references.
- **Tolerance:** as DATES (curve values, fit, jackknife SE).

## Phase B: admixtools2 parity

1. **qpGraph weight optimizer.** R's Mersenne Twister with `set.seed`
   scrambling for the starts, and L-BFGS-B written from the published
   algorithm (Byrd, Lu, Nocedal and Zhu 1995) with R's finite-difference
   gradient (`ndeps = 1e-3`) and `fnscale`. Goal: the 3 strict xfails in
   `test_find_graphs_admixtools2.py` pass. If the from-scratch version
   cannot match R's exactly, record by how much.
2. **qpadm_rotate, qpwave_pairs**: reference tests; make `qpadm_multi`
   public; add `qpadm_models` and `qpadm_p`.
3. **Resampling**: `qpgraph_resample_snps/inds/multi`, `compare_fits`,
   `qpadm_resample_*`, `qpwave_resample_*`, `qpdstat_resample_*`.
4. **Graph search parity**: admixtools2's evolutionary `find_graphs` and its
   constraints (event order, zero / nonzero f4, number of admixtures).
5. **The rest**: `write_f2` / `extract_f2_subset`, DOT output, parsers for
   classic ADMIXTOOLS output, `qpfstats`, `qpf4diff`, `test_cladality`,
   `lazadm`, msprime simulation helpers.
6. **VCF and PED input**: f2 from them checked against admixtools2.

## Phase C: finish

1. Plots (graphs, fits, PCA, ADMIXTURE bar plots, heatmaps).
2. Cleanup: clippy warnings, TESTING.md counts, README.
3. bioconda recipe; containers.
4. nf-core: modules for every new command, test data uploaded to
   nf-core/test-datasets (provenance checked), then the PR.

## Decisions (2026-10-24)

- A3: both KING-robust (PLINK 2 reference) and READ-style kinship for
  ancient DNA (READv2 installed as its reference).
- A4: dadi installed in the venv as the projection reference.
- A6: ALDER 1.03 and an ADMIXTOOLS release with Rolloff downloaded and built
  as references.
- C1: matplotlib is an optional extra (`popmix[plot]`); tables are always
  written.
