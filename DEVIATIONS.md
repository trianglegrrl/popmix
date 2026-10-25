# Deviations from Reference Implementations

popmix is checked against reference tools (admixtools2, DATES) on real data.
Sometimes the reference has a bug or an artefact that changes results. popmix
then corrects it by default, but it can still reproduce the reference exactly,
so parity stays provable. This file lists every such correction and sets the
rules for adding one.

## Fixed in popmix

| ID | Reference bug | popmix default | Reproduce with |
|----|---------------|----------------|----------------|
| DATES-0001 | Mesh FFT lags wrap around | Lags with no real pairs are zero | `compat="dates-v753"` |
| DATES-0002 | Last FFT lag left unscaled | Every lag scaled | `compat="dates-v753"` |
| DATES-0003 | Fits the printed (rounded) values | Full precision | `compat="dates-v753"` |
| ADMIXTOOLS2-0001 | qpAdm `popdrop` nested path compares siblings, not nested models | A real nested path from the full model | `compat="admixtools2-2.0.10"` |
| ADMIXTOOLS2-0002 | Constrained qpAdm with ≥ 3 sources stops with an error | Returns the full result | nothing to reproduce |
| ADMIXTOOLS2-0003 | `joint_sfs()` drops SNPs whose allele count is not exactly an integer in floating point | Counts them | `compat="admixtools2-2.0.10"` (`joint_sfs`) |
| ADMIXTOOLS2-0004 | `write_dot()` renames nodes on edges (removes "." and "-") but declares them unrenamed | Quoted names, unchanged | `compat="admixtools2-2.0.10"` (`write_dot`) |
| ADMIXTOOLS2-0005 | `parse_qpgraph_graphfile()` returns node labels instead of population names, and fails by default | Population names | `compat="admixtools2-2.0.10"` (`parse_qpgraph_graphfile`) |
| ADMIXTOOLS2-0006 | PLINK files with chromosome codes X and Y stop with an error (`auto_only = TRUE`) | Reads them; X and Y are not autosomes | nothing to reproduce |
| ADMIXTOOLS2-0007 | A sample with no heterozygous call among the first 1,000 SNPs is pseudo-haploid, which catches diploid genomes | Pseudo-haploid only with no heterozygous call at any SNP | `pseudohaploid="admixtools2"` (`--pseudohaploid admixtools2`) |
| EIGENSOFT-0001 | smartpca reports the Tracy-Widom density as the p-value beyond its table | The tail | `compat="eigensoft-7.2.1"` |
| PLINK2-0001 | `--make-king-table` labels its two het-hom counts the wrong way round | Counts named as they are | nothing to reproduce (labels only) |
| ALDER-0001 | The inter-chromosome affine term counts pairs with SNPs that have no weight | Pairs it sums over | `compat="alder-1.03"` |
| READV2-0001 | Each jackknife window also counts the last SNP of the window before | Disjoint windows | `compat="readv2-2.01"` |
| READV2-0002 | The block jackknife never leaves out the last window | Every window left out once | `compat="readv2-2.01"` |
| READV2-0003 | The "kinship coefficient" 95% interval is that of the normalized P0 | The kinship coefficient's interval | `compat="readv2-2.01"` |

## The process

1. **Parity first.** popmix must reproduce the reference exactly before
   anything is corrected. The golden tests run in compat mode at the
   tolerances in CLAUDE.md. A correction is never a way to explain away a
   mismatch that isn't understood: find the cause in the reference source,
   reproduce it, and only then decide whether to correct it.

2. **One flag per deviation.** Each tool has a compat struct with one boolean
   per quirk (`popmix_core::dates::DatesCompat`,
   `popmix_core::Admixtools2Compat`).
   - `Default` / `CORRECTED` corrects everything; this is what users get.
   - A constant named after the reference version (`DatesCompat::V753`)
     turns every flag on.
   - Python and the CLI take a single `compat="<tool>-<version>"` /
     `--compat <tool>-<version>` switch that selects that constant.

   Individual flags exist so tests can switch one quirk at a time.

3. **Four kinds of test for each deviation:**

   | Test | Shows | Example (DATES-0001) |
   |------|-------|----------------------|
   | Golden | compat mode still matches the reference output | `test_dates_v753.py` `ld_table`, `full_fit`, ... |
   | Isolation | switching only this flag changes only what it should, and everything else still equals the reference | `test_circular_fft_lags_flag_only_changes_wrapped_bins`; live `default_ld_table` |
   | Oracle | the corrected behaviour is right, judged by something independent of the quirk: another mode, brute force, an analytic case, or simulation truth | `test_corrected_mesh_matches_exact_past_the_fft_window`; live `default_vs_exact_mode` |
   | Live numbers | how much it changes real results, recorded below | table under DATES-0001 |

   Write these tests first (RED), like any other change. Isolation tests
   should predict exactly which outputs change, not just that something does.

4. **A registry entry here**, giving the ID, what the reference does and
   where in its source, what popmix does by default, its impact, and its
   tests. IDs are `<TOOL>-NNNN` and are never reused. Cite the ID in code
   comments and test docstrings.

5. **Decide matched behaviours too.** When a reference behaviour looks odd
   but popmix deliberately matches it by default, record it under "Matched,
   not corrected" with the reasoning, so it isn't rediscovered later.

## DATES v753

`DatesCompat` in `crates/popmix-core/src/dates/compat.rs`;
`compat="dates-v753"` in Python, `--compat dates-v753` on the CLI.

### DATES-0001: mesh FFT lags wrap around

- **DATES:** In mesh mode (`qbin > 0`) the per-chromosome cross-correlation
  is a circular FFT of length n = 2^(⌈log₂ m⌉+1) for a window of m mesh
  cells (`fftsubs.c`). Lags beyond n − m also pick up the pairs at the
  negative lag n − d: short-range pairs counted again at long distance.
- **popmix default:** Lags at or beyond m have no real pairs and are zero.
  Because n ≥ 2m, no lag below m can wrap.
- **Flag:** `circular_fft_lags`.
- **Affects:** mesh mode only, when `maxdis` is longer than n − m for some
  chromosome. It changes a tail of the LD table, from the first wrapped lag to
  `maxdis`. Exact mode (`qbin = 0`, the default) is unaffected.
- **Live numbers** (full-data fit; jackknife estimate ± SE):

  | Config | DATES exact mode (oracle) | DATES mesh = compat | popmix corrected mesh |
  |--------|---------------------------|---------------------|-----------------------|
  | chr20_22_q10 (tail from 91.7 cM) | 3.165; 1.217 ± 3.648 | 2.890; 0.751 ± 3.411 | 3.160; 1.221 ± 3.645 |
  | simulation_q10 (tail from 91.7 cM) | 36.686; 35.680 ± 6.333 | 36.747; 35.738 ± 6.346 | 36.750; 35.750 ± 6.345 |
  | aadr_french_q10 (maxdis 0.5) | not run | 140.595 ± 25.037 | unchanged, no lag wraps |

  In the chr20–22 tail, DATES's mesh gets pair counts up to 113% off from
  exact mode. The corrected mesh is at most 5% off, close to the mesh's
  normal binning error (2.9% elsewhere).
- **Tests:**
  - `ld_decay.rs`: `test_corrected_mesh_matches_exact_past_the_fft_window`
    (oracle), `test_v753_mesh_reproduces_dates_quirks`,
    `test_circular_fft_lags_flag_only_changes_wrapped_bins` (isolation).
  - `tests/golden/test_dates_v753.py`: `default_ld_table` (isolation, live),
    `default_vs_exact_mode` (oracle, live).
  - `tests/python/test_dates.py::test_compat_changes_only_mesh_tail`.

### DATES-0002: last FFT lag left unscaled

- **DATES:** `fftauto` divides its output by n for lags 0 … outmax − 1 but
  not for lag outmax = min(n − 1, diffmax), which comes out n times too large.
- **popmix default:** Every lag is scaled.
- **Flag:** `unscaled_last_fft_lag`.
- **Affects:** Only one bin, and only when DATES-0001 is also in effect.
  Lag outmax holds real pairs only when wrapped pairs land there. Otherwise
  it is lag diffmax, which falls in an unreported bin. The live numbers are
  part of DATES-0001's.
- **Tests:** `ld_decay.rs`: `test_unscaled_last_lag_flag_only_changes_its_bin`
  (isolation: exactly the predicted bin changes).

### DATES-0003: fitting printed values

- **DATES:** `dates` writes the LD table as text (correlation `%12.6f`), and
  `dates_expfit` fits what it reads back. The jackknife likewise reads the
  per-chromosome estimates from `.jin` at `%9.3f`.
- **popmix default:** Fits and runs the jackknife at full precision.
- **Flag:** `fit_printed_values`.
- **Affects:** The estimates, by rounding noise. On AADR French.DG the
  full-data fit moves by 0.002 generations.
- **Tests:** The golden `full_fit`, `jackknife_by_chrom` and
  `jackknife_result` checks pass only with it on (compat). `default_result`
  checks that the full-precision result stays within DATES's reported
  precision wherever no other correction applies.

### Matched, not corrected

- **`mincount` is ignored.** DATES reads `mincount` but strips it ("mincount
  stripped" in `dates.c`): a SNP is used once each source has one called
  sample. popmix defaults `min_source_count` to 1 to match. Requiring more
  samples would be a defensible change. It would be a new deviation with its
  own flag, needing an oracle (e.g. simulation truth) to show it helps.
- **Mesh drops lag 0.** Pairs within one mesh cell (under binsize/qbin
  apart) are never counted, because lag 0 also holds each SNP paired with
  itself. This affects only the first bin, which is below the fit range.
- **Uninitialized `numbins`.** This crashes DATES's exact mode at random and
  doesn't change results. It is fixed in the reference build
  (`tests/golden/dates/dates_numbins_fix.patch`), not reproduced.

## admixtools2

`Admixtools2Compat` in `crates/popmix-core/src/qpadm/compat.rs`;
`compat="admixtools2-2.0.10"` in Python (`qpadm`). The CLI's qpadm output
has no popdrop table, so it has no `--compat` switch yet.

Apart from the entries below, popmix follows admixtools2 2.0.10 defaults:
`maxmiss = 0`, `poly_only = "f2"`, blocks built after filtering, blocks
weighted by length in the jackknife, FST on its own SNP set and blocks, and
qpWave/qpAdm fits by the C++ code path (`cpp = TRUE`, the default).

### ADMIXTOOLS2-0001: popdrop's nested path compares siblings

- **admixtools2:** `drop_pops` (`R/qpadm.R`) fits every subset of the
  sources, then is meant to walk a nested path down from the full model: at
  each rank, among the children of the path so far (one more source
  dropped), take the feasible model with the lowest χ² (`top_n(1,
  -jitter(chisq))`, the jitter only breaking ties), and test each step
  against the one before. Two mistakes stop that:
  - the path starts with *every* model at full rank − 1, not the full
    model;
  - `child_patterns(nested$pat)` reads only `str_locate_all(...)[[1]]`, the
    children of the first pattern, which are at full rank − 1 again. So
    nothing below that rank is ever added.

  With two or more sources, `best` marks exactly the models that drop one
  source. `dofdiff`, `chisqdiff` and `p_nested` compare those siblings with
  each other: dofdiff is always 0, so p_nested is always 0 or 1.
- **The intent is clear from both admixtools.** The R loop and its
  `top_n(1, -chisq)` over feasible children say what was meant. Classic
  ADMIXTOOLS (`qpAdm.c`, "best pat") does the same thing per level: the
  feasible pattern with the best p-value, compared with the level above on
  1 dof and flagged "not nested" when it isn't. It has its own bug: the
  `nfeas[wt] == 1` branch resets the best on every pattern until a second
  feasible one turns up, so an infeasible pattern with p = 0 can become the
  incumbent. Run on Druze = Sardinian + Dinka + Han (8.0.2, AADR v62), it
  picks 100 (χ² 790) over 001 (χ² 12.7) at level 1, then reports level 2 as
  "not nested". popmix doesn't target classic ADMIXTOOLS; this is evidence
  of the intent only.
- **popmix default:** The path starts at the full model. Each step drops
  one more source, taking the feasible child with the lowest χ² (the first
  in `power_set` order on an exact tie, in place of `jitter`). It goes down
  to a single source, or stops early if no child is feasible. Each step's
  row compares it with the model it came from: `dofdiff` = 1, `chisqdiff` =
  the χ² increase, and `p_nested` its upper tail. This is a likelihood-ratio
  test of whether the source just dropped is needed. Requiring children
  makes every comparison truly nested, which classic ADMIXTOOLS's per-level
  choice doesn't guarantee.
- **Flag:** `popdrop_sibling_path`.
- **Affects:** Only `best`, `dofdiff`, `chisqdiff` and `p_nested`. Every
  model fit (χ², p, weights, feasibility) is unchanged. With one source
  there is no path: admixtools2 marks nothing, popmix marks the full model.
- **Live numbers** (`best` models; χ² difference / dof, p):

  | Model | admixtools2 2.0.10 | popmix corrected |
  |-------|--------------------|------------------|
  | Druze = Sardinian + Dinka + Han | 001 vs 010: −3.54/0, p 1; 010 vs 100: −792.85/0, p 1 | 000 → 001 (drop Han): 6.30/1, p 0.012 → 011 (Sardinian alone): 85.99/1, p 1.8e-20 |
  | French = Han + Yoruba | 01 vs 10: −7418.82/0, p 1 | 00 → 01: 613.82/1, p 1.7e-135 |
  | Papuan = Han + French | 01 vs 10: −70.35/0, p 1 | 00 → 01: 296.53/1, p 1.9e-66 |
  | French = Sardinian + Yamnaya | 01 vs 10: −0.50/0, p 1 | 00 → 01: 0.18/1, p 0.67 |
  | same, `maxmiss = 0.25` | 01 vs 10: 0.06/0, **p 0** | 00 → 10: 0.00/1, p 0.98 |

  In the last row, admixtools2 calls a χ² difference of 0.06 significant
  (p = 0, because dofdiff is 0). The constrained models give the same paths
  as their unconstrained twins.
- **Tests:**
  - `fit.rs`: `test_corrected_popdrop_path_is_nested_and_greedy` (oracle,
    hand-built case: skips an infeasible best child and a better non-child),
    `test_corrected_popdrop_path_stops_without_a_feasible_child`,
    `test_v2_0_10_popdrop_compares_the_one_drop_siblings` (compat),
    `test_popdrop_compat_flag_only_changes_the_nested_columns` (isolation).
  - `tests/golden/test_fstats_admixtools2.py`: `qpadm_popdrop` (golden,
    compat mode); `qpadm_popdrop_corrected` (live isolation: only the nested
    columns differ from compat; live oracle: they follow the nested path
    walked independently over admixtools2's own model fits).
  - `tests/python/test_qpadm_compat.py` (the `compat` argument).

### ADMIXTOOLS2-0002: constrained qpAdm stops with three or more sources

- **admixtools2:** `qpadm(constrained = TRUE)` errors with "matrix D in
  quadratic function is not positive definite". `drop_ranks` refits every
  lower rank through `qpadm_fit` → `cpp_qpadm_weights`, which always
  computes weights. At rank r < sources − 1, x = [A 1]′ has only r + 1 rows,
  so D = x′x (sources × sources) is singular. The unconstrained `solve`
  tolerates that; quadprog doesn't. Those weights are then thrown away
  (`addweights = FALSE`): rankdrop's f4rank, dof, χ² and p come from the
  rank fit alone.
- **popmix default:** The rank tests never compute weights. Constrained
  weights are solved only at full rank, for the model and for every popdrop
  subset, where the problem is well posed. So popmix returns weights, SEs,
  rankdrop and popdrop.
- **Flag:** none. There is no admixtools2 output to reproduce.
- **Verification:** the `dg_11pop` reference for Druze = Sardinian + Dinka +
  Han (constrained) is built from the admixtools2 internals `qpadm()` uses:
  - `cpp_qpadm_weights` and `cpp_get_weights_covariance` for the weights
    and SEs;
  - `drop_pops(constrained = TRUE)` for popdrop, which runs fine;
  - `drop_ranks` run unconstrained for rankdrop, which is identical by the
    argument above.

  popmix matches all four: weights and SEs to 1e-11 relative, the model χ²
  to 2e-12, the popdrop χ² to 4e-11. Constrained
  weights 0.9713 / 0.0287 / 0 against unconstrained 1.0087 / 0.0179 /
  −0.0266.
- **Tests:** golden `qpadm_weights`, `qpadm_fit`, `qpadm_popdrop`,
  `qpadm_popdrop_corrected` on `dg_11pop`; `fit.rs`
  `test_constrained_popdrop_with_three_sources`.

### ADMIXTOOLS2-0003: joint_sfs() drops SNPs whose count is not exactly an integer

- **admixtools2:** `joint_sfs()` (`R/fstats.R`) takes the allele count of
  each population as `afs * counts`, frequency times sample size, groups the
  SNPs by those values and joins them onto the integer grid
  `expand_grid(0:max_count, ...)`. Frequency times sample size is not always
  an exact integer in floating point (23/42 * 42 = 23.000000000000004), and
  such a SNP matches no row of the grid, so it is silently left out.
- **Effect:** 930 of 21,162 SNPs (4.4%) of `family_3pop` and 67,979 of
  1,233,013 (5.5%) of `dg_3pop` are dropped, depending only on the sample
  sizes and counts; the spectrum is not just smaller but distorted, since
  which counts round badly depends on the count.
- **popmix:** counts within 1e-6 of an integer are rounded. Half counts (a
  heterozygous call in a pseudo-haploid sample) are still left out, and
  reported.
- **Also documented, not a bug in the numbers:** `joint_sfs(pref = ...)`
  calls `geno_to_afs()`, which does not exist in 2.0.10; the reference
  generator calls `packedancestrymap_to_afs()` and passes its result.

### ADMIXTOOLS2-0004: write_dot() renames nodes on edges but not where it declares them

- **admixtools2:** `write_dot()` (`R/io.R`) writes each edge with "." and "-"
  removed from both node names, because DOT does not allow them in an
  unquoted name, and declares the internal nodes (`[shape = point]`) under
  their unchanged names.
- **Effect:** population labels change (Han.DG is drawn as HanDG); two
  populations that differ only in those characters (A.B and AB) become one
  node; and an internal node with a "." in its name (n.1) is declared as
  something else (Graphviz reads `n.1` as two stray nodes, `n` and `.1`)
  while its edges use `n1`, which is drawn as a labelled ellipse instead of
  a point.
- **popmix:** quotes every name and keeps it.

### ADMIXTOOLS2-0005: parse_qpgraph_graphfile() keeps the labels, and fails by default

- **admixtools2:** a qpGraph graph file names each population's node with
  `label <population> <label>` and uses the labels in its edges.
  `parse_qpgraph_graphfile()` (`R/wrappers.R`) builds its renaming map the
  wrong way round (population to label), so it never replaces a label by a
  population, whatever `uselabels` says (`parse_qpgraph_output_edges()`
  reads qpGraph's output, which lists label first, correctly). With its
  default `split_multi = TRUE` it then stops with "names do not match
  previous names" in 2.0.10, even when nothing is multifurcating.
- **Effect:** the graph's leaves are the labels (Chimp, not Chimp.REF), so
  the graph cannot be fitted to f2 statistics of the populations it names;
  with the defaults, no graph file can be read.
- **popmix:** returns population names (labels with `uselabels=True`).
- **Reference:** generated with `split_multi = FALSE`.

### ADMIXTOOLS2-0006: PLINK chromosome codes X and Y stop f2_from_geno()

- **admixtools2:** with `auto_only = TRUE` (the default), `discard_snps()`
  (`R/io.R`) parses chromosome names as numbers and stops with "Could not
  parse chromosome numbers!" if any is not a number, so a .bim with PLINK's
  usual codes X and Y (as PLINK 2 and popmix's `convertf` write them)
  cannot be read; with numeric codes (23, 24) it works.
- **popmix:** reads X, Y, XY and MT as non-autosomes.
- **Reference:** generated from a copy of the PLINK files with numeric
  codes (`plink2 --output-chr 26`).

### ADMIXTOOLS2-0007: diploid genomes taken as pseudo-haploid

- **admixtools2:** a sample is pseudo-haploid if it has no heterozygous call
  among the first 1,000 SNPs (`cpp_packedancestrymap_ploidy`, and its PLINK
  counterpart), and its heterozygous calls then count half an allele. A
  diploid genome whose first calls happen to be homozygous is caught: in
  AADR v62.0, 9 of 17,629 samples, all high-coverage genomes, each with
  959 to 984 calls among those SNPs: HGDP00068.DG (Balochi.DG),
  HGDP00614.DG (BedouinA.DG), HGDP00612.DG (BedouinB.DG), HGDP00581.DG
  (Ignore_Druze.DG), NA18606.DG (CHB.DG), HG04056.DG (ITU.DG), and
  HG02691.DG, HG02699.DG and HG03237.DG (PJL.DG). No ancient sample is
  affected.
- **popmix:** a sample is pseudo-haploid only if it has no heterozygous call
  at any SNP (`pseudohaploid="detect"`, for `extract_f2`, `qpdstat`,
  `qp3pop`, `sample_ploidy`, `read_genotypes` and the CLI commands that read
  genotypes; `sfs` already did this). `pseudohaploid="admixtools2"` keeps the
  first-1,000 rule.
- **Impact:** f2 of a population holding such a sample. With NA18606.DG among
  CHB.DG's samples, f2(CHB.DG, Han.DG) is 0.000385318 by default and
  0.000384351 with the admixtools2 rule (0.25%, 0.05 SE); a population of
  that sample alone would change far more. Populations without it change
  only through the shared SNP set (`poly_only`), and in the AADR test not at
  all. No admixtools2 golden holds such a sample.
- **Tests:** isolation and oracle on small data
  (`tests/python/test_extract_f2_formats.py`: a diploid sample made
  homozygous over the first 1,000 SNPs is diploid by default, and only
  statistics with its population change), live on AADR
  (`tests/golden/test_ploidy_rule.py`: only NA18606.DG is reclassified, only
  CHB.DG's statistics change, and the default equals counting every sample
  diploid).

### Matched, not corrected

- **qpAdm weight SEs scale by sqrt(blocks - 1).** admixtools2's C++ path
  (`cpp_get_weights_covariance`, the default) and its R path
  (`get_weights_covariance`) scale the leave-one-out weights differently,
  by sqrt(n - 1) against (n - 1)/sqrt(n). popmix matches the default C++ path.
- **qpGraph admixture weights: the same optimizer, apart in rounding.**
  admixtools2 optimizes them with R's `optim(method = "L-BFGS-B")`
  (`fnscale` 1e-6, `maxit` 1e4, finite-difference gradients with step 1e-3)
  from `numstart` starts drawn by `runif` after `set.seed(seed)`. popmix
  does the same: R's Mersenne Twister with `set.seed`'s scrambling
  (`rcompat::RRng`, the same numbers as R for every seed) and L-BFGS-B
  written from its published description, with R's finite differences,
  scaling and stopping rules (`rcompat::optim_lbfgsb`). On bounded test
  problems it evaluates the objective the same number of times at the same
  points as R to 1e-8 and ends within 1e-10, with the same convergence code
  and message (`crates/popmix-core/tests/rcompat_golden.rs`). It is not
  bit-for-bit: R takes its dot products from its BLAS (OpenBLAS here, whose
  summation order depends on memory alignment), so the two differ in the
  last bit.

  Everything that follows from the weights is reproduced exactly: the f3
  statistics and their covariance, the edge lengths (the same quadratic
  program `qpgraph` passes to quadprog, solved exactly by an active-set
  method), the fitted f3 and the score, all within 3e-11 of admixtools2's at
  the same weights. On the qpGraph golden graphs every start (the same
  starts as admixtools2's) ends at admixtools2's score to 3e-12 and its
  weights to 3e-7, so the best fit and the low/high range over the starts
  match too.

  Over every graph admixtools2 `generate_all_graphs()` makes for 4 and 5
  `dg_6pop` populations with one or two admixtures (2871 graphs,
  `all_graphs.json`):

  | graphs | same score (1e-6) | popmix lower | admixtools2 lower |
  |---|---|---|---|
  | 4 populations, 1 admixture (84) | 84 | 0 | 0 |
  | 5 populations, 1 admixture (1485) | 1485 | 0 | 0 |
  | 4 populations, 2 admixtures (1302) | 1219 | 50 | 33 |

  The 83 two-admixture graphs that differ are ones whose score is nearly
  flat along a ridge of local minima; the score's own rounding (1e-12 to
  1e-11 apart from admixtools2's, from inverting the f3 covariance and the
  quadratic program) changes a finite-difference gradient at its eighth
  digit, and a path that starts identical to admixtools2's ends in a
  neighbouring minimum. The best graph is the same in every set (2.382203,
  0.00406303, 0.00114228). Before popmix used R's starts and L-BFGS-B, its
  own quasi-Newton search from its own starts agreed on 59, 1310 and 825 of
  these graphs.
- **The evolutionary graph search draws its own random numbers.**
  `find_graphs(method="evolutionary")` makes admixtools2's moves
  (negative-edge rearrangements, subtree prune and regraft, leaf swaps,
  admixture moves and flips, root placement, adding and removing admixture
  edges) on the best graph so far, with its selection, deduplication,
  constraints and stopping rules, but its random choices come from Python's
  generator, so a run takes a different path from admixtools2's with the
  same seed. Compared over runs (`test_graph_search_admixtools2.py`): with
  5 populations and 1 admixture, popmix reaches the optimum in at least half
  of 30 runs and admixtools2 in 5 of 8; with 4 populations and 2
  admixtures, popmix's best run reaches the exhaustive optimum (0.001137),
  better than any of admixtools2's 8 runs (0.001239), but its median run
  ends higher (0.001655 against 0.00142). Both fit a similar number of
  graphs (52 and 56, 67 and 80).
- **admixtools2 2.0.10 needs igraph 1.x for qpGraph.** It passes vertex
  sequences to `get.edge.ids()`, which igraph 2 rejects. The reference
  generator shims that one call (same vertex ids, as a plain vector) and
  records it; this affects only how the reference is produced.
- **joint_sfs() counts SNPs with missing data at their observed count.** A
  SNP where 40 of 42 alleles are observed enters the 42-allele spectrum as if
  its count were among 42. popmix's `joint_sfs` matches; `sfs` with a
  projection (dadi's method) is the spectrum for data with missing calls.
- **X and Y counted with autosomal ploidy.** admixtools2 reads every SNP and
  counts a diploid male's X as two alleles. popmix's `joint_sfs` does the
  same; `sfs` uses autosomes by default.
- **FST drops missing blocks per pair.** `fst(f2_dir, pop1, pop2)` reads one
  pair and leaves out that pair's undefined blocks; `fst(prefix)` over many
  populations leaves out a block for every pair if any pair is undefined
  there. popmix's `fst(a, b)` matches the former. The two agree whenever no
  block is undefined, as in every reference config.

- **f2 blocks with no data for some pair are dropped for every pair.** With
  `maxmiss > 0` a block can hold no SNP where two populations both have
  data; `f2_from_geno(remove_na = TRUE)` drops it for all pairs, and popmix
  does the same (`ancient_maxmiss` drops 2 of 713 blocks). Because of this no
  missing block reaches qpWave/qpAdm, so `est_to_loo_nafix`'s zero-filling
  never applies in this pipeline.

### Interface differences

These change no number that both tools compute.

- **qpdstat D mode on f2 blocks is an error.** admixtools2 ignores
  `f4mode = FALSE` when given f2 blocks and returns f4; D needs genotypes.
  popmix raises an error rather than return f4 under the name D.
- **Statistics that cannot be estimated report n = 0.** admixtools2 returns
  NA for every column, `n` included (for example f3 with a single
  pseudo-haploid target, whose bias correction is 0/0 at every SNP); popmix
  returns NaN for the statistic and 0 for `n`.
- **qp3pop on f2 blocks keeps the combination order.** admixtools2 returns
  it sorted by population (a side effect of dplyr `summarize`), but keeps the
  order from genotypes and in qpdstat; popmix keeps it in every case.
- **qpadm_p at intermediate ranks has no weights.** admixtools2 solves a
  singular system for the weights of a fit at a rank between 0 and n - 1
  (armadillo warns and returns an approximate solution) and reports
  `feasible` from them; popmix returns no weights there and `feasible`
  true, as admixtools2 does when it has none. The chi-squared, degrees of
  freedom and p-value at every rank agree.
- **Tables are data frames, populations are list columns.** `qpwave_pairs`,
  `qpadm_rotate`, `rotate_models`, `qpadm_multi(full_results = FALSE)`,
  `qpadm_models` and the resampling functions return pandas data frames
  with admixtools2's rows and columns; `qpadm_p` returns one row as a dict.
  Resampling block indices (`sel`) are 0-based.
- **PLINK genetic positions in centimorgans are recognized.** admixtools2
  (and EIGENSOFT, and the Reich lab's PLINK files) read and write the third
  column of a .bim or .map in Morgans; PLINK 2 documents centimorgans.
  popmix reads Morgans too, but treats a file whose largest position is
  over 10 (no chromosome is 10 Morgans long) as centimorgans, which
  admixtools2 would take 100 times too long; it writes Morgans.
- **qpAdm and qpWave from a genotype prefix use f2 blocks.** Given a
  genotype prefix, admixtools2's `qpadm()` and `qpwave()` compute each f4
  from the genotypes on the SNPs its four populations share; `popmix
  qpadm --input` (and `qpwave`) extract f2 blocks for the model's
  populations first, on the SNPs they all share, which is admixtools2's
  `qpadm(f2_from_geno(prefix, pops))`. popmix agrees with that exactly; the
  two admixtools2 paths differ (on 16 AADR v66.p1 populations, weights
  0.8206 and 0.8130).
- **Leaving out individuals recomputes the usual f2.** admixtools2's
  `*_resample_inds` read per-individual and per-pair counts precomputed by
  `extract_counts()` and combine them into f2 (`indpairs_to_f2blocks`), an
  estimator of its own. popmix's take a genotype prefix and recompute f2
  without each individual as `extract_f2` (admixtools2 `f2_from_geno(inds,
  pops)`) does, so they agree with admixtools2's statistics on those f2
  blocks.
- **msprime scripts: which parent gets `admix_default`.** admixtools2 gives
  it to the parent whose name sorts first in R's locale (`min()`); popmix
  sorts by code point, which agrees for names of one case. With the default
  of 0.5 it makes no difference.
- **Resampling seeds are explicit.** admixtools2's bootstrap draws from R's
  global random state, and runs qpGraph in parallel streams (furrr) unless
  given a `seed`. popmix takes the bootstrap seed as `boot_seed` (the same
  draws as admixtools2 after `set.seed(boot_seed)`) and passes `seed` to
  qpgraph.

## EIGENSOFT 7.2.1 (smartpca)

`EigensoftCompat` in `crates/popmix-core/src/pca/tracy_widom.rs`;
`compat="eigensoft-7.2.1"` in Python (`smartpca`), `--compat
eigensoft-7.2.1` on the CLI (`pca`).

### EIGENSOFT-0001: Tracy-Widom density reported as the p-value beyond the table

- **smartpca:** `gettw` (`src/nicksrc/statsubs.c`) interpolates the
  Tracy-Widom tail and density from a table on [-8, 8]. Beyond the table it
  uses the Margetis-Edelman asymptotic forms, but assigns them the wrong way
  round: `*tailp = twdensx (x); *densp = twtailx (x);`. Every eigenvalue
  with a TW statistic above 8, which is any strongly significant one, gets
  the asymptotic density as its p-value.
- **Effect:** the density is the tail times sqrt(x), so the p-value is too
  large by sqrt(TW): a factor of 5.1 for the first eigenvalue of
  `family_all` (TW 26.3; 7.0e-41 printed, 1.3e-41 correct). Such p-values
  are far below any threshold, so no conclusion changes.
- **popmix:** the asymptotic tail. `test_corrected_p_values_differ_only_beyond_the_table`
  checks that only these p-values change, by exactly 1/sqrt(TW).

### Matched, not corrected

- **Pseudo-haploid samples count as diploid.** smartpca reads genotypes 0, 1
  and 2 whatever the sample's ploidy, so a pseudo-haploid sample coded 0/2
  enters the mean and the normalization as a homozygote. popmix matches.
- **Eigenvector signs.** They are whatever the eigensolver returns (LAPACK
  `dspev` in smartpca, nalgebra in popmix), so they can differ between the
  two; the tests compare each eigenvector up to sign.

## PLINK 2 (KING-robust kinship)

Checked against PLINK v2.0.0-a.7LM (31 Dec 2025).

### PLINK2-0001: het-hom counts labeled the wrong way round

- **PLINK 2:** `--make-king-table counts cols=...,ibs1` reports
  `HET1_HOM2` and `HET2_HOM1`. `HET1_HOM2` holds the SNPs where IID2 is
  heterozygous and IID1 homozygous, the reverse of its name: for the first
  pair of `family` (IID1 NA00001, IID2 NA00000) it is 4210, and NA00000 is
  the heterozygote at 4210 of those SNPs, NA00001 at 4360.
- **Effect:** none on the kinship, which uses only their sum.
- **popmix:** `het1_hom2` is het in `id1`, hom in `id2`. The golden test
  compares it with PLINK's `HET2_HOM1`.

### Matched, not corrected

- **No heterozygote, kinship minus infinity.** KING-robust divides by the
  smaller heterozygote count; a sample with none (any pseudo-haploid sample)
  gives minus infinity, as in PLINK 2. popmix matches; READ-style kinship is
  the estimator for such samples.

## READv2 (READ-style kinship)

`ReadCompat` in `crates/popmix-core/src/kinship/read.rs`;
`compat="readv2-2.01"` in Python (`read_kinship`), `--compat readv2-2.01`
on the CLI (`kinship --method read`). Checked against READv2 2.01
(GuntherLab/READv2 commit 73739be).

### READV2-0001: jackknife windows overlap by one SNP

- **READv2:** windows are cut where `floor(position / window_size)` changes,
  at indices `i` (the last SNP of a window), and taken as
  `arr[pre_ind:i+1]` with `pre_ind = i` afterwards, so every window but the
  first starts with the previous window's last SNP. The same happens to the
  20 Mb windows that tell parent-offspring from siblings.
- **Effect:** about one SNP in 2000 counts twice in the window SEs and the
  first-degree window shares; the genome-wide P0 and the kinship are not
  affected.
- **popmix:** disjoint windows.

### READV2-0002: the jackknife leaves out each window but the last

- **READv2:** `blockJN_SE` computes the leave-one-out means for
  `range(num_blocks - 1)`, so the last window is never left out, but scales
  the variance by `(num_blocks - 1) / num_blocks` as if all were.
- **Effect:** the SE misses one of about 600 replicates. With READV2-0001,
  correcting both changes the SE by 1.6% for the median `aadr_gurgy` pair
  and by up to 28%; no relationship or first-degree type changes there.
- **popmix:** every window left out once.

### READV2-0003: the kinship interval is the normalized P0's

- **READv2:** `KinshipCoefficient_95percCI[min,max]` is the normalized P0
  plus and minus 1.96 SE, next to `KinshipCoefficient`, which is 1 minus the
  normalized P0. For a first-degree pair with kinship 0.272 it reports
  [0.714, 0.741].
- **popmix:** the interval of the kinship coefficient, 1 minus that interval
  reversed ([0.259, 0.286]).

### Matched, not corrected

- **Unweighted windows in the jackknife.** The SE is the jackknife of the
  plain mean of the 5 Mb windows' P0, not of the genome-wide P0 (which
  weights windows by their SNPs); popmix matches.
- **Unrelated pairs measured against the second-degree threshold.** For a
  pair classified unrelated, `Zdown` is its distance in SEs to the
  second-degree boundary (0.90625), not to the third-degree range (0.953125),
  presumably because third-degree pairs cannot be told apart reliably;
  popmix matches.
- **Windows by base-pair bin only.** READv2 cuts windows where the bin
  changes and ignores chromosome boundaries; on human data the bin changes
  at every chromosome boundary anyway. popmix also cuts at chromosome
  boundaries, except in compat mode.

## ADMIXTURE 1.3.0

Checked against ADMIXTURE 1.3.0 (closed source; the behaviour below is
inferred from its output).

### Matched, not corrected

- **The likelihood.** ADMIXTURE's log-likelihood has no binomial constant
  (heterozygotes are not counted twice), and P is the frequency of the
  allele the genotype files do not count. popmix matches both: at
  ADMIXTURE's printed Q and P its log-likelihood is ADMIXTURE's to the
  precision lost by printing them.
- **Bounds.** Q and P stay within [1e-5, 1 - 1e-5], as ADMIXTURE's, so
  supervised reference samples have Q 0.99999, not 1.
- **Cross-validation error.** It is half the mean deviance of the held-out
  genotypes, as ADMIXTURE reports it. ADMIXTURE's folds come from its own
  random number generator, which popmix cannot reproduce, so the errors agree
  only as closely as two random splits do: within 0.9% and 1.2% on the small
  data (K = 2, 3), 0.08% on `aadr_1kg_k3`.
- **Supervised cluster order.** Reference populations are numbered in the
  order their first sample appears, as ADMIXTURE numbers them.

### Differences in the optimizer, not deviations

popmix maximizes the same likelihood by block relaxation (a Newton step on
each sample's Q, then on each SNP's P, as ADMIXTURE) accelerated with
SQUAREM, not ADMIXTURE's quasi-Newton acceleration, and from its own random
start. Where the maximum is well defined the two agree (`family_k2`,
`aadr_1kg_k3`: Q to 1e-4). Where it is not, they can stop at different
local maxima: at K = 3 on `family` (two source populations and admixed
samples), ADMIXTURE's run ends at -1044310.79 and popmix's at -1044149.28.
With P fixed (projection) the problem is convex and the two agree to the
printed precision.

## ALDER 1.03

`AlderCompat` in `crates/popmix-core/src/alder/mod.rs`; `compat="alder-1.03"`
in Python (`alder`), `--compat alder-1.03` on the CLI.

### ALDER-0001: the affine term counts pairs it does not sum over

- **ALDER:** the inter-chromosome affine term (`compute_inter_chrom_affine`)
  is a sum over pairs of SNPs on different chromosomes divided by their
  number. The sum skips SNPs without a weight (a reference with no data
  there, `snp_ignore`), but the count (`AffineData.count`, the SNPs with no
  missing data in the test population) includes them.
- **Effect:** none when every SNP has a weight (`family`, `aadr_asw`). With
  a pseudo-haploid ancient reference (`aadr_asw_lbk`, ASW with YRI and
  Germany_EN_LBK) the affine term is 2.2% too small and the date 0.6% too
  young (5.126 against 5.156 generations).
- **popmix:** counts only the SNPs it sums over.

### Matched, not corrected

- **amp_tot.** The total amplitude is amp_exp plus half the affine term, as
  ALDER reports it.
- **Genotypes as called.** Pseudo-haploid samples count as homozygous
  diploids, as ALDER reads them.
- **Printed counts.** ALDER's raw output truncates its (FFT-computed) bin
  counts, and prints them with `%Ld`, a `long double` format, given a
  `long long`, so the inter-chromosome pair count comes out wrapped to 32
  bits. This affects only its output; the golden tests allow for it.

## ROLLOFF (ADMIXTOOLS 7.0.2)

Checked against rolloff from ADMIXTOOLS 7.0.2, with its curves fit by
rexpfit.r's method and its weighted jackknife by dowtjack, as expfit.sh.

### Matched, not corrected

- **Bins labeled by their upper edge.** A pair at distance d goes in bin
  floor(d / binsize), written (and fit) at (bin + 1) binsize.
- **Run mode 0 weights.** (pA - pB) / sqrt(mean frequency), as rolloff
  (the source marks it "??"); run mode 1, pA - pB, is the default and the
  one tested.
- **No weight centering.** The source's comments say weights are forced to
  mean 0, and `clipwt` was meant to zero tiny weights, but neither does
  anything; popmix matches.
- **The fit is unbounded.** rexpfit.r starts nls from a DEoptim search with
  the decay in (0, 500] but nls itself is unbounded; on noisy curves
  (the 4 admixed samples of `family`) it leaves the range (decay 998).
  popmix's least-squares fit keeps the decay within (0, 500], so the two
  agree where the curve determines a date (`aadr_asw`). There nls stops at
  its convergence tolerance on a flat ridge of the affine term and the
  amplitude: popmix's fits have a smaller squared error in every case (all
  data and each chromosome left out), and their decays agree within 7.7e-6,
  the jackknife within 3.4e-6.
- **expfit.sh's jackknife assumes chromosomes 1-22.** wtjack.pl matches
  per-chromosome fits to chromosomes 1 to 22 in order; with a subset of
  chromosomes its weights and estimates are misaligned. popmix weights each
  chromosome present by its SNPs, as wtjack.pl intends, and the reference
  is generated the same way.


## mergeit (AdmixTools 8.0.2)

`popmix merge` is checked against mergeit (tests/golden/formats/
generate_merge.sh); `plink --bmerge` stops at the first SNP whose alleles
need reconciling and so serves only to show the SNPs PLINK would ask to
`--flip`.

### Matched, not corrected

- **Strands.** A SNP whose alleles in file B are the complement of file A's
  (A/G against T/C or C/T) is the same SNP on the other strand; it is kept
  and its genotypes counted against file A's alleles.
- **A/T and C/G SNPs are left out** (mergeit's `strandcheck`), because their
  strand cannot be told from the alleles. So is a SNP with a third allele.
- **Monomorphic SNPs.** The unknown allele (PLINK `0`, EIGENSTRAT `X`, VCF
  `.`) is taken from the other file.

### Interface differences

- **SNPs are paired by position, not name.** mergeit and PLINK pair SNPs by
  name, so the same SNP under two names (rs ID and chr:pos) is not merged.
  popmix pairs by chromosome and position; where several SNPs of file B
  share a position it takes the one with file A's name, and with none, it
  leaves the SNP out and counts it.
- **Every SNP left out is counted and printed** by reason.

## convertf (AdmixTools 8.0.2): reading .snp and .bim

Checked black-box with AdmixTools 8 convertf.

### Matched, not corrected

- **Six columns.** A `.snp` line needs its two alleles (convertf: "6
  columns required").
- **One character per allele.** An allele written with more characters
  (an indel, `AT`) is read as its first, as convertf reads it; popmix's
  formats hold single-base alleles. VCF indels are skipped instead (A7).

### Interface differences

- **A genetic position must be a number.** convertf reads `abc` as 0 (C
  `atof`), which puts the SNP at the start of its chromosome's map; popmix
  stops with the file and line. A map of all zeros is still read, and f2
  blocks then fall back to 2 Mb of physical distance, as in admixtools2.

## PLINK 2: reading VCF

Checked black-box with plink2 `--vcf` and AdmixTools 8 convertf
(tests/golden/formats/generate_multiallelic.sh, generate_haploid.sh).

### Matched, not corrected

- **Haploid calls are homozygous.** `1` (male chrX, for instance) counts
  as two copies of the ALT allele, `0` as none.
- **FILTER is not read.** A record that fails its filters is kept, as
  plink2 keeps it without `--var-filter`; filter the VCF first to drop it.
- **Only biallelic SNPs.** Multiallelic records, indels and symbolic
  alleles are skipped and counted (plink2 `--snps-only just-acgt
  --max-alleles 2`).
