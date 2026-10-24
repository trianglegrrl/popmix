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

### Matched, not corrected

- **Ploidy from the first 1,000 SNPs.** A sample is treated as
  pseudo-haploid if it has no heterozygous call among the first 1,000 SNPs
  (`cpp_packedancestrymap_ploidy`). A low-coverage diploid sample with no
  calls there would be misclassified. popmix matches this.
- **qpAdm weight SEs scale by sqrt(blocks - 1).** admixtools2's C++ path
  (`cpp_get_weights_covariance`, the default) and its R path
  (`get_weights_covariance`) scale the leave-one-out weights differently,
  by sqrt(n - 1) against (n - 1)/sqrt(n). popmix matches the default C++ path.
- **qpGraph admixture weights come from popmix's own optimizer.**
  admixtools2 optimizes them with R's `optim(method = "L-BFGS-B")`
  (finite-difference gradients with step 1e-3, default `factr`) from
  `numstart` random starts. popmix uses a bounded quasi-Newton search with
  central differences, converged more tightly, from its own seeded starts.
  Everything that follows from the weights is reproduced exactly: the f3
  statistics and their covariance, the edge lengths (the same quadratic
  program `qpgraph` passes to quadprog, solved exactly by an active-set
  method), the fitted f3 and the score, all within 3e-12 at the same
  weights. At the optimum the scores agree to 3e-12 too. The weights agree
  to 1e-7 where the score pins them down, and otherwise within admixtools2's
  own spread across its starts. On `dg_6pop` `admix1` the score is so flat
  that double-precision rounding (3e-10 on 258) leaves the weight undetermined
  over about 3e-4: admixtools2's ten starts end between 0.35922 and 0.35948,
  and popmix gives 0.35944. The `low`/`high` columns (the range over the
  starts) depend on the starts and are not compared.

  Over every graph admixtools2 `generate_all_graphs()` makes for 4 and 5
  `dg_6pop` populations with one or two admixtures (2871 graphs,
  `all_graphs.json`), the objective still agrees exactly: fitted at
  admixtools2's best weights, every graph scores admixtools2's score to
  1e-10. The optimized scores do not always agree. Many of these graphs are
  unidentifiable, and their score has several local minima in the weights,
  some a single point on the boundary or a basin narrower than 0.01, so
  the two optimizers, from different random starts, end in different minima:

  | graphs | same score (1e-6) | popmix lower | admixtools2 lower |
  |---|---|---|---|
  | 4 populations, 1 admixture (84) | 59 | 11 | 14 |
  | 5 populations, 1 admixture (1485) | 1310 | 27 | 148 |
  | 4 populations, 2 admixtures (1302) | 825 | 316 | 161 |

  The best graph of each set scores no worse in popmix (2.382203 against
  2.382203; 0.00192 against 0.00406; 0.00095 against 0.00114). The golden
  test `test_qpgraph_optimized_score_of_every_graph` is a strict xfail that
  records this. Matching graph by graph would need admixtools2's starts
  (R's Mersenne Twister after `set.seed`) and R's L-BFGS-B, finite
  differences included.
- **admixtools2 2.0.10 needs igraph 1.x for qpGraph.** It passes vertex
  sequences to `get.edge.ids()`, which igraph 2 rejects. The reference
  generator shims that one call (same vertex ids, as a plain vector) and
  records it; this affects only how the reference is produced.
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
