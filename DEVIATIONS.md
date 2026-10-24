# Deviations from Reference Implementations

popmix is checked against reference tools (admixtools2, DATES) on real data.
Sometimes the reference has a bug or an artefact that changes results. popmix
then corrects it by default, but it can still reproduce the reference exactly,
so parity stays provable. This file lists every such correction and sets the
rules for adding one.

## The process

1. **Parity first.** popmix must reproduce the reference exactly before
   anything is corrected. The golden tests run in compat mode at the
   tolerances in CLAUDE.md. A correction is never a way to explain away a
   mismatch that isn't understood: find the cause in the reference source,
   reproduce it, and only then decide whether to correct it.

2. **One flag per deviation.** Each tool has a compat struct with one boolean
   per quirk (`popmix_core::dates::DatesCompat`; admixtools2 will get its own).
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

No deviations. popmix follows admixtools2 2.0.10 defaults: `maxmiss = 0`,
`poly_only = "f2"`, blocks built after filtering, blocks weighted by length
in the jackknife, FST on its own SNP set and blocks, and qpWave/qpAdm fits by
the C++ code path (`cpp = TRUE`, the default).

### Matched, not corrected

- **Ploidy from the first 1,000 SNPs.** A sample is treated as
  pseudo-haploid if it has no heterozygous call among the first 1,000 SNPs
  (`cpp_packedancestrymap_ploidy`). A low-coverage diploid sample with no
  calls there would be misclassified. popmix matches this.
- **qpAdm weight SEs scale by sqrt(blocks - 1).** admixtools2's C++ path
  (`cpp_get_weights_covariance`, the default) and its R path
  (`get_weights_covariance`) scale the leave-one-out weights differently,
  by sqrt(n - 1) against (n - 1)/sqrt(n). popmix matches the default C++ path.
- **FST drops missing blocks per pair.** `fst(f2_dir, pop1, pop2)` reads one
  pair and leaves out that pair's undefined blocks; `fst(prefix)` over many
  populations leaves out a block for every pair if any pair is undefined
  there. popmix's `fst(a, b)` matches the former. The two agree whenever no
  block is undefined, as in every reference config.

### Not yet verified

- **`maxmiss > 0`.** With missing data, admixtools2 sets missing block values
  to 0 before the qpWave/qpAdm jackknife (`est_to_loo_nafix`, with a warning),
  where popmix drops the block. Nothing with `maxmiss > 0` has been checked
  against admixtools2 yet. If the behaviours differ, deciding between them is
  a deviation decision under the process above.
- **qpAdm `popdrop`.** admixtools2 also fits every subset of the sources
  (`popdrop`); popmix does not compute it yet. The reference files record it
  (`qpadm.json`) for when it is added.
- **Constrained qpAdm away from the boundary.** The only constrained reference
  model has weights pinned at (1, 0). popmix solves the same quadratic program
  exactly, but no reference has an interior or multi-source constrained
  solution yet.
