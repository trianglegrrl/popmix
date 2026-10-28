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
| DATES-0004 | Counts every sample as diploid, so a pseudo-haploid call is two alleles | A pseudo-haploid call is one allele | `compat="dates-v753"` |
| ADMIXTOOLS2-0001 | qpAdm `popdrop` nested path compares siblings, not nested models | A real nested path from the full model | `compat="admixtools2-2.0.10"`; `qpadm --compat admixtools2-2.0.10` |
| ADMIXTOOLS2-0002 | Constrained qpAdm with ≥ 3 sources stops with an error | Returns the full result | nothing to reproduce |
| ADMIXTOOLS2-0003 | `joint_sfs()` drops SNPs whose allele count is not exactly an integer in floating point | Counts them | `compat="admixtools2-2.0.10"` (`joint_sfs`) |
| ADMIXTOOLS2-0004 | `write_dot()` renames nodes on edges (removes "." and "-") but declares them unrenamed | Quoted names, unchanged | `compat="admixtools2-2.0.10"` (`write_dot`); `qpgraph --dot ... --compat admixtools2-2.0.10` |
| ADMIXTOOLS2-0006 | PLINK files with chromosome codes X and Y stop with an error (`auto_only = TRUE`) | Reads them; X and Y are not autosomes | nothing to reproduce |
| ADMIXTOOLS2-0007 | A sample with no heterozygous call among the first 1,000 SNPs is pseudo-haploid, which catches diploid genomes | Pseudo-haploid only with no heterozygous call at any SNP | `pseudohaploid="admixtools2"` (`--pseudohaploid admixtools2`) |
| ADMIXTOOLS2-0009 | qpGraph's L-BFGS-B stops short of a minimum where a drift edge reaches 0 (a kink in the score) | The best start's weights polished by a line search on each weight | `compat="admixtools2-2.0.10"` (`qpgraph`, `find_graphs`); `qpgraph --compat admixtools2-2.0.10` |
| ADMIXTOOLS2-0008 | A statistic with no standard error gets a z anyway: rounding noise for f4(A, A; B, C), Inf from `compare_fits` | Exactly 0 where it is 0 by construction; z and p NaN when the SE is 0 | nothing to reproduce (rounding noise) |
| EIGENSOFT-0001 | smartpca reports the Tracy-Widom density as the p-value beyond its table | The tail | `compat="eigensoft-7.2.1"` |
| PLINK2-0001 | `--make-king-table` labels its two het-hom counts the wrong way round | Counts named as they are | nothing to reproduce (labels only) |
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
  cells (established from its output). Lags beyond n − m also pick up the pairs at the
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
  | aadr_french_q10 (maxdis 0.5) | not run | 130.737 ± 16.669 | unchanged, no lag wraps |

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

- **DATES:** Its FFT correlation divides its output by n for lags 0 … outmax − 1 but
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
  full-data fit moves by 0.008 generations (AADR v66.p1; 0.002 on v62.0),
  and the jackknife, which sums the per-chromosome differences, by 0.16.
- **Tests:** The golden `full_fit`, `jackknife_by_chrom` and
  `jackknife_result` checks pass only with it on (compat). `default_result`
  checks that the full-precision result stays within DATES's reported
  precision wherever no other correction applies and this one moves the
  result less than that precision; on AADR French.DG it moves it more, and
  the check is skipped with the measured shift.

### DATES-0004: every sample counted as diploid

- **DATES:** Source allele frequencies count every sample as diploid: a
  call of 0 or 2 from a pseudo-haploid sample is two copies of one allele.
  Established black-box (benchmarks/cases/dating.py's DATES case, French =
  Sardinian + Russia_Samara_EBA_Yamnaya on AADR v66.p1): popmix counted
  0.046% more pairs than DATES and fit 153.8 generations against 152.9.
  Russia_Samara_EBA_Yamnaya holds 45 pseudo-haploid samples and one diploid
  genome, I6729.DG. Leaving I6729.DG out, or adding a chromosome-23 SNP
  heterozygous in every sample so that popmix takes them all as diploid,
  makes popmix's LD table and fits identical to DATES's.
- **popmix default:** A pseudo-haploid call is one allele, as in every
  other popmix statistic (and admixtools2), so each observed allele weighs
  the same. Counting it twice weighs a pseudo-haploid sample like a diploid
  genome, though it carries half the information. In a population of one
  ploidy the two give the same frequencies.
- **Flag:** `every_sample_diploid`.
- **Affects:** Only source populations that mix pseudo-haploid and diploid
  samples. The frequencies change, so the SNPs that pass the frequency
  filters, the weights and every LD bin.
- **Live numbers** (bench16, chromosomes 1-22, maxdis 0.5, qbin 10; full-data
  fit; jackknife estimate ± SE): DATES = compat 152.925; 153.339 ± 22.306
  (972,605 SNPs); popmix default 153.843; 154.470 ± 22.536 (972,942 SNPs).
- **Tests:**
  - `dating.rs`: `every_sample_diploid_flag_counts_pseudohaploid_calls_as_two_alleles`
    (isolation: the flag alone gives exactly the result with every sample
    marked diploid). `genotype_matrix.rs`:
    `a_pseudohaploid_call_counts_one_allele` (the corrected count).
  - `tests/golden/test_dates_v753.py`, config `yamnaya_chr21_22_q10`
    (chromosomes 21-22 of 74 AADR v66.p1 samples, with I6729.DG among 12
    Yamnaya; `generate_yamnaya_subset.sh`): every compat check, at DATES's
    printed precision.

### Matched, not corrected

- **`mincount` is ignored.** DATES reads `mincount` but ignores it (its output is
  byte-identical with `mincount` 1, 4 or 40): a SNP is used once each source has one called
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
`compat="admixtools2-2.0.10"` in Python (`qpadm`), `--compat
admixtools2-2.0.10` on the CLI (`popmix qpadm`, whose text and JSON show
the popdrop table).

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
- **Impact:** f2 of a population holding such a sample. On AADR v66.p1, with
  NA18606.DG among CHB.DG's samples, f2(CHB.DG, Han.DG) is 0.000420873 by
  default and 0.000419073 with the admixtools2 rule (0.43%, 0.07 SE); with
  HGDP00581.DG in Druze.DG, f2(Druze.DG, Han.DG) is 0.0385412 and 0.0384821
  (0.15%, 0.16 SE). A population of that sample alone would change far more.
  Populations without it change only through the shared SNP set
  (`poly_only`), and in the AADR tests not at all. v66.p1 labels
  HGDP00581.DG Druze (v62.0 left it out, as Ignore_Druze.DG), so the
  dg_11pop admixtools2 goldens hold it; they compare with
  `pseudohaploid="admixtools2"`.
- **Tests:** isolation and oracle on small data
  (`tests/python/test_extract_f2_formats.py`: a diploid sample made
  homozygous over the first 1,000 SNPs is diploid by default, and only
  statistics with its population change), live on AADR
  (`tests/golden/test_ploidy_rule.py`: for NA18606.DG and HGDP00581.DG, only
  that genome is reclassified, only its population's statistics change, and
  the default equals counting every sample diploid).

### ADMIXTOOLS2-0008: a z without a standard error

- **admixtools2:** `f4()` computes f4(A, A; B, C) and f4(A, B; C, C), which
  are 0 by construction, as sums of f2 that cancel only up to rounding: on
  the dg_6pop f2, f4(Mbuti.DG, Han.DG; French.DG, French.DG) is
  -2.9e-20 with SE 1.1e-16 and z -0.0003, so the z is the ratio of two
  rounding errors. `compare_fits()` gives z = Inf and p = 0 when every
  bootstrap replicate has the same score difference (SE 0).
- **popmix:** such an f4 is exactly 0 in every block, so its SE is 0; and
  wherever the SE is 0 (or not finite), z and p are NaN, meaning "no
  information", not "no effect" (z = 0, as popmix's Rust f3 and f4 gave)
  nor "certain" (z = Inf, as its genotype-path statistics gave for a
  constant D). This covers f3 and f4 on f2 blocks, `qpdstat` and `qp3pop`
  on genotypes, `f4ratio` and `compare_fits`.
- **Impact:** only statistics whose SE is 0 or rounding noise: a population
  repeated where the statistic cancels, or a statistic constant across
  blocks. Others are unchanged.
- **Tests:** `tests/python/test_z_without_se.py` (f4 with pop1 = pop2 and
  pop3 = pop4 on the dg_6pop f2, the `qpdstat` CLI, `compare_fits`); Rust
  `f3_f4::tests::z_is_nan_without_an_se`,
  `geno_fstats::tests::z_and_p_are_nan_without_an_se`.

### ADMIXTOOLS2-0009: qpGraph stops short of a zero-bound kink

- **admixtools2:** `qpgraph()` optimizes the admixture weights with R's
  `optim(method = "L-BFGS-B")` and finite-difference gradients. Where a
  fitted drift edge reaches its bound of 0 the score has a kink, and the
  difference quotient straddles it: from every start the search stops
  short. In AADR v66.p1 ancient_5pop admix1 (s -> Sardinian.DG at 0) both
  tools stop at weight 0.4502-0.4503, score 0.677125-0.677129, where a
  fine grid of fixed-weight fits puts the minimum at 0.45087, 0.677059.
- **popmix:** the best start's weights are then polished: each weight in
  turn by a golden-section search over a window of 0.01 around it, kept
  only if the score falls, until a round gains less than 1e-12 relative
  (`qpgraph/polish.rs`). The starts in `opt`, and the weights' low/high
  range over them, stay as L-BFGS-B left them.
- **Impact:** only fits that stalled. Of the one-admixture graphs of the
  find_graphs goldens, 17 of the 4- and 5-leaf graphs score lower than
  admixtools2's and none higher; ancient_5pop admix1 reaches 0.677059.
  Fits L-BFGS-B finishes change only by rounding.
- **Tests:** `qpgraph/polish.rs` (a kinked function from the stall point),
  `tests/golden/test_fstats_admixtools2.py` `qpgraph_polished` (never worse
  than admixtools2's; ancient_5pop admix1 at the minimum) with
  `qpgraph_fit` under compat, `tests/golden/test_find_graphs_admixtools2.py`
  (parity under compat, never worse by default),
  `tests/python/test_cli_consistency.py` (`--compat` on the CLI).

### Matched, not corrected

- **`keepsnps` overrides every other SNP filter.** As its documentation says
  ("Overrides other SNP filtering options"), admixtools2's `extract_f2`
  keeps every listed SNP whatever `maxmiss`, `minmaf`, `maxmaf`,
  `transitions` and `transversions` say, and on chromosome X too despite
  `auto_only` (measured: 9,444 SNPs, chromosome X included, in
  `tests/golden/admixtools2/v2.0.10/snp_filters/keepsnps_maf.json` with
  `minmaf = 0.2, transitions = FALSE`); only polymorphism still applies.
  popmix's `keepsnps` (`--keep-snps`) matches, so a list made for
  admixtools2 gives the same blocks. `excludesnps` (`--exclude-snps`) is
  popmix's own, and applies with `keepsnps` (TASKS.md H1).
- **f2 and f3 with a single pseudo-haploid population.** Its one allele per
  SNP leaves its sampling bias p(1 - p)/(n - 1) unestimable, so f2 with it
  and f3 with it as the target keep its sampling variance: on the dg_6pop
  f2, f3(Chimp.REF; Mbuti.DG, Han.DG) is 0.158 with z 344, while `qp3pop`
  on genotypes gives NaN with n = 0 for such a target. popmix matches both
  values and, on f2 blocks it extracted, warns (`Warning::OneAllele`,
  TASKS.md G7); admixtools2-format f2 directories do not store what tells
  (`tests/python/test_one_allele_warning.py`).
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
  | 4 populations, 2 admixtures (1302) | 1222 | 41 | 39 |

  On AADR v66.p1, the 80 two-admixture graphs that differ are ones whose
  score is nearly flat along a ridge of local minima; the score's own
  rounding (1e-12 to 1e-11 apart from admixtools2's, from inverting the f3
  covariance and the quadratic program) changes a finite-difference
  gradient at its eighth digit, and a path that starts identical to admixtools2's ends in a
  neighbouring minimum. The best graph is the same in every set (2.116214,
  0.00397629, 0.00113061). Before popmix used R's starts and L-BFGS-B, its
  own quasi-Newton search from its own starts agreed on 59, 1310 and 825 of
  these graphs (on AADR v62.0).
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
- **Projected samples are shrunk toward the origin.** A PCA sample's
  coordinates include its own noise, which helped fit the axes; a projected
  sample's do not. Both smartpca's default projection and popmix's place
  projected samples closer to the origin than PCA samples of their own
  population. With six 1000 Genomes populations in AADR v66.p1 (CEU, TSI,
  FIN, YRI, CHB, GIH; about 290 PCA samples), a projected half of each
  reaches 98% of its PCA half's spread on PC1, 95% on PC2, and 74%, 55%,
  62% and 28% on PC3-PC6 (`tests/golden/test_pca_holdout_aadr.py`); on the
  small test data, 82% on PC1. popmix keeps the projection as it is, so its
  coordinates stay smartpca's, and adds `holdout` (`--holdout`): every
  second sample of each PCA population is projected too, projected samples
  are compared with reference samples placed the same way, and each
  component's shrinkage is measured from the halves, with a warning below
  0.9. smartpca's `shrinkmode` is not documented, so popmix does not
  reproduce it.
  A model-based correction was tried and removed before release: dividing
  each significant component by ρ = 1 − Σᵢ nᵢuᵢ²/λ, nᵢ a PCA sample's
  noise estimated as its diagonal less its share of the Tracy-Widom-
  significant components. It came within 5% on PC1 of the small data, but on
  the AADR panel it estimated 0.997 for PC5 and PC6 against a measured 0.62
  and 0.28, and undercorrected PC1 and PC3. Every significant component
  takes some of each sample's own noise with it, so the noise estimate came
  out too small. Estimating the noise from population mates instead gave
  0.79 against 0.82 on the small data. Measuring the shrinkage needs no
  model.

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
  data and each chromosome left out). On AADR v62.0 their decays agreed
  within 7.7e-6 and the jackknife within 3.4e-6; on v66.p1 the
  chromosome-left-out decays agree within 6.4e-6, but the all-data decay
  differs by 1.3e-5 (6.512929 against nls's 6.512847, the squared errors
  equal to 6 digits) and the jackknife by 6.4e-5, beyond the golden
  test's 1e-5.
- **expfit.sh's jackknife assumes chromosomes 1-22.** wtjack.pl matches
  per-chromosome fits to chromosomes 1 to 22 in order; with a subset of
  chromosomes its weights and estimates are misaligned. popmix weights each
  chromosome present by its SNPs, as wtjack.pl intends, and the reference
  is generated the same way.


## ADMIXTOOLS 7.0.2 (qpAdm)

popmix's qpAdm follows admixtools2 (`qpadm()`), which is what its goldens
check. The benchmark case `classic_qpadm` (benchmarks/cases/classic.py) runs
classic qpAdm (version 1520, `allsnps: NO`, `inbreed: NO`) on the same
genotypes, with `--pseudohaploid diploid`, and the two differ beyond
qpAdm's printed precision. This is a difference of method, not a bug in
either, and popmix does not reproduce it.

### Differences, not reproduced

- **Which SNPs are used.** Established black-box on bench16 (Spain_EN =
  Turkey_N + Serbia_IronGates_Mesolithic, 7 right): classic qpAdm keeps a
  SNP when every basis statistic f4(L0, Li; R0, Rj) (L0 the target, R0 the
  first right population) can be evaluated, where a statistic with one
  factor exactly 0 (pL0 = pLi, or pR0 = pRj) counts as 0 even when the
  other factor's population has no data there. Monomorphic SNPs are kept.
  admixtools2 and popmix keep the SNPs every population has data at, less
  the monomorphic ones. The rule gives classic qpAdm's "numsnps used"
  exactly: 803,891 on the model (popmix 675,980; 725,784 with monomorphic
  SNPs), 1,074,466 with Israel_Natufian (the population with least
  coverage) left out of right, and 12,585 on chromosome 22 alone; and of
  320 SNPs of chromosome 22 that only Israel_Natufian lacks, added one at a
  time to the complete ones, it predicts which 52 qpAdm counts. Such a SNP
  adds nothing to an f4 sum but counts in its block, so it shrinks that
  block's f4 means. The rule also selects SNPs by their genotypes (left
  populations fixed for one allele), and it is not expressible with f2
  statistics, which need every population's frequency.
- **Standard errors and chi-squared.** With the SNP sets made the same (classic
  qpAdm given `badsnpname` for every SNP some population lacks: 725,784
  SNPs; popmix with `poly_only=False`), the weights agree to the printed
  precision (0.814 / 0.186; jackknife mean 0.813668, popmix 0.813619), but
  the standard error is 0.021 (error variance 434e-6) against popmix's and
  admixtools2's 0.0218, and chi-squared 12.582 against 12.631 (12.678 with
  `oracle: NO`, `diagplus: 0`). qpAdm's own jackknife and covariance
  estimate (oracle shrinkage) differ from admixtools2's.

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
