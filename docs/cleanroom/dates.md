# Clean-room specification: DATES (`crates/popmix-core/src/dates/`)

This specification describes, in mathematical and behavioural terms, what
popmix's DATES module computes: SNP weights, the weighted LD decay table in
both of DATES's modes (every pair, and the FFT mesh), the exponential fit,
the leave-one-chromosome-out weighted jackknife, and the four DATES v753
behaviours that popmix corrects by default (DATES-0001 to DATES-0004).

It was written from Narasimhan et al. 2019 (Science 365, eaat7487,
supplementary methods of DATES), Chintalapati, Patterson & Moorjani 2022
(eLife 11:e77625), the DATES README and parameter-file documentation,
Busing, Meijer & van der Leeden 1999 (Statistics and Computing 9:3-8, the
weighted delete-m jackknife), `DEVIATIONS.md` (DATES-0001 to DATES-0004),
and DATES v753 run as a black box (its outputs are the golden references in
`tests/golden/dates/v753/`). DATES's source carries no licence; the module
is written from this document and the tests alone, by someone who has not
read DATES's source or the earlier module. It contains no code.

Two helpers are already clean-room and stay as they are:
`mesh_cells` and `MeshLagSums` (with `new`, `merge`, `to_bins`) in
`ld_decay.rs`, specified in `docs/cleanroom/dates-mesh.md`. The
implementer may keep that code verbatim (it is the only part of the old
module to keep) or rewrite it from that document; either way its
behaviour is as described there and summarised in section 4.

Record the provenance (this document and the papers above) in each new
file's module doc comment (CLAUDE.md, "Provenance").

## 1. Overview

DATES dates an admixture event from the decay, with genetic distance, of
the ancestry covariance in an admixed (test) population. Inputs: a genotype
matrix, the name of the test population, and two source populations.

1. **Weights.** For each SNP, w = p1 − p2, the difference of the source
   allele frequencies; SNPs that are not informative enough are dropped
   (section 3).
2. **Per-individual residuals.** For each test individual, regress its
   genotypes on the source frequencies to get its admixture proportion α,
   and take the weighted residual y = (x − α p1 − (1 − α) p2) · w at each
   SNP (section 4.2).
3. **LD table.** For SNP pairs on the same chromosome, accumulate y1·y2 and
   friends into distance bins, summed over test individuals; report a
   correlation per bin (section 4).
4. **Fit.** Fit A·r^k + K to the correlations from `lovalfit` on; the
   number of generations is −ln(r)/binsize (section 5).
5. **Jackknife.** Repeat the fit with each chromosome left out and combine
   with the weighted block jackknife, each chromosome weighted by its SNP
   count (section 6).

Units: the public options are in centimorgans (cM); internally binsize and
maxdis are in Morgans, binsize = `bin_width_cm` / 100 and maxdis =
`max_distance_cm` / 100, each computed by that division in f64. SNP
genetic positions are `Snp::position_morgans` (Morgans).

## 2. Public API to preserve exactly

`mod.rs` declares the private submodules `compat`, `dating`,
`exponential_fit`, `ld_decay`, `weights` (names may change, the re-exports
may not) and re-exports:

```
pub use compat::DatesCompat;
pub use dating::{dates, DatesOptions, DatesResult, JackknifeBlock};
pub use exponential_fit::{fit_exponential, weighted_jackknife, DecayFit, FitOptions};
pub use ld_decay::{compute_ld_decay, compute_ld_decay_by_chromosome, ChromosomeLd, LdBin, LdDecay, LdOptions};
pub use weights::{compute_weights, WeightResult};
```

`crates/popmix-core/src/lib.rs` re-exports at the crate root:
`compute_ld_decay, compute_weights, dates, fit_exponential,
weighted_jackknife, DatesOptions, DatesResult, DecayFit, LdBin, LdOptions,
WeightResult` (it stays unchanged).

(The listings in this section are signatures, the interface the callers
compile against, not implementation.)

### Types

All structs below derive `Debug, Clone`; all fields are `pub`.

- `DatesCompat` — derives `Debug, Clone, Copy, Default, PartialEq, Eq`.
  Four `bool` fields, in this order: `circular_fft_lags` (DATES-0001),
  `unscaled_last_fft_lag` (DATES-0002), `fit_printed_values` (DATES-0003),
  `every_sample_diploid` (DATES-0004). Associated constants
  `CORRECTED` (all false; equal to `Default`) and `V753` (all true).
- `LdOptions { bin_width_cm: f64, max_distance_cm: f64,
  min_pairs_per_bin: usize, qbin: usize, compat: DatesCompat }`.
  `Default`: 0.1, 100.0, 0, 0, `DatesCompat::CORRECTED`.
  Method `pub fn validate(&self) -> Result<()>` (section 4.1).
- `LdBin { distance_cm: f64, weighted_ld: f64, z_score: f64,
  correlation: f64, regression_coef: f64, n_pairs: usize, std_err: f64 }`.
  Built with struct literals by the golden test, so no other fields.
- `LdDecay { bins: Vec<LdBin>, chromosomes: Vec<ChromosomeLd> }`.
- `ChromosomeLd { chromosome: u8, n_snps: usize, leave_one_out: Vec<LdBin> }`.
- `FitOptions { min_distance_cm: f64, max_distance_cm: f64,
  use_affine: bool, bin_size_cm: f64 }`. `Default`: 0.45, 100.0, true, 0.1.
- `DecayFit { amplitude: f64, decay_rate: f64, intercept: f64,
  r_squared: f64, decay_rate_se: f64, rss: f64 }`.
- `WeightResult { weights: Vec<f64>, n_valid: usize, valid_indices: Vec<usize> }`.
- `DatesOptions { ld_options: LdOptions, fit_options: FitOptions,
  generation_time: f64, min_source_count: u32, run_jackknife: bool,
  compat: DatesCompat }`. `Default`: `LdOptions::default()`,
  `FitOptions::default()`, 29.0, 1, true, `DatesCompat::CORRECTED`.
- `DatesResult { test_pop: PopulationId, source1: PopulationId,
  source2: PopulationId, generations: f64, generations_se: f64,
  years_bp: f64, years_bp_se: f64, ld_bins: Vec<LdBin>, fit: DecayFit,
  n_snps: usize, weights: WeightResult, n_chromosomes: usize,
  jackknife_by_chromosome: Vec<JackknifeBlock> }`.
- `JackknifeBlock { chromosome: u8, weight: f64, generations: f64 }`.

### Functions

```
pub fn compute_weights(af: &AlleleFrequencies, source1: &str, source2: &str,
                       min_count: Option<u32>) -> Result<WeightResult>
pub fn compute_ld_decay(geno: &GenotypeMatrix, af: &AlleleFrequencies, test_pop: &str,
                        weights: &[f64], source1_freqs: &[f64], source2_freqs: &[f64],
                        options: Option<LdOptions>) -> Result<Vec<LdBin>>
pub fn compute_ld_decay_by_chromosome(/* same arguments */) -> Result<LdDecay>
pub fn fit_exponential(ld_bins: &[LdBin], options: Option<FitOptions>) -> Result<DecayFit>
pub fn weighted_jackknife(full_estimate: f64, estimates: &[(f64, f64)]) -> Result<(f64, f64)>
pub fn dates(geno: &GenotypeMatrix, af: Option<&AlleleFrequencies>, test_pop: &str,
             source1: &str, source2: &str, options: Option<DatesOptions>) -> Result<DatesResult>
```

`Result` is `crate::Result`, errors `crate::error::PopmixError`.
`compute_ld_decay` returns `compute_ld_decay_by_chromosome(..)?.bins`. In
the LD functions `af` is accepted and unused (kept for API compatibility).

### Callers

- `crates/popmix-python/src/lib.rs` (`dates` pyfunction, `PyLdBin`,
  `PyDecayFit`, `PyDatesResult`): builds `LdOptions::default()` and sets
  `qbin`, `bin_width_cm`, `max_distance_cm`; builds
  `dates::FitOptions::default()` and sets `bin_size_cm` (= the bin width),
  `min_distance_cm` (lovalfit), `max_distance_cm` (= maxdis in cM); builds
  `DatesOptions { ld_options, fit_options, generation_time, run_jackknife,
  compat, ..Default::default() }`; maps `compat=None` to
  `DatesCompat::CORRECTED` and `"dates-v753"` to `DatesCompat::V753`
  (anything else is a `ValueError` naming "dates-v753"); calls
  `popmix_core::dates(&geno, None, ...)` and reads every field of
  `DatesResult` (`test_pop`, `source1`, `source2` via `to_string()`),
  every field of `LdBin` and `DecayFit`, and `chromosome`, `weight`,
  `generations` of each `JackknifeBlock`. Errors map to Python exceptions
  by variant (`Validation`, `PopulationNotFound`, ... → `ValueError`).
- `crates/popmix-python/src/rolloff.rs`: calls
  `popmix_core::dates::weighted_jackknife(full, &blocks)` for ROLLOFF's
  jackknife (`tests/python/test_rolloff.py` depends on it).
- `python/popmix/__init__.py` re-exports `dates`, `LdBin`, `DecayFit`,
  `DatesResult`; `python/popmix/cli/__init__.py` (`popmix dates`) passes
  `bin_width_cm`, `max_distance_cm`, `min_distance_cm`, `qbin`, `compat`
  and reads `fit.amplitude`, `fit.decay_rate`, `fit.r_squared`,
  `generations(_se)`, `years_bp(_se)`, `n_snps`, `ld_bins`.
- Crate helpers outside the module that the module uses and that stay:
  `AlleleFrequencies::from_genotypes_counting(geno, Some(&pops), all_diploid)`,
  `AlleleFrequencies::pop_index`, `.frequencies[[snp, pop]]`,
  `.counts[[snp, pop]]` (observed allele count, `u32`), `.n_snps()`;
  `GenotypeMatrix::population_indices`, `snps()`, `snp_raw(snp)`,
  `n_snps()`; `crate::options::check_bins`, `check_positive`;
  `crate::threading::ordered_sum` (deterministic ordered reduction, see
  section 4.6); the `rustfft` dependency.

## 3. Weights (`compute_weights`)

Look up both sources in `af`; a missing one is
`PopulationNotFound(name)`. `min_count` defaults to 4 when `None`
(`dates` always passes `Some(min_source_count)`, default 1).

For every SNP i of `af` (in order), with c1, c2 the observed allele counts
and p1, p2 the frequencies of source1 and source2, the SNP is **valid**
when all of:

- c1 ≥ min_count and c2 ≥ min_count;
- p1 and p2 are finite;
- |p1 − p2| ≥ 0.001;
- with m = (p1 + p2) · 0.5, m · (1 − m) ≥ 0.001.

`weights[i]` = p1 − p2 for valid SNPs and NaN otherwise (length = number
of SNPs); `valid_indices` lists the valid SNPs in increasing order;
`n_valid` is its length. Swapping the sources negates every weight.

The `mincount` par option of DATES is read but ignored by DATES v753 (a
SNP is used once each source has one called sample); popmix matches this
with `min_source_count = 1` ("Matched, not corrected" in DEVIATIONS.md).

## 4. LD decay (`compute_ld_decay_by_chromosome`)

### 4.1 Validation and sizes

1. `options` defaults to `LdOptions::default()`; call `validate()`:
   `check_bins("bin_width_cm", bin_width_cm, "max_distance_cm",
   max_distance_cm)` (both finite and > 0, at most 1e6 bins; messages name
   the option), then if `qbin · max_distance_cm / bin_width_cm` (f64,
   left to right) exceeds 1e7, a `Validation` error whose message starts
   "qbin {qbin} makes {cells:.0} mesh cells; at most 10000000 are allowed".
   (`tests/python/test_bad_options.py` matches "bin_width_cm" and "qbin".)
2. The test population's samples are `geno.population_indices(test_pop)`,
   in sample order; none is `PopulationNotFound(test_pop)`.
3. Number of reported bins: n_bins = round(maxdis / binsize) (round half
   away from zero). Bin k (k = 0 … n_bins − 1) holds pairs at distance in
   [k · binsize, (k + 1) · binsize). Pairs that would fall at bin index
   n_bins or beyond are never reported anywhere.

### 4.2 SNPs used and residuals

Only chromosomes 1 … 22 (`Chromosome::Autosome(n)` with 1 ≤ n ≤ 22) take
part; SNPs on any other chromosome are skipped (but still count in the
mesh axis, 4.4).

- **Chromosomes present**: those of 1 … 22 having at least one SNP in the
  input, in increasing order. For each, its **SNP count** is the number of
  its SNPs in the input before any filtering.
- **SNPs used**: SNP i on a chromosome 1 … 22 is used when `weights[i]` is
  finite with |weights[i]| ≥ 0.001 and, with m = (source1_freqs[i] +
  source2_freqs[i]) · 0.5, m(1 − m) ≥ 0.001. (The same filter as section 3,
  applied again because the function is public and takes arbitrary
  weights.)
- If any used SNP has a non-finite `position_morgans`, fail with
  `Validation("SNP {id} has no genetic position ({pos}); DATES needs one for
  every SNP")` for the first such SNP (in chromosome order, then input
  order). The message must contain the SNP id.
- Within each chromosome, the used SNPs are ordered by increasing
  `position_morgans`; ties keep input order (a stable sort).

For each test individual j (raw genotype g ∈ {0, 1, 2} = count of the
counted allele; a raw value > 2 is missing):

- x = g / 2 for each used SNP where j is called (whatever the sample's
  ploidy: a pseudo-haploid call 0 or 2 gives 0 or 1).
- α = Σ (x − p2) · w / Σ w², both sums over the SNPs where j is called (all
  chromosomes, in chromosome order then position order), with w the weight
  and p2 the source2 frequency. If Σ w² ≤ 1e-15, α = 0.5.
- At each called used SNP, the residual value is
  y = (x − (α · p1 + (1 − α) · p2)) · w.

Missing calls contribute nothing (no value, no pairs).

### 4.3 Exact mode (`qbin = 0`, the DATES default)

For each individual and chromosome, visit every pair of its valued SNPs
(i before j in the ordering of 4.2) with distance
dis = pos_j − pos_i < maxdis. Since positions are sorted, the scan from i
stops at the first j with dis ≥ maxdis. The pair's bin is
trunc(dis / binsize). Into that bin add one pair with **y1 = the earlier
SNP's value (i), y2 = the later SNP's value (j)**:

  s0 += 1, s1 += y1, s2 += y2, s11 += y1², s12 += y1·y2, s22 += y2².

Pairs are added in the order i ascending, then j ascending.

### 4.4 Mesh mode (`qbin > 0`)

The mesh trades exactness for speed: values are placed on cells of width
c = binsize / qbin, and per-lag sums over cell pairs are computed by FFT.

1. **Cells.** `mesh_cells(geno.snps(), binsize, qbin)` over **all** SNPs of
   the genotype matrix in input order (not just used ones; see
   `dates-mesh.md`: the axis advances by the position difference within a
   chromosome and by exactly 5 Morgans at each chromosome change).
2. **Maximum lag.** diffmax = round(qbin · maxdis / binsize) (f64, left to
   right; round half away from zero). Lag sums cover lags 0 … diffmax.
3. **Per individual and chromosome.** Take the individual's valued SNPs
   on that chromosome. If there are none, all lag sums are zero. Otherwise
   let lo and hi be the smallest and largest cell among them and
   m = hi − lo + 1 the window; if m ≤ 1 all lag sums are zero. Otherwise
   build three arrays of length m indexed by cell − lo: the number of
   values, the sum of y, and the sum of y² in each cell (values added in
   the order of 4.2).
4. **FFT length.** n = 2^(⌈log₂ m⌉ + 1), with log₂ m computed in f64 and
   the result rounded to an integer (so m = 50 gives n = 128, m = 64 gives
   n = 128). n ≥ 2m always. The **same n is used whatever the compat
   flags**, so that lags unaffected by a flag come out bit-identical with
   and without it (the golden `default_ld_table` check and
   `tests/python/test_dates.py` compare them with exact equality).
5. **Circular cross-correlations.** For arrays a, b zero-padded to n, the
   circular cross-correlation is C_ab[d] = Σ_i a[i] · b[(i + d) mod n],
   d = 0 … n − 1, computed as inverse-FFT(conj(FFT a) · FFT b), real part,
   times 1/n (unnormalised forward and inverse transforms, as `rustfft`
   gives). With N, Y, Q the count, sum and sum-of-squares arrays, compute
   C_NN, C_NY, C_YY and C_NQ.
6. **Lag sums** (the six of `dates-mesh.md`), for lag d:

   | sum | value at lag d |
   |---|---|
   | pairs | C_NN[d] |
   | later | C_NY[d] |
   | earlier | C_NY[(n − d) mod n] |
   | cross | C_YY[d] |
   | later_sq | C_NQ[d] |
   | earlier_sq | C_NQ[(n − d) mod n] |

   Which lags are filled:
   - **Corrected (default):** lags 0 … min(diffmax, m − 1). Every real
     pair has lag at most m − 1, and since n ≥ 2m no such lag wraps. All
     other lags are zero.
   - **`circular_fft_lags` (DATES-0001):** lags 0 … min(diffmax, n − 1),
     with the circular values above. Lags d > n − m then also include
     the pairs at negative lag n − d (short-range pairs counted again at
     long distance). Lags ≥ n are zero.
   - **`unscaled_last_fft_lag` (DATES-0002):** at the single lag
     d* = min(n − 1, diffmax), and only if that lag is among the lags
     filled, **pairs** and **cross** (only these two) are multiplied by n,
     i.e. left without the 1/n scaling. The four mixed sums (later,
     earlier, later_sq, earlier_sq) are unaffected. Without
     `circular_fft_lags` the filled lags reach d* only when diffmax ≤ m − 1,
     and lag diffmax falls in bin n_bins or near it, normally unreported;
     so on its own this flag is normally invisible.
7. **Accumulate over individuals** at the level of lag sums: per chromosome,
   the individuals' lag sums are added element-wise in sample order
   (section 4.6). Only then convert each chromosome's summed lag sums to
   bins with `to_bins(binsize, qbin, n_bins_or_more)` (`dates-mesh.md`):
   lags d ≥ 1 with pairs ≥ 0.5, bin trunc(d · (binsize / qbin) / binsize);
   pairs→s0, later→s1, earlier→s2, cross→s12, later_sq→s11,
   earlier_sq→s22. So in mesh mode **y1 is the later SNP of a pair** (the
   opposite of exact mode). This matters for s11 and s22 and so for the
   regression column; both are as DATES reports them in each mode.

Lag 0 is never used: pairs within one cell (less than binsize / qbin
apart) are not counted, matching DATES ("Matched, not corrected").

### 4.5 From sums to the table

Per chromosome present, the result of 4.3 or 4.4 is an array of bin sums
(s0, s1, s2, s11, s12, s22). The **total** is the element-wise sum over
chromosomes in increasing chromosome order, starting from zero. The
**leave-one-out** sums for chromosome c are total − sums(c), element-wise
(a subtraction, not a re-summation).

A bin's sums become an `LdBin` (bin index k):

- `distance_cm` = (k + 1) · `bin_width_cm` (f64; the bin's upper edge, as
  DATES prints it).
- If s0 < 0.5: `weighted_ld` = `correlation` = `z_score` = 0. Otherwise,
  with v11 = s11/s0, v12 = s12/s0, v22 = s22/s0:
  `weighted_ld` = v12; `correlation` = v12 / sqrt(v11 · v22 + 1e-20);
  `z_score` = sqrt(s0) · correlation. (No mean is subtracted.)
- `regression_coef` = s12 / (s11 + 1e-20), whatever s0.
- `n_pairs` = s0 rounded to the nearest integer (half away from zero).
- `std_err` = 1 / sqrt(s0 − 3) if s0 > 3, else NaN.

Reported bins: indices 0 … n_bins − 1, in order, dropping any bin with
s0 < `min_pairs_per_bin` (as f64). With the default 0 every bin is
reported, including empty ones as zeros, as DATES does.

Output: `LdDecay.bins` = the reported total; `LdDecay.chromosomes` = one
`ChromosomeLd` per chromosome present, in increasing order, with its SNP
count (4.2) and its reported leave-one-out bins. With one chromosome, its
leave-one-out bins all have `n_pairs` 0.

Relation to DATES's `.out` file (`tests/golden/dates/v753/*/dates.out`):
columns are distance (cM, `%9.3f`), `weighted_ld` (`%15.9f`),
`regression_coef` (`%12.6f`), `correlation` (`%12.6f`), `n_pairs`; the
`.out` of DATES has no z-score or std_err column.

### 4.6 Summation order and threads

Results must not depend on the thread count, to the last bit
(`bins_do_not_depend_on_the_thread_count`). Each test individual's
contribution (per-chromosome bin sums in exact mode, per-chromosome lag sums
in mesh mode) is computed independently (in parallel is fine), and the
contributions are added into an accumulator that starts at zero, one
individual at a time, in sample order. `crate::threading::ordered_sum`
does exactly this (it evaluates batches of `current_num_threads()`
individuals in parallel and adds them in order). Within one individual the
order is that of 4.3 / 4.4.

## 5. Exponential fit (`fit_exponential`)

DATES's `dates_expfit` (one exponential, affine mode, as run by DATES with
`afffit: YES`) fits the correlation column. popmix minimises the same
objective to full precision.

**Bins fitted.** Options default to `FitOptions::default()`. Take
`ld_bins` in the order given, skip leading bins while
`distance_cm < min_distance_cm`, then take bins while
`distance_cm ≤ max_distance_cm`; the fitted values y_0 … y_{M−1} are their
`correlation`s. The exponent is the position in this list, not the
distance (bins dropped by `min_pairs_per_bin` are simply absent).

- M < 4: `InsufficientData { context: "exponential fit (bins in fitting
  range)", required: 4, actual: M }`.
- Any non-finite y: `InvalidValue { context: "exponential fit", details:
  "non-finite correlation in fitted bin {index}" }`.

**Model.** y_k ≈ A · r^(k+1) + K, k = 0 … M − 1, with 0 < r < 1 (r is the
per-bin decay; the first fitted bin has x_0 = r). For a fixed r, with
x_k = r^(k+1) computed by repeated multiplication (x_0 = r, x_k =
x_{k−1} · r), A and K are the ordinary least-squares solution:

- affine (`use_affine`): with S_x = Σx, S_y = Σy, S_xx = Σx², S_xy = Σxy
  (sums in index order) and det = M·S_xx − S_x², A = (M·S_xy − S_x·S_y)/det
  and K = (S_xx·S_y − S_x·S_xy)/det; if |det| < 1e-300 the objective at r
  is +∞;
- not affine: K = 0, A = S_xy/S_xx (A = 0 if S_xx ≤ 1e-300).

The objective is the mean squared error MSE(r) = Σ (y_k − A x_k − K)² / M.

**Search.** Minimise MSE over r ∈ [e^(−20), e^(−1e-7)] (decay per bin
u = −ln r from 1e-7 to 20) and return the global minimum to full
precision. The current method, which the tests are known to pass with:
evaluate MSE on 4000 values of r, r = exp(−u) with ln u equally spaced
from ln 1e-7 to ln 20 (so r decreases along the grid); every grid point
whose MSE is not greater than either neighbour's (the two ends compared
with +∞ outside) is refined by golden-section search on the bracket
between its two grid neighbours (the point itself at an end), stopping
when the bracket is no wider than 1e-15 · |upper end| or after 200
iterations, and returning the better of the two interior points; the
refined r with the smallest MSE wins (the earliest on the grid in a tie).
If no finite r results: `InvalidValue { context: "exponential fit",
details: "no minimum of the fit error found" }`. Any other optimiser is
acceptable if it finds the same global minimum to about 1e-9 relative in
the decay rate (see the tolerances in section 8).

**Result.** With r̂ the minimiser and A, K at r̂, step = `bin_size_cm` /
100:

- `decay_rate` = −ln(r̂) / step (generations);
- `amplitude` = A · r̂ (the fitted excess at the first fitted bin);
- `intercept` = K (0 if not affine);
- `rss` = MSE(r̂) (the mean, despite the name);
- `r_squared` = 1 − MSE · M / SS_tot with SS_tot = Σ (y − ȳ)², ȳ = Σy/M;
  0 if SS_tot is not > 0;
- `decay_rate_se` = NaN (filled in by `dates` when the jackknife runs).

## 6. Weighted jackknife (`weighted_jackknife`)

The weighted delete-one-block jackknife of Busing et al. 1999, as DATES
uses it over chromosomes (DATES's `dowtjack` program, whose inputs and
outputs are the `.jin` and `.jout` files).

Input: θ̂ = `full_estimate` and, per block, (w_k, θ_k) = (weight, estimate
with block k left out). Blocks with w_k < 1e-6 are ignored. With g blocks
left: g ≤ 1 is `InsufficientData { context: "weighted jackknife (blocks
with positive weight)", required: 2, actual: g }`. Otherwise, with
N = Σ w_k:

- estimate θ_J = Σ_k (θ̂ − θ_k) + Σ_k w_k θ_k / N;
- h_k = N / w_k; pseudo-value τ_k = h_k θ̂ − (h_k − 1) θ_k;
- variance = (1/g) Σ_k (τ_k − θ_J)² / (h_k − 1);
- returns (θ_J, sqrt(variance)) if variance > 0, else (θ_J, 0).

With equal weights this is the ordinary delete-one jackknife:
θ_J = g θ̂ − (g − 1) mean(θ_k), variance = (g − 1)/g · Σ (θ_k − mean)².
Ignored zero-weight blocks change nothing, bit for bit.

## 7. The analysis (`dates`)

1. Options default to `DatesOptions::default()`.
   `check_positive("generation_time", generation_time)` first (finite and
   > 0; the message names `generation_time`).
2. Allele frequencies: `af` if given, used as is. Otherwise
   `AlleleFrequencies::from_genotypes_counting(geno, Some(&[test_pop,
   source1, source2]), compat.every_sample_diploid)`.
3. `compute_weights(af, source1, source2, Some(min_source_count))`. If
   `n_valid` is 0: `InsufficientData { context: "No SNPs with valid
   weights", required: 100, actual: 0 }`.
4. Source frequency vectors (one value per SNP of `geno`) from
   `af.frequencies` at the two sources' indices.
5. `compute_ld_decay_by_chromosome` with `opts.ld_options` but its `compat`
   replaced by `opts.compat` (the analysis-wide flags override
   `ld_options.compat`). If the reported total has no bins:
   `InsufficientData { context: "No LD bins computed", required: 10,
   actual: 0 }`.
6. Full-data fit of the total bins (section 5 with `opts.fit_options`).
   With `fit_printed_values` (DATES-0003), the correlations are first
   replaced by their values rounded to 6 decimals.
7. If `run_jackknife`:
   - for every chromosome present (in parallel is fine; keep increasing
     chromosome order), fit its leave-one-out bins the same way (rounded
     to 6 decimals under DATES-0003). Any failing fit fails the whole call.
     `JackknifeBlock { chromosome, weight: SNP count as f64, generations:
     that fit's decay_rate }` (unrounded).
   - `weighted_jackknife(θ̂, [(weight, generations) per block])` with
     θ̂ = the full fit's `decay_rate`; under DATES-0003 θ̂ and each block's
     generations are first rounded to 3 decimals.
   - `generations`, `generations_se` = the jackknife estimate and SE;
     `fit.decay_rate_se` = the SE; `fit.decay_rate` stays the unrounded
     full-data fit; `jackknife_by_chromosome` = the blocks.
   Otherwise `generations` = the full fit's `decay_rate`,
   `generations_se` = NaN, and no blocks.
8. `years_bp` = generations · generation_time; `years_bp_se` =
   generations_se · generation_time if generations_se is finite, else NaN.
9. `n_snps` = `weights.n_valid` (all chromosomes); `n_chromosomes` = the
   number of distinct autosomal chromosomes (`Snp::is_autosomal`) among
   all SNPs of `geno`; `ld_bins` = the total bins; `weights` = the
   `WeightResult`; names as `PopulationId::new(..)`.

Rounding "to d decimals" means the correctly rounded decimal of the binary
value, as C's `printf("%.{d}f")` gives on glibc and as Rust's `{:.d}`
formatting gives (format, then parse back).

Errors propagate in that order: a missing source is reported by step 3, a
missing test population by step 5.

## 8. Compat flags (DATES-0001 … DATES-0004)

Each flag reproduces one DATES v753 behaviour that popmix corrects by
default; `DatesCompat::V753` turns them all on and must then reproduce
DATES's output (DEVIATIONS.md has the evidence and live numbers).

| Flag | DATES v753 | popmix default | What it changes |
|---|---|---|---|
| `circular_fft_lags` (0001) | Mesh lags past n − m pick up wrapped short-range pairs (circular FFT) | Lags past the window are zero | Mesh only; a tail of the table from the first wrapped lag to maxdis |
| `unscaled_last_fft_lag` (0002) | Pair and cross sums at lag min(n − 1, diffmax) left n times too large | Every lag scaled | Mesh only; one bin, visible only together with 0001 |
| `fit_printed_values` (0003) | Fits correlations as printed (`%12.6f`); jackknife on estimates as printed (`%9.3f`) | Full precision | The estimates, by rounding noise |
| `every_sample_diploid` (0004) | Every sample counted diploid in source frequencies: a pseudo-haploid call is two alleles | A pseudo-haploid call is one allele | Only sources mixing ploidies: their frequencies, so the SNPs used, weights and every bin |

0001 and 0002 apply wherever `LdOptions.compat` is read (section 4.4);
0003 and 0004 only in `dates` (0004 only when `dates` computes the
frequencies itself).

## 9. Tests that exercise the module

Tolerances fixed by DATES's printed precision (half a unit in the last
digit, plus a little): covariance 5e-10, regression and correlation 5e-7,
generations 5e-4; for fits add DATES's optimiser tolerance 1e-6 in r,
i.e. 1e-6 / binsize generations (DATES stops early, popmix goes to the
true minimum). Pair counts must match exactly.

### Rust unit tests in the module

These live in the old files' `#[cfg(test)]` modules (test code, which may
be handed to the implementer as is; otherwise rebuild them from this).
Fixtures use `Sample::new(id, PopulationId::new(pop))`, `Snp::new(id,
Chromosome::Autosome(c), bp, morgans, 'A', 'G')`,
`GenotypeMatrix::new(data (snps × samples, u8), samples, snps)` and
`AlleleFrequencies::from_genotypes`.

`weights.rs` (3 samples each in pop1 and pop2, 10 SNPs on chr 1;
pop1 genotypes 2, 2, 1, pop2 0, 1, 0 at every SNP):
- `test_weight_computation`: `min_count` 2 → all 10 valid, weights ≈ 4/6.
- `test_weight_signs`: swapping the sources negates the weights (1e-10).
- `test_missing_population`: unknown source → error.

`ld_decay.rs` (fixture: 10 "test", 5 "source1", 5 "source2" samples; 50
SNPs on chr 1, bp i·10000, position i·spacing Morgans; genotype
(snp + sample) mod 3; weights p1 − p2 from the frequencies):
- `bins_do_not_depend_on_the_thread_count`: 64 test + 8 + 8 source
  samples, 200 SNPs (chr 1 + i/100, position (i mod 100)·0.0031),
  xorshift genotypes; bin width 0.05 cM, max 2 cM, `min_pairs_per_bin` 1;
  for qbin 0 and 10, the bits of `weighted_ld` and `correlation` are the
  same with 1, 2, 3 and 8 threads (each run in its own rayon pool).
- `test_ld_decay_basic`: spacing 0.005, bin 1 cM, max 25 cM, min pairs 1:
  non-empty, every bin distance > 0 and n_pairs > 0.
- `test_mesh_matches_exact_when_snps_sit_on_mesh_cells`: spacing 1/128 M,
  bin width 100/128 cM, max 10 bins, qbin 2 vs 0: same length, same
  n_pairs, weighted_ld within 1e-12, correlation within 1e-9.
- Mesh quirk tests, spacing 1/128 M, bin width 100/128 cM, 130 bins,
  `min_pairs_per_bin` 0, qbin 1 (window m = 50, n = 128; real lags ≤ 49;
  first wrapped bin 79 = n − m + 1; unscaled lag/bin 127 = n − 1); "differ"
  means n_pairs differ, or weighted_ld by > 1e-12, or correlation by > 1e-9:
  - `test_corrected_mesh_matches_exact_past_the_fft_window`: corrected
    mesh = exact mode in every bin.
  - `test_v753_mesh_reproduces_dates_quirks`: V753 mesh vs exact: the
    first differing bin is 79, and bin 127 differs.
  - `test_circular_fft_lags_flag_only_changes_wrapped_bins`: only
    `circular_fft_lags` vs corrected: some bins differ, all ≥ 79.
  - `test_unscaled_last_lag_flag_only_changes_its_bin`: circular vs
    circular + unscaled: exactly bin 127 differs.
- `test_leave_one_out_with_single_chromosome_is_empty`: default options,
  one chromosome: one `ChromosomeLd` (chromosome 1, n_snps 50) whose bins
  all have n_pairs 0.
- `test_ld_decay_missing_population`: unknown test population → error.
- `a_snp_without_a_genetic_position_is_an_error`: SNP 7's position NaN →
  error whose message contains "rs7".

`exponential_fit.rs` (synthetic bins at 0.5 … 30 cM every 0.1 cM,
correlation A·r^i + K with r = exp(−n · 0.001), i from 0):
- `test_exponential_fit_recovery`: n = 50, A = 0.001, K = 1e-5, fit from
  0.45 to 30 cM: decay_rate within 1e-6 relative, amplitude within 1e-9,
  intercept within 1e-12, r² > 0.999999.
- `test_exponential_fit_different_decay_rates`: n ∈ {1, 20, 50, 100, 200}
  recovered within 1e-6 relative.
- `test_exponential_fit_insufficient_data`: one bin → error.
- `test_weighted_jackknife_equal_weights_is_ordinary_jackknife`: within
  1e-12 of the ordinary jackknife formulas (section 6).
- `test_weighted_jackknife_ignores_zero_weight_blocks`: zero-weight blocks
  give exactly the same result; one positive block → error.

`dating.rs` (20 "test", 10 "source1", 10 "source2" samples; 3 × 100 SNPs
on chr 1-3, positions i·0.003 M; patterned genotypes):
- `test_dates_basic`: `min_pairs_per_bin` 10, no jackknife: n_snps > 0,
  generations > 0, years_bp > 0, bins non-empty.
- `every_sample_diploid_flag_counts_pseudohaploid_calls_as_two_alleles`
  (DATES-0004 isolation): a source of mixed ploidy (after
  `detect_ploidy`); with only `every_sample_diploid` the weights
  (`valid_indices` and the weights' `Debug` text) and every bin's
  (n_pairs, weighted_ld, correlation) equal exactly those of the same
  data with `set_all_diploid()` and default flags, and differ from the
  default run's weights.
- `test_dates_missing_population`: unknown test population → error.

### Other tests

- `crates/popmix-core/tests/dates_v753_golden.rs`: feeds DATES's own
  `.out` (bins with distance, weighted_ld, correlation, n_pairs; other
  fields NaN) to `fit_exponential` (min 0.45 cM, max ∞, affine, bin
  0.1 cM) and checks the decay rate against DATES's full-data estimate
  within 5e-4 + 1e-3; feeds `.jin` (weight, estimate; zero-weight rows
  included) and the full estimate to `weighted_jackknife` and checks the
  `.jout` estimate and SE within 5e-4. Configs chr20_22_q0, chr20_22_q10,
  simulation_q0, simulation_q10, aadr_french_q10.
- `tests/golden/test_dates_v753.py`: runs `pm.dates` with each config's
  par values (bin width, maxdis, lovalfit, qbin) and checks, with
  `compat="dates-v753"`: the LD table (`ld_table`: n_pairs exact,
  weighted_ld, regression_coef, correlation to printed precision, bins
  matched by distance rounded to 3 decimals), the full fit
  (`fit.decay_rate`), the per-chromosome blocks (same chromosome set as
  the `.jin` rows with weight > 0, weights exactly equal, estimates within
  the fit tolerance), and the jackknife result (within (g + 2) · fit tol +
  5e-4). And for the default: `default_ld_table` (exact mode: bit-identical
  to compat; mesh: the bins that differ from compat, compared with exact
  equality of (n_pairs, weighted_ld, correlation), form a tail of the
  table, and every other bin matches DATES); `default_result` (where no
  bin differs, the jackknife result matches DATES); `default_vs_exact_mode`
  (corrected mesh pair counts within 1/qbin of DATES's exact-mode counts,
  first bin exempt). Quick configs chr20_22_q0, chr20_22_q10,
  yamnaya_chr21_22_q10 (mixed ploidy, DATES-0004); slow ones simulation_q0,
  simulation_q10, aadr_french_q10.
- `tests/golden/test_dates.py`: reference JSON format; `pm.dates`
  signature (`geno, test_pop, source1, source2, generation_time,
  min_distance_cm, max_distance_cm, bin_width_cm, run_jackknife`) and the
  `DatesResult`, `LdBin`, `DecayFit` classes; chr20-22 default run (first
  bin n_pairs > 100000, weighted_ld > 0); `test_dates_vs_reference`: the
  default (exact mode, maxdis 100 cM) table equals
  `chr20_22_reference.json` (pair counts exact, correlation within 1e-6).
- `tests/python/test_dates.py`: unknown compat rejected with a message
  naming "dates-v753"; with qbin 10 and no jackknife, compat changes some
  bins (by n_pairs or correlation), all beyond 90 cM.
- `tests/python/test_bad_options.py`: bin width 0, −0.1, NaN, 1e-12 →
  `ValueError` matching "bin_width_cm"; qbin 10⁹ → matching "qbin";
  generation_time 0 or NaN → matching "generation_time".
- `tests/python/test_cli.py` (`TestDatesCommand`): text output mentions
  generations and amplitude; JSON has generations, amplitude, ld_bins; TSV
  has a header; `--qbin 10 --compat dates-v753` reproduces chr20_22_q10
  (decay rate within 1.5e-3, jackknife estimate and SE within 0.008);
  missing input fails.
- `tests/python/test_rolloff.py`: ROLLOFF's jackknife via
  `weighted_jackknife`.
- `tests/nextflow/test_nfcore_modules.py::test_dates`: the nf-core module's
  JSON equals `popmix dates --format json` (determinism).
