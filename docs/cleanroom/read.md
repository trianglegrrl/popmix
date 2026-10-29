# Clean-room specification: READ kinship (`crates/popmix-core/src/kinship/read.rs`)

This specification describes, in mathematical and behavioural terms, what
popmix's READ-style kinship estimator computes: the pairwise mismatch rate
P0, its windows and block jackknife, the normalisation, the classification
into degrees of relationship, the z scores, the parent-offspring / sibling
call, and the kinship coefficient with its interval, in popmix's corrected
form (`ReadCompat::CORRECTED`, the default) and in READv2 2.01's exact form
(`ReadCompat::V2_01`).

It was written from:

- Monroy Kuhn, Jakobsson and Günther 2018, "Estimating genetic kin
  relationships in prehistoric populations", PLOS ONE 13:e0195491,
  doi:10.1371/journal.pone.0195491 (READ: P0 per window, normalisation by
  the P0 of unrelated pairs, thresholds at the mid-points between the
  expected normalized P0 of each degree).
- Alaçamlı et al. 2024, "READv2: advanced and user-friendly detection of
  biological relatedness in archaeogenomics", Genome Biology 25:216,
  doi:10.1186/s13059-024-03350-3, Methods ("READ re-implementation"): the
  genome-wide P0 with a 5 Mb block jackknife, the third-degree class
  between 0.90625 and 0.953125 called only with at least 3000 "expected
  mismatches", the 20 Mb windows and the 0.3 / 0.35 / 0.6 proportions
  that tell parent-offspring from siblings with at least 10 000 expected
  mismatches, heterozygous sites as half a mismatch, and θ = 1 − normalized
  P0.
- READv2's README (https://github.com/GuntherLab/READv2, README.md only)
  and `READ2.py --help` / `-v` (version string `v2.01`): options and
  defaults (median normalisation, `--window_size` 5000000, `--window_est`,
  `--2pow`, `--norm_value`), the output columns.
- READv2 commit 73739bef (the one in
  `~/.local/share/popmix-references/readv2/COMMIT`) run as a black box on
  hand-made PLINK datasets (section 12 lists every run) and its outputs in
  the goldens `tests/golden/readv2/73739be/{family,aadr_gurgy}/read.json`.
  The model in this document reproduces both goldens: P0 and overlap SNPs
  exactly, the share of first-degree windows exactly, and the jackknife SE
  to 3e-14 relative (checked on every 37th pair of `family`, 85 sampled
  pairs of `aadr_gurgy`, and all 78 first-degree pairs of `aadr_gurgy`,
  whose types all match).

READv2 is GPL-3; this document contains no READv2 code and was written
without reading READv2's source, popmix's earlier `read.rs`, or
DEVIATIONS.md's READV2 entries. The implementer writes the module from this
document and the tests alone. Record the provenance (this document and the
papers above) in the module doc comment (CLAUDE.md, "Provenance").

## 1. Overview

For every pair of samples (i, j), over the SNPs called in both:

1. **P0**: the mean per-SNP mismatch (0, ½ or 1), either over all those SNPs
   at once (genome-wide, the default) or as the mean of per-window P0s
   (`window_estimate`).
2. **SE**: from windows of `window_size` base pairs: a jackknife of the
   mean of the window P0s (genome-wide mode), or the standard error of the
   mean of the window P0s (window mode).
3. **Normalisation**: one value for the whole run (the median, mean or
   maximum of P0 over the pairs, or a given value); normalized P0 =
   P0 / norm.
4. **Classification** of normalized P0 by four thresholds into identical,
   first, second, third degree or unrelated, third degree only with enough
   data.
5. **z scores**: distances, in standard errors, to the class's two
   boundaries.
6. **First-degree type**: from windows of `first_degree_window_size`
   (20 Mb), the share of windows whose normalized P0 is outside the
   first/second-degree range; parent-offspring, siblings or neither.
7. **Kinship** θ = 1 − normalized P0 and a 95 % interval.

## 2. Public API to preserve exactly

`kinship/mod.rs` (unchanged) declares `pub mod read;` and re-exports
`read_kinship, read_thresholds_2pow, FirstDegreeType, ReadCompat, ReadNorm,
ReadOptions, ReadPair, ReadRelationship, ReadResult, READ_THRESHOLDS`;
`crates/popmix-core/src/lib.rs` re-exports `read_kinship, ReadCompat,
ReadNorm, ReadOptions, ReadResult` at the crate root (unchanged). (The
listings here are the interface callers compile against, not
implementation.)

### Types (all fields `pub`)

- `ReadCompat` — derives `Debug, Clone, Copy, PartialEq, Eq, Default`.
  Four `bool` fields in this order: `windows_overlap`,
  `jackknife_skips_last_window`, `ci_of_normalized_p0`,
  `windows_ignore_chromosomes` (meanings in sections 4, 6 and 10).
  Associated constants `CORRECTED` (all four false; equal to `Default`) and
  `V2_01` (all four true).
- `ReadNorm` — derives `Debug, Clone, Copy, PartialEq`; variants `Median`,
  `Mean`, `Max`, `Value(f64)`.
- `ReadOptions` — derives `Debug, Clone`; fields `window_size: u64`,
  `first_degree_window_size: u64`, `norm: ReadNorm`,
  `window_estimate: bool`, `thresholds: [f64; 4]`, `compat: ReadCompat`.
  `impl Default`: `window_size` 5 000 000, `first_degree_window_size`
  20 000 000, `norm` `Median`, `window_estimate` false, `thresholds`
  `READ_THRESHOLDS`, `compat` `ReadCompat::CORRECTED` (READv2's documented
  defaults; the Python binding builds `ReadOptions::default()` and sets
  `window_size`, `norm`, `window_estimate`, `thresholds`, `compat`, never
  `first_degree_window_size`).
- `pub const READ_THRESHOLDS: [f64; 4] = [0.953125, 0.90625, 0.8125, 0.625];`
- `pub fn read_thresholds_2pow() -> [f64; 4]` — section 7.2.
- `ReadRelationship` — derives `Debug, Clone, Copy, PartialEq, Eq`;
  variants `Unrelated`, `UnrelatedOrThirdDegree`, `ThirdDegree`,
  `SecondDegree`, `FirstDegree`, `IdenticalOrTwins`; `pub fn label(self)
  -> &'static str` (section 7.3).
- `FirstDegreeType` — derives `Debug, Clone, Copy, PartialEq, Eq`;
  variants `ParentOffspring`, `Siblings`, `NotApplicable`; `pub fn
  label(self) -> &'static str` (section 9).
- `ReadPair` — derives `Debug, Clone, PartialEq`; fields in this order:
  `i: usize`, `j: usize`, `relationship: ReadRelationship`,
  `z_upper: Option<f64>`, `z_lower: Option<f64>`, `normalized_p0: f64`,
  `normalized_p0_se: f64`, `p0: f64`, `p0_se: f64`,
  `first_degree_type: Option<FirstDegreeType>`,
  `first_degree_windows_outside: f64`, `overlap_snps: u64`,
  `effective_snps: f64`, `kinship: f64`, `kinship_ci: (f64, f64)`.
- `ReadResult` — derives `Debug, Clone`; fields `norm_value: f64`,
  `pairs: Vec<ReadPair>`.

Keep the doc comments on these items as the declarations carry them
(`ReadCompat`'s fields, `ReadOptions`' fields, `ReadPair`'s fields); they
summarise the behaviour below.

### Function

```
pub fn read_kinship(g: &BitGenotypes, positions: &[u64], chroms: &[u32],
                    options: &ReadOptions) -> Result<ReadResult>
```

`Result` is `crate::error::Result`, errors `crate::error::PopmixError`.
`BitGenotypes` is `crate::kinship::king::BitGenotypes` (unchanged): for
sample s, `masks(s)` gives three bit sets (het, hom0, hom2) of
`n_snps().div_ceil(64)` words each; SNP k is bit k % 64 of word k / 64;
bits past `n_snps()` are zero. A SNP is **called** in s when its bit is set
in any of the three; otherwise it is missing.

### Private helpers the unit tests call

`read/tests.rs` (kept unchanged, `use super::*`) calls three private
functions, which must exist with these signatures:

- `windows(positions: &[u64], chroms: &[u32], size: u64, compat: ReadCompat)
  -> Vec<(usize, usize)>` — section 4.
- `jackknife_se(values: &[f64], skip_last: bool) -> f64` — section 6.1.
- `mean(values: &[f64]) -> f64` — the arithmetic mean, Σ in order / n.

Other private helpers are up to the implementer.

### Callers

- `crates/popmix-python/src/kinship.rs` (`read_kinship_py`): maps `norm`
  ("median"/"mean"/"max", or `norm_value` given → `Value(v)`), thresholds
  ("read" → `READ_THRESHOLDS`, "2pow" → `read_thresholds_2pow()`), compat
  (None → `CORRECTED`, "readv2-2.01" → `V2_01`); passes
  `positions` = each SNP's `position_bp` and `chroms` = a run index that
  starts at 0 and increases by one at each SNP whose chromosome differs from
  the previous SNP's (so only "differs from the previous SNP" is
  meaningful); reads every field of `ReadResult` / `ReadPair`, mapping
  `None` z scores to NaN and `None` type to "". After the call, for
  `Median` or `Mean` it warns when `pairs.len() < 10`
  (`warnings::check_read_pairs`; not this module's job).
- `python/popmix/kinship.py` and the CLI (`popmix kinship --method read`)
  call the binding.

## 3. Notation and per-SNP mismatch

n samples, m SNPs, SNP k at `positions[k]` = x_k (base pairs) and
`chroms[k]` = c_k. For a pair (i, j) and SNP k:

- called in both: o_k = 1 if k is called in i and in j, else 0;
- mismatch d_k (only where o_k = 1):
  - ½ if either sample is heterozygous at k (het-het, het-hom: whatever
    the other call);
  - 1 if one is hom0 and the other hom2;
  - 0 if both are hom0 or both hom2.
  Where o_k = 0, d_k = 0 and the SNP counts nowhere.

For an index range S of SNPs: D(S) = Σ_{k∈S} d_k, N(S) = Σ_{k∈S} o_k and
P0(S) = D(S) / N(S) (undefined when N(S) = 0).

Established by: unit test `identical_samples_are_identical_and_halves_are_counted`
(the heterozygote counts ½); paper Methods (all three het cases ½); run
e7 (het-hom, missing) and e11 (het-het) reproduce READv2's P0 to the last
digit with this rule; goldens (P0 exact on 85 + 32 checked pairs).

## 4. Windows (`windows`)

Input: positions, chroms (equal length L), a size w ≥ 1, the compat.

1. Bin of SNP k: b_k = ⌊x_k / w⌋ (integer division on the absolute
   position; not relative to the first SNP, not x − 1).
2. SNP order is the input order; SNPs are never sorted. A new window
   starts at SNP k > 0 when b_k ≠ b_{k−1}, or — unless
   `compat.windows_ignore_chromosomes` — when c_k ≠ c_{k−1}. (So bins that
   recur after a change start new windows: bins 0, 0, 1, 0 give three
   windows.) With starts s_0 = 0 < s_1 < … < s_{W−1} and e_r = s_{r+1} − 1
   (e_{W−1} = L − 1), the plain windows are [s_r, e_r].
3. With `compat.windows_overlap`, every window but the first also takes
   the last SNP of the window before it: window r ≥ 1 is
   [s_r − 1, e_r]. (Each SNP that ends a window, other than the last
   window, is then in two windows.)
4. Return the windows as inclusive index pairs `(start, end)`, in order.
   L = 0 gives an empty vector.

The same function makes both the jackknife windows (`window_size`) and the
first-degree windows (`first_degree_window_size`), with the same compat
flags.

Examples (from `read/tests.rs`): positions 1, 2, 6, 7, 11, one
chromosome, w = 5: CORRECTED `[(0,1),(2,3),(4,4)]`, V2_01
`[(0,1),(1,3),(3,4)]`; positions 1, 2, 3, 4 on chromosomes 1, 1, 2, 2,
w = 100: CORRECTED `[(0,1),(2,3)]`, V2_01 `[(0,3)]`.

On human data, positions restart near 0 on each chromosome, so the bin
changes at (almost) every chromosome boundary anyway; READv2's 5 Mb windows
on the AADR set number 604, its 20 Mb windows 164, both reproduced by this
rule with chromosomes ignored.

Established by: run e6 (bins ⌊x/w⌋ with SNPs at 99, 100, 101, 199, 200;
the chromosome case 10, 20 | 30, 40, 150 | 160, 50, 60 on chromosomes
1, 1, 2, 2, 2, 3, 3, 3 gives windows {0–3}, {4, 5}, {6, 7}, which only
"ignore chromosomes, split on bin change" explains; starting position
150 not used as origin); runs e4 and e5 (overlap: a single mismatch at the
last SNP of a window moves the SE, at any other SNP of the window it does
not; section 12 has the numbers); run e8 (overlap also in the 20 Mb
windows); goldens (window counts 28 / 10 on `family`, 604 / 164 on
`aadr_gurgy`).

## 5. P0 and overlap per pair

Let the jackknife windows be V = `windows(positions, chroms,
window_size, compat)`. For a pair, the **kept windows** are those r with
N(V_r) > 0 (a window where the pair has no SNP called in both is dropped
everywhere: from the mean, the jackknife and its g); p_r = P0(V_r) for the
kept windows, in order; g = number of kept windows.

- **Genome-wide mode** (`window_estimate` false, the default):
  - p0 = D(all SNPs) / N(all SNPs), each SNP counted once (windows play no
    part; overlap does not change it).
  - overlap_snps = N(all SNPs).
- **Window mode** (`window_estimate` true; READv1's method):
  - p0 = mean of p_r over the kept windows.
  - overlap_snps = Σ_r N(V_r) over all windows: with `windows_overlap` a
    SNP that is in two windows is counted twice (READv2 reports, e.g.,
    311 overlap SNPs for 300 SNPs in 12 windows); without it this is
    N(all SNPs).

A pair with N(all SNPs) = 0 (no SNP called in both) is not reported at all
and takes no part in the normalisation (section 7.1).

Established by: e7 (overlap = SNPs called in both; P0 the genome-wide
ratio; a window with nothing called in both dropped: its SE is that of the
remaining windows), e11 and e16w (window mode: p0 = mean of the kept
overlapping windows, overlap counted with duplicates, 311 for 300 SNPs and
62 for 53), e10 and e17b (a pair with no SNP in common is left out of the
output, with a warning, and the norm is computed without it).

## 6. Standard error of P0

### 6.1 `jackknife_se(values, skip_last)`

With x_1 … x_g the values (g = len):

- S = Σ_k x_k (summed in order), and for every k the leave-one-out mean
  ℓ_k = (S − x_k) / (g − 1).
- The set used: K = {1 … g} or, with `skip_last`, K = {1 … g − 1}.
- ℓ̄ = mean of ℓ_k over K; Q = Σ_{k∈K} (ℓ_k − ℓ̄)² (Q = 0 when K is empty).
- Return sqrt((g − 1) / g · Q). Note g, not |K|, in the factor, also with
  `skip_last`.
- g = 0 returns NaN. g = 1 returns NaN without `skip_last` (ℓ_1 = 0/0) and
  0 with it (K empty). g = 2 with `skip_last` returns 0.

Without `skip_last` this equals sd(x) / sqrt(g) (sd with g − 1 in the
denominator): the standard error of the mean of the values. Use this
order of operations: then equal values give exactly 0 (unit test
`jackknife_of_equal_windows_is_zero_and_skipping_one_changes_it`: [0.2,
0.2, 0.2] → 0.0 exactly; [0.1, 0.2, 0.4, 0.3] → within 1e-15 of
sd / 2, and `skip_last` changes it by more than 1e-3).

### 6.2 p0_se

- **Genome-wide mode**: p0_se = `jackknife_se(p, compat.jackknife_skips_last_window)`
  over the kept window P0s p_1 … p_g. This is a jackknife of the **mean of
  the window P0s** (unweighted), not of the genome-wide ratio, in both
  compat settings; READv2 does exactly this, and CORRECTED keeps it,
  changing only the windows and the left-out set. With V2_01 each window
  but the last is left out in turn (the last window's own P0 still enters
  every ℓ_k through S), and the factor is (g − 1)/g with g the number of
  kept windows including the last.
- **Window mode**: p0_se = `jackknife_se(p, false)` = sd(p) / sqrt(g), in
  both compat settings (READv2 leaves no window out of this one: the skip
  flag does not apply).
- normalized_p0_se = p0_se / norm_value.

Established by: runs e2 and e4 (a mismatch only in the last window gives
SE 0; e2a gives 0.07216878364870322 = 0.125 / √3, which fixes the factor
(g − 1)/g with g = 3 over the two left-out windows); e5 (windows of 2, 3,
4, 5 SNPs: single mismatches give 1/(6√2), 1/√96, 1/√288, 7/2400 under the
root, 1/√450 and 0 — exactly the window-mean jackknife and not a
leave-one-window-out of D/N, which section 12 shows does not fit); e3 and
e1 (random data, matched to 1e-15); e11, e16w (window mode = sd/√g over
the kept windows); e14 (one window: READv2 gives SE 0; two windows: 0);
goldens (SE to 3e-14 on both).

## 7. Normalisation, thresholds and classification

### 7.1 Norm value

Over the reported pairs (those with N(all SNPs) > 0), in the order of
section 11:

- `Median`: the median of their p0 (an even count averages the two middle
  values);
- `Mean`: their mean;
- `Max`: their maximum;
- `Value(v)`: v.

`norm_value` in the result is this value. Then for each pair
normalized_p0 = p0 / norm_value (one division; READv2's normalized P0 is
bit-for-bit P0 / norm in both goldens) and
effective_snps = overlap_snps (as f64) · norm_value (READv2's
"NSNPsXNorm", the paper's "expected mismatches").

Established by: e10 (median of 9 pairs, mean, max, with the no-overlap
pair excluded from each), e15 (median of 6 pairs = mean of the 3rd and
4th), `tests/python/test_kinship.py` (norm_value = median of p0;
Value(0.25) gives normalized = p0 / 0.25 exactly), goldens (norm_value =
median of the means file; normalized P0 = P0/norm exactly; NSNPsXNorm =
overlap · norm exactly).

### 7.2 Thresholds

`options.thresholds` = [t0, t1, t2, t3]: the normalized-P0 boundaries
unrelated / third, third / second, second / first and first / identical.

- `READ_THRESHOLDS` = [0.953125, 0.90625, 0.8125, 0.625] (the mid-points
  between the expected normalized P0 of 1, 0.9375, 0.875, 0.75 and 0.5).
- `read_thresholds_2pow()` (READv2 `--2pow`) = [1 − 2^−4.5, 1 − 2^−3.5,
  1 − 2^−2.5, 1 − 2^−1.5] ≈ [0.9558058261758408, 0.9116116523516815,
  0.8232233047033631, 0.6464466094067263] — the geometric mid-points of
  the kinship values 1/32 … 1/2 (computed as `1.0 - 2f64.powf(e)`).

Established by: paper Methods (READ thresholds), e9p (`--2pow`: the
boundaries recovered from READv2's z scores as x ± z·SE agree with these
four values to 1e-14 relative, and the class changes between the 2⁻⁶-step
values that bracket each).

### 7.3 Classification

With x = normalized_p0 and E = effective_snps:

| condition | relationship | `label()` |
|---|---|---|
| x > t0 | `Unrelated` | "Unrelated" |
| t1 ≤ x ≤ t0 and E ≥ 3000 | `ThirdDegree` | "Third Degree" |
| t1 ≤ x ≤ t0 and E < 3000 | `UnrelatedOrThirdDegree` | "Unrelated/Consistent with Third Degree" |
| t2 ≤ x < t1 | `SecondDegree` | "Second Degree" |
| t3 ≤ x < t2 | `FirstDegree` | "First Degree" |
| x < t3 | `IdenticalOrTwins` | "IdenticalTwins/SameIndividual" |

Note the inclusive t0 (x = 0.953125 exactly is third degree) and the
inclusive lower ends elsewhere. The labels are READv2's `Rel` strings
exactly (the golden compares them).

Established by: e9 (normalized P0 = k/64 exactly for k = 0 … 65 with
E = 3008: 0.609375 identical, 0.625 first, 0.796875 first, 0.8125 second,
0.890625 second, 0.90625 third, 0.953125 third, 0.96875 unrelated), e16
(E = 3000 exactly → third; 2999.5 → "Unrelated/Consistent with Third
Degree"), e7, e10 (small E → the "consistent" class), goldens (every
label; "consistent" pairs have E up to 2894, third-degree pairs from
3178).

## 8. z scores

se_n = normalized_p0_se. Each class has an upper boundary U (towards less
related) and a lower boundary B (towards more related):

| relationship | U | B |
|---|---|---|
| `Unrelated` | none | t1 |
| `UnrelatedOrThirdDegree` | t0 | t1 |
| `ThirdDegree` | t0 | t1 |
| `SecondDegree` | t1 | t2 |
| `FirstDegree` | t2 | t3 |
| `IdenticalOrTwins` | t3 | none |

z_upper = Some((U − x) / se_n), or None where U is none; z_lower =
Some((x − B) / se_n), or None where B is none. Note that unrelated pairs
measure z_lower to t1, not t0 (READv2 does, with or without enough data
for third degree). The divisions are plain IEEE: se_n = 0 gives ±inf or
NaN (0/0) inside `Some`, as READv2 writes inf / NaN; None is only for the
missing boundary (READv2's NA). With the 2pow thresholds the same table
applies with their values.

Established by: goldens (for every class, the z scores of all pairs equal
these formulas to 1e-9 with se_n from the means file; unrelated Zup and
identical Zdown are NA), e1 and e9 (unrelated with E < 3000: z_lower to
0.90625), e9 before the SE was nonzero (se_n = 0: inf, and NaN for x on
the boundary), e14 (one window: z_lower inf).

## 9. First-degree windows and type

Let F = `windows(positions, chroms, first_degree_window_size, compat)`
(the same compat flags; `window_size` plays no part: READv2's 20 Mb
windows do not change with `--window_size`). For the pair, keep the
windows with N(F_r) > 0, and for each y_r = P0(F_r) / norm_value. A window
is **inside** when t3 ≤ y_r ≤ t1 (note: t1 inclusive here, unlike the
genome-wide second-degree class), else **outside** (identical, third
degree or unrelated). With the 2pow thresholds, their t3 and t1.

first_degree_windows_outside = q = (number of outside kept windows) /
(number of kept windows). It is computed and reported for **every** pair,
whatever its class (READv2's `Perc_Win_1stdeg_P0`, a proportion in [0, 1]
despite its name). Every window is counted: the jackknife's skip-last does
not apply here.

first_degree_type:

- None unless relationship is `FirstDegree`;
- None if E < 10 000 (E = 10 000 is enough);
- Some(`ParentOffspring`) if q ≤ 0.3;
- Some(`Siblings`) if 0.35 ≤ q < 0.6;
- Some(`NotApplicable`) if q ≥ 0.6;
- None if 0.3 < q < 0.35.

`label()`: "Parent-offspring", "Siblings", "NotApplicable". (READv2 writes
"N/A" for every None case — the golden reads it as missing and the Python
binding writes "" — and the literal "NotApplicable" for q ≥ 0.6.)

Established by: paper Methods (20 Mb windows, proportion of windows
identical / unrelated / third degree, < 0.3, 0.35–0.6, 10 000 expected
mismatches); e12a (20 windows of 1000 SNPs, E = 10 000: q = 0.25 and 0.3
→ Parent-offspring, 0.35 and 0.5 → Siblings, 0.6 and 0.65 →
NotApplicable; E = 9999 and 9999.5 → N/A whatever q); e12b (40 windows:
0.275, 0.3 → Parent-offspring, 0.325 → N/A, 0.35, 0.575 → Siblings, 0.6 →
NotApplicable); e13 (a window at y = 0.625 and at y = 0.90625 is inside,
0.59375 and 0.9375 outside; with `--2pow`, 0.625 outside and 0.65625 and
0.90625 inside; a window with nothing called in both is dropped from the
denominator: 2/3, while one rescued by its overlap SNP counts: 3/4); e8
(overlap in the 20 Mb windows; all windows counted, none skipped); the
`aadr_gurgy` golden (164 windows, denominators 164, 163, 162 … where a pair
has empty windows; all 78 first-degree types reproduced; q reproduced
exactly), the `family` golden (q reported for unrelated, second- and
third-degree pairs).

## 10. Kinship and its interval

- kinship = 1 − normalized_p0 (READv2's θ; note it is 1 − normalized P0,
  not twice that).
- The 95 % interval uses the constant 1.96 and se_n:
  - CORRECTED (`ci_of_normalized_p0` false): (kinship − 1.96 · se_n,
    kinship + 1.96 · se_n), an interval around the kinship coefficient.
  - V2_01 (`ci_of_normalized_p0` true): (normalized_p0 − 1.96 · se_n,
    normalized_p0 + 1.96 · se_n): READv2 labels it the kinship
    coefficient's interval but centres it on the normalized P0.

Established by: goldens (CI = x ∓ 1.96 · se_n exactly, 0 difference, on
all pairs of both; kinship = 1 − x exactly), unit test
`kinship_interval_is_around_the_kinship_unless_compat`, golden test
`test_corrected_changes_only_the_errors` (CORRECTED: p0 and kinship as
V2_01 to 1e-14, some p0_se different, kinship strictly inside the
interval).

## 11. Pairs reported, order, and the overall procedure

`read_kinship`:

1. Validate (section 13).
2. V = `windows(positions, chroms, window_size, compat)`,
   F = `windows(positions, chroms, first_degree_window_size, compat)`.
3. For every pair (i, j), i < j, in the order (0, 1), (0, 2), …,
   (0, n−1), (1, 2), …, (n−2, n−1): if N(all SNPs) = 0 skip it; otherwise
   compute D and N over all SNPs, D and N per window of V and of F, then
   p0, p0_se, overlap_snps (sections 5–6). Pairs are independent: compute
   them in parallel if wished, but collect them in this order, and keep
   every per-pair sum in a fixed order so that results do not depend on
   the thread count (`tests/python/test_thread_determinism.py` compares
   bytes across thread counts).
4. norm_value from the collected pairs (7.1); check it (13).
5. For each pair: normalized_p0, normalized_p0_se, effective_snps,
   relationship, z scores, q, first_degree_type, kinship, kinship_ci.
6. Return `ReadResult { norm_value, pairs }`, the pairs in the order of
   step 3.

READv2 itself writes its pairs labelled "id1,id2" with id1 the sample
earlier in the .fam file, sorted as strings ("a,Ca" < "a,c" < "b,Ca" …);
popmix returns sample order and the golden test re-indexes by the label,
so this order is not reproduced (e10).

Samples with no called SNP at all simply give no pairs. Missing genotypes
never count as mismatches or as SNPs; they only remove the SNP from that
pair (and from that pair's windows).

READv2 also prints a warning for pairs with heterozygous sites and for
completely identical pairs; popmix has no such warnings in this module.

Performance: the AADR golden has 90 samples, 1.2 million SNPs and 4005
pairs. Work from the bit sets: for a pair, per 64-SNP word, called-in-both
= (het_i|hom0_i|hom2_i) & (het_j|hom0_j|hom2_j), full mismatches =
(hom0_i & hom2_j) | (hom2_i & hom0_j), half mismatches = called-in-both &
(het_i | het_j); window sums from popcounts of these words masked to the
window's bit range. D is a multiple of ½ and N an integer, so keep them as
integer counts (2D) and divide once: P0 per window and genome-wide are then
exact ratios, independent of summation order.

## 12. Black-box runs

Scripts and READv2 outputs: `/tmp/claude-1000/-home-a-popmix/e86dce51-440b-4ab8-83f2-5e3779939540/scratchpad/c9/spec/`
(`scripts/` the generators, `rb.py` writes PLINK .bed/.bim/.fam and runs
`READ2.py -i data` in `runs/<name>/`; `hyp.py` the model). All pseudo-
haploid calls are 0 or 2; "norm 0.5" means `-n value --norm_value 0.5`.
Recreate them as tests: each line below is a small fixture.

- **e1** — 200 SNPs at 1 Mb + 1 spacing, 3 random samples, defaults.
  Pairs: P0 0.46 / 0.52 / 0.55, SE 0.03454620077970638 /
  0.025480982009181707 / 0.03360398558469755 (V2_01 model: 6.3e-16
  relative). Second degree has z 0.3256508601839845 / 1.085502867279946.
- **e2** — norm 0.5, w = 10, SNPs at 1–4, 11–13, 21–25; sample 0 all 0.
  e2a mismatches at 1, 11, 12, 25: P0 1/3, SE 0.07216878364870322; e2b at
  11, 12, 25: P0 0.25, SE 0.14433756729740643; e2c at 4, 11, 12, 25: SE
  0.14433756729740643; e2d at 4 only: SE 0. (A leave-one-window-out of the
  genome-wide ratio gives 0.0625 for e2a, not 0.0722.)
- **e3** — windows of 2, 3, 4, 5, 6 SNPs (w = 100), random mismatches:
  SE 0.1097345281425435, 0.05809475019311122, 0.1355544171172596; window
  mean jackknife, overlap, skip last: equal to 1e-15. Every other
  combination (no overlap, no skip, D/N jackknife) is off by 5–40 %.
- **e4** — 4 windows × 3 SNPs (1–3, 101–103, 201–203, 301–303), w = 100,
  norm 0.5, one mismatch at each SNP in turn: SE 0.07856742013183861
  (SNPs 0, 1), 0.07081971546656644 (SNP 2, the last of window 0),
  0.05892556509887896 (SNPs 3–8), 0 (SNPs 9–11, the last window).
- **e5** — windows of 2, 3, 4, 5 SNPs, one mismatch in turn: 0.11785113…
  (SNP 0), 0.10206207261596575 (SNP 1, end of window 0),
  0.05892556509887896 (2, 3), 0.05400617248673217 (4, end of window 1),
  0.04714045207910317 (5–8: the end of window 2 is in window 3 = the last,
  which is never left out, so it is not special), 0 (9–13).
- **e6** — single-mismatch probes (one run, sample k+1 differs from sample
  0 at SNP k), w = 100: positions 99, 100, 101, 150, 199, 200, 250, 301;
  the chromosome case above; positions starting at 150. Only bins ⌊x/w⌋
  fit all; only "ignore chromosomes" fits the chromosome case.
- **e7** — 15 SNPs in 5 windows, 4 samples with missing calls and
  heterozygotes; a pair with every SNP of a window (including its overlap
  SNP) missing: SE 0.1571348402636772 = the jackknife over the 4 kept
  windows.
- **e8** — 20 Mb windows with 2 SNPs each (pos 1, 2, 20M+1, …): a mismatch
  at the last SNP of window 0 makes window 1 first-degree (q = 2/3): overlap
  applies; three windows, all counted.
- **e9 / e9p** — 6016 SNPs at 16 kb spacing, norm 0.5, 47·k mismatches so
  that x = k/64 exactly, E = 3008; e9p the same with `--2pow`.
- **e10** — 5 samples with .fam ids zed, b, a, Ca, c; Ca and c never both
  called. Pair labels "zed,b", "b,a" (fam order), sorted as strings;
  (Ca, c) absent; median 0.495 of 9 pairs, mean 4.44/9, max 0.535.
- **e11** — `--window_est`, 300 SNPs over 60 Mb, het-het at one SNP:
  P0 = mean of the 12 overlapping windows, SE = sd/√12, overlap 311 / 310.
- **e12a / e12b** — first-degree type grid (section 9).
- **e13** — first-degree window boundaries (section 9).
- **e14** — one window: SE 0, z_lower inf, CI (x, x); two windows: SE 0;
  `--norm_value 0` and `-0.5`: READv2 stops ("no valid value provided").
- **e15** — 6 pairs: median of even count; three identical samples: READv2
  stops ("Normalization value is zero").
- **e16** — E = 3000 → Third Degree, 2999.5 → consistent; e16w window mode
  with a fully missing window and duplicate-counted overlap (62 for 53).
- **e17** — one sample: READv2 stops ("At least two individuals are
  required"); two samples never both called, norm value given: empty
  output.
- **Goldens** — the model reproduces `family` (49 samples, 1176 pairs) and
  `aadr_gurgy` (90 samples, 4005 pairs) as stated at the top
  (scripts `fam.py`, `aadr.py`).

CORRECTED cannot be run in READv2: its values follow from the formulas
with the four flags false. Two derived examples for unit tests (from the
model, not from READv2): e2a with CORRECTED has windows {1–4}, {11–13},
{21–25}, window P0s 1/4, 2/3, 1/5 and SE = `jackknife_se([0.25, 2/3,
0.2], false)` = sd/√3 = 0.14792807728549273; e4 with a mismatch only in the last
window has CORRECTED SE = sd([0, 0, 0, 1/3])/2 = 1/12.

## 13. Errors

Checked before any pair is computed, in this order:

- `positions.len()` or `chroms.len()` ≠ `g.n_snps()`:
  `DimensionMismatch { context: "READ kinship positions" (or "chroms"),
  expected: n_snps, actual: len }`.
- `window_size` = 0 or `first_degree_window_size` = 0:
  `InvalidParameter { name: "window_size" (or
  "first_degree_window_size"), details: "must be at least 1 base pair" }`.
- `norm` = `Value(v)` with v not finite or v ≤ 0:
  `InvalidParameter { name: "norm_value", details: "must be a positive
  number, got {v}" }` (READv2 refuses 0 and negative values, e14).
- thresholds not all finite or not strictly decreasing (t0 > t1 > t2 >
  t3): `InvalidParameter { name: "thresholds", … }`.
- fewer than two samples: `TooFew { context: "READ kinship", what:
  "samples", required: 2, actual: n }` (READv2 stops, e17).

After the pairs:

- `Median`, `Mean` or `Max` with no reported pair (no two samples share a
  called SNP): `InvalidData { details: "READ: no pair of samples has a SNP
  called in both, so there is no P0 to normalize by" }`. With `Value(v)`
  the result is `ReadResult { norm_value: v, pairs: [] }` (READv2 writes
  an empty table, e17b).
- computed norm_value = 0 (or not finite): `InvalidData { details:
  "READ: the normalization value is {v}: are most samples duplicates?" }`
  (READv2 stops, e15).

No other input is an error: SNPs in any order, any chromosome codes, and
samples with no calls are all accepted.

## 14. What could not be pinned down by black-box runs

- CORRECTED's numbers are popmix's own (no reference produces them); they
  are fixed by sections 4–10 with the flags false, and the tests check
  only their properties (unchanged P0 and kinship, different SE, interval
  around the kinship).
- The exact floating-point order of READv2's jackknife sums: the model
  agrees to 3e-14 relative, so popmix's SE differs from READv2's in the
  last digits; the golden tolerance (1e-10) covers this. P0, overlap,
  q, labels and norm agree exactly.
- Whether READv2 compares q with ≤ 0.3 or computes 1 − q ≥ 0.7 (and
  similarly at 0.35 and 0.6) cannot be told apart at proportions that are
  exact in binary; q = count / total compared as in section 9 matches
  every run and golden.
- Window mode's first-degree share and third-degree / 10 000 rules were
  checked only through E (duplicate-counted overlap in V2_01); no run
  probed a window-mode first-degree pair.
