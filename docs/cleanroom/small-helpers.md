# Clean-room specification: the remaining helpers (C2)

Five small pieces, each described by behaviour only. Implement each from
this document, the papers and documentation named, and black-box runs of the
reference tools. No code was consulted to write the implementations; none is
quoted here.

## 1. qpF4ratio's jackknife (`crates/popmix-core/src/qpadm/f4ratio.rs`)

Method: Reich et al. 2009 (Nature 461:489, supplement) and Patterson et al.
2012 (Genetics 192:1065): α = f4(A, O; X, B) / f4(A, O; C, B), a ratio of
means, jackknifed over blocks. Keep the public items (`F4RatioResult`,
`F4RatioOptions` with `min_denominator` default 1e-6, `f4ratio`) and their
errors.

- Per block k: numerator f4 n_k and denominator f4 d_k from the f2 blocks
  (the f4 identity in CLAUDE.md).
- Blocks with |d_k| < `min_denominator` are left out of everything (as if
  absent); so are blocks where either value is not finite.
- The estimate: θ = (Σ m_k n_k / Σ m_k) / (Σ m_k d_k / Σ m_k) over the kept
  blocks (m_k block lengths).
- The leave-one-out values: the ratio of the leave-one-out numerator and
  the leave-one-out denominator, each computed as in `leave_one_out` of
  `fstats::jackknife` on the kept blocks.
- α and its standard error: `jackknife_with_total(loo, lengths, θ)`; z =
  α / SE. `n_blocks` is the number of kept blocks. Fewer than two kept blocks
  is the error the current tests expect.

Goldens: `tests/golden/test_fstats_admixtools2.py` (f4ratio, rtol 1e-6),
plus the module's own tests.

## 2. Genotype f-statistics' jackknife (`crates/popmix-core/src/fstats/geno_fstats.rs`)

Replace the module's two private helpers (`est_to_loo`, `summarize`) with
calls to the shared jackknife (`fstats::jackknife`, rewritten from
`docs/cleanroom/jackknife.md`), adding a function there if one is missing.
Behaviour (statistics per block are means over the SNPs a statistic uses; a
block's weight m_k is the number of those SNPs, 0 when it uses none):

- leave-one-out values of a block mean: as `leave_one_out`, over the blocks
  with m_k > 0 (a block with m_k = 0 gets a non-finite value);
- for D and f3, which are ratios, the leave-one-out value is the ratio of
  the numerator's and the denominator's leave-one-out values;
- the total θ is the weighted mean of the finite leave-one-out values with
  weights 1 − m_k / n, n = Σ m_k over all blocks (for a block mean this is
  the full-data mean exactly; for a ratio it is not, and the goldens need
  this definition);
- estimate and SE: `jackknife_with_total(loo, m, θ)` over the blocks whose
  leave-one-out value is finite; z = est / SE; p two-sided normal.

Goldens: `tests/golden/test_fstats_admixtools2.py` (qpdstat/qp3pop from
genotypes) and the module's tests. Rewrite the module doc's reference to
admixtools2 internals as a description of the behaviour.

## 3. R's uniform generator (`crates/popmix-core/src/rcompat/rng.rs`)

Reproduces R's default generator (`RNGkind()`: "Mersenne-Twister",
"Inversion", "Rejection") so popmix's random starts match admixtools2's.
Sources: Matsumoto & Nishimura 1998 (ACM TOMACS 8:3) and their reference
implementation `mt19937ar.c` (2002, BSD-3: credit it in NOTICE); R's
documentation `?Random`, `?set.seed`, `?sample`, `?runif` and the R 3.6.0
NEWS entry on `sample.kind = "Rejection"`; R run as a black box (you may
print `.Random.seed` after `set.seed(s)` and compare draws). Do not read R's
C sources. Keep the API: `RRng::set_seed(i32)`, `unif_rand`, `runif(n)`,
`unif_index(n)`, `sample_replace(n, size)` (0-based indices).

Behaviour established from R 4.5.2:

- **Seeding.** Start from the seed as an unsigned 32-bit integer s. Scramble
  it: 50 times, s ← 69069·s + 1 (mod 2^32). Then fill R's 625 state words:
  for each, s ← 69069·s + 1 and the word is s. The first word is R's
  position counter, not part of the MT state: after seeding it is set to
  624 (so the first draw regenerates the 624-word state). The other 624
  words are the MT19937 state. (`.Random.seed` after `set.seed` shows all
  of this: its first element encodes the kind, the next 625 these words.)
- **A draw** (`unif_rand`): the next MT19937 tempered output y, as in
  `mt19937ar.c`'s `genrand_int32`, scaled to (0, 1) as y · 2.3283064365386963e-10
  (2^−32), then mapped into the open interval: values <= 0 become
  half of 1/(2^32 − 1), values >= 1 become 1 minus that. (The clamp's
  constants were known from R's source, read before the provenance rule;
  R's documentation says only that the values lie in the open interval.
  Only y = 0 reaches it, since y · 2^−32 < 1.)
- **`runif(n)`**: n draws.
- **`sample` with rejection** (`unif_index(n)`, uniform on 0..n−1): with
  b = ceil(log2(n)) bits, repeatedly form a b-bit integer from draws, 16 bits
  at a time (take floor(draw · 65536) as the next 16 bits, most significant
  first, stopping once at least b bits are collected, then keep the low b
  bits), until it is < n. n <= 0 is invalid.
- `sample_replace(n, size)`: `size` independent `unif_index(n)`.

Goldens: `crates/popmix-core/tests/rcompat_golden.rs` (runif and sample
against R, exact equality).

## 4. Python helpers

### `compare_fits` (`python/popmix/resample.py`)

Compares two vectors of scores from the same bootstrap replicates (Maier et
al. 2023, eLife, "comparing graphs"; admixtools2 `?compare_fits`). Keep its
signature and return keys. Differences d = scores1 − scores2, dropping
non-finite and zero differences; n the number left.

- `diff`: the mean of d; `se`: its standard deviation over replicates;
  `z` = diff / se; `p` the two-sided normal p of z (check the existing keys
  and the goldens for the exact set).
- `frac` values: the fractions of d below and above 0.
- `p_emp_nocorr` = 2 · min(those fractions);
- `p_emp` = 2 · min((count(d > 0) + 1) / (n + 1), (count(d < 0) + 1) / (n + 1)).

Goldens: `tests/golden/test_resample_admixtools2.py`.

### The model sets of `qpadm_models.py`

`rotate_models`, `qpadm_rotate` and `qpwave_pairs` produce model lists
whose order is observable in the goldens
(`tests/golden/test_qpadm_models_admixtools2.py`). Rewrite their
combinatorics in your own code:

- the subsets of a list in R's `power_set` order: by size, smallest first,
  each size in lexicographic order of positions (the empty set excluded);
- `rotate_models(pops, rightfix, ...)`: for every left set L (non-empty
  subset, at most min(|pops|, floor((|pops| + |rightfix| − 1) / 2)) members)
  the right set is the rest of `pops` followed by `rightfix`;
- `qpadm_rotate`: as above, but the right set is `rightfix` followed by the
  rest.

Confirm each order against the goldens.

### `fstat_combinations` (`python/popmix/dstats.py`)

For f3 the combinations of three populations, each in three rotations
(a, b, c), (b, c, a), (c, a, b): all first rotations, then all second, then
all third. For f4 the combinations of four, each in its three splits
(a, b; c, d), (a, c; b, d), (a, d; b, c), the three splits of one
combination together. The combinations themselves in lexicographic order.
Confirm against the goldens and the module's tests.

## 5. Exhaustive graph enumeration (`crates/popmix-core/src/qpgraph/graph_search.rs`)

`all_graphs(leaves, nadmix)`: every admixture graph on `leaves` with exactly
`nadmix` admixture nodes, each once, up to isomorphism as given by
`graph_key`, with the same **set** as admixtools2 2.0.10's
`generate_all_graphs` (its documentation: no graphs with fewer admixture
events, "most graphs which are unidentifiable" excluded, for example a node
expanded to a loop, no duplicates by `graph_hash`). The order is popmix's
own (make it deterministic); the golden test now compares sets.

- Graphs: rooted DAGs, a root of out-degree 2; leaves of in-degree 1 and
  out-degree 0; internal non-admixture nodes in-degree 1, out-degree 2;
  admixture nodes in-degree 2, out-degree 1. Built from rooted binary trees
  (counts (2n − 3)!! for n leaves: 1, 3, 15, 105, 945) by inserting
  admixture edges, each from a new node on one edge to a new admixture node
  on another, keeping the graph acyclic.
- Which candidate graphs admixtools2 excludes is not documented beyond the
  sentence above: find the rule from the golden set
  (`tests/golden/admixtools2/v2.0.10/dg_6pop/all_graphs.json`, and the
  counts 84 for 4 leaves and 1 admixture, 1485 for 5 and 1, 1302 for 4 and
  2), by comparing your candidate set with it. Describe the rule you find in
  the module doc, in terms of graph structure.
- `graph_key(edges)` stays as it is (its definition is in the golden
  generator's header, `generate_all_graphs.R`); so do the edge-list output
  (`edge_list`, base first) and `find_graphs`.

Goldens: `tests/golden/test_find_graphs_admixtools2.py` (sets of keys,
scores at fixed weights) and the module's tests (counts, base first,
internal labels).
