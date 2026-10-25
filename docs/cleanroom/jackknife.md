# Clean-room specification: the weighted block jackknife

`crates/popmix-core/src/fstats/jackknife.rs`. Written from Busing, Meijer
and van der Leeden 1999 (Statistics and Computing 9:3, "Delete-m jackknife
for unequal m"), as used for f-statistics by Patterson et al. 2012, and
from admixtools2's black-box outputs (every f-statistic golden under
`tests/golden/admixtools2/`, compared at rtol 1e-10 for estimates and 1e-6
for standard errors). No code.

## Notation

Blocks k = 1..g with lengths m_k (SNPs), n = Σ m_k, h_k = n / m_k. For a
statistic that is a weighted mean of per-block values x_k, the full-data
value is θ = Σ m_k x_k / n and the leave-one-block-out value is
θ_(k) = (θ − x_k m_k / n) / (1 − m_k / n).

The jackknife estimate (Busing et al., eq. 5):

  θ_J = g θ − Σ_k (1 − m_k / n) θ_(k),

and the variance (eq. 7), through the pseudo-values
τ_k = h_k θ − (h_k − 1) θ_(k):

  var = (1/g) Σ_k (τ_k − θ_J)² / (h_k − 1).

For a block mean θ_J = θ exactly; the formulas are written for statistics
that are not block means too, where θ (`tot`) is computed separately.

## Functions (keep the names and signatures; callers use them)

- `block_jackknife(values: &[f64], lengths: &[f64]) -> Option<(f64, f64)>`:
  estimate and standard error (sqrt of var) of a block mean. Blocks whose
  value is not finite are left out entirely (as if absent); `None` if none
  remains.
- `leave_one_out(values, lengths) -> Vec<f64>`: θ_(k) for every block; the
  total θ uses only finite blocks (n their total length); a non-finite
  block's own value stays non-finite.
- `leave_one_out_nafix(values, lengths) -> Vec<f64>`: as above but every
  block gets a value: θ from the finite blocks only (sum over them of
  m_k x_k divided by their total length), then each block's θ_(k) with its
  value taken as 0 if not finite and with m_k / n relative to the total
  length of **all** blocks.
- `jackknife_with_total(loo, lengths, tot) -> Option<(f64, f64)>`: θ_J and
  sqrt(var) from given leave-one-out values and a given θ = `tot`; blocks
  with a non-finite θ_(k) are left out; `None` if none remains.
- `block_jackknife_cov(stats: &[Vec<f64>], lengths) -> Option<(Vec<f64>, Array2<f64>)>`:
  several block-mean statistics at once (`stats[s][k]`): keep the blocks
  where every statistic is finite (`None` if fewer than two); for each
  statistic θ_(k) as above and θ_J; the covariance of statistics s and t is
  (1/g) Σ_k u_sk u_tk with u_sk = (τ_sk − θ_J,s) / sqrt(h_k − 1). Returns
  the θ_J and the covariance matrix.

Keep the existing tests in the file's test module; they are written against
these definitions.
