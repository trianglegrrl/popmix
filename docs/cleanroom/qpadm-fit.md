# Clean-room specification: the qpWave/qpAdm fit

This specification describes, in mathematical and behavioural terms, the
fit behind popmix's qpWave, qpAdm, `qpadm_p`, rotation and pairs
(`crates/popmix-core/src/qpadm/fit.rs`). It was written from the published
method (Patterson et al. 2012, Genetics, appendix; Haak et al. 2015, Nature,
SI 10; Harney et al. 2021, Genetics), admixtools2's documentation, and
numerical facts established by running admixtools2 2.0.10 as a black box
(its outputs are the golden references under `tests/golden/admixtools2/`).
It contains no code. The implementation is written from this document and
the tests alone, by someone who has not read admixtools2's source or the
earlier `fit.rs`.

## Inputs

- An f4 matrix X with m rows and n columns: row i, column j is
  f4(L0, L_{i+1}; R0, R_{j+1}) for the left populations L0..Lm and right
  populations R0..Rn (the callers build it; m >= 1, n >= 1).
- Q, the (m·n) x (m·n) jackknife covariance of the entries of X, indexed
  **row-major**: the entry of X at row i, column j is element i·n + j of
  vec(X). Every vectorisation below is row-major.
- A fudge factor f (default 1e-4).

## 1. Weight matrix

W = (Q + f·tr(Q)·I)⁻¹: add f times the trace of Q to every diagonal
element, then invert. A singular result is a `SingularMatrix` error
("f4 covariance matrix").

## 2. Rank-r fit

For a rank r with 1 <= r <= min(m, n), find A (m x r) and B (r x n)
minimising

  S(A, B) = (x − vec(A B))ᵀ W (x − vec(A B)),   x = vec(X).

If r exceeds the number of singular values of X (min(m, n)), it is a
`Validation` error whose message reads
"rank {r} exceeds the f4 matrix's {k} singular values".

Procedure (the numbers depend on it, so it is part of the behaviour):

1. Start: B = the r right singular vectors of X with the largest singular
   values (as the rows of B, largest first); A = X Bᵀ.
2. Then 20 rounds, each an A-step followed by a B-step:
   - A-step: with B fixed, S is a quadratic in the m·r entries of A. Form
     its normal equations M a = c, add f·tr(M) to every diagonal element of
     M, and solve. The entries of a are A's entries **row by row**.
   - B-step: with A fixed, likewise for the r·n entries of B (normal matrix
     M, add f·tr(M) to its diagonal, solve), B's entries row by row.
   A singular normal matrix is a `SingularMatrix` error ("qpWave/qpAdm fit").
3. Rank 0: A is the m x 1 zero matrix and B the 1 x n zero matrix (no
   iterations).

The fit is a value with fields `a` (A) and `b` (B), returned by
`fit_rank(x, w, r, f)`.

## 3. Fit statistic

chi² = (x − vec(A B))ᵀ W (x − vec(A B)), with (m − r)(n − r) degrees of
freedom (`evaluate_fit(x, w, fit, r)` returns both).

p-value: the upper tail of the chi-squared distribution with that many
degrees of freedom (`chisq_p_value`); NaN when the degrees of freedom are 0
(a saturated model is not tested). For nested comparisons a separate
`chisq_upper_tail(x, df)` is R's `pchisq(x, df, lower.tail = FALSE)`:
with df 0 it is 1 for x < 0 and 0 otherwise; for x <= 0 it is 1.

## 4. qpAdm weights

In qpAdm the target is L0 and the sources L1..Lm, and the model has rank
r = m − 1. The weights w (one per source, summing to 1) are those for which
the sources' rows of the fitted A combine to zero: wᵀA = 0 with Σw = 1.

Behaviourally: let C be the (r + 1) x m matrix whose first r rows are Aᵀ
and whose last row is all ones, and e the vector (0, ..., 0, 1) of length
r + 1. Unconstrained: w solves the normal equations (CᵀC) w = Cᵀe; then w
is divided by its sum. Constrained: w minimises ½ wᵀ(CᵀC)w − (Cᵀe)ᵀw
subject to Σw = 1 and 0 <= w_i <= 1 (a convex quadratic program; any exact
method), then divided by its sum; no feasible point is an `Infeasible`
error ("constrained qpAdm weights"). With rank 0 (one source) the weight
is [1].

`qpadm_weights(x, w, r, f, constrained)` fits rank r and returns the
weights; `weights_from_fit(fit, r, constrained)` does the weights step for
an existing fit; `simplex_qp(D, d)` is the quadratic program
(minimise ½ wᵀDw − dᵀw subject to Σw = 1, 0 <= w <= 1).

## 5. Weight standard errors

Given the leave-one-block-out f4 matrices X₁..X_N (the jackknife
replicates), compute the weights of each with the **same** W, then for each
source the standard error is sqrt(Σ_k (w_k − w̄)²) over the N replicates
(w̄ their mean). With N < 2 the standard errors are NaN.
(`qpadm_weight_std_errors(replicates, w, r, f, constrained)`.)

## 6. Models dropping sources (popdrop)

`drop_pops(x, w, f, constrained, compat)` fits every non-empty subset of
the m sources. Subsets are taken by size (1, 2, ..., m), each size in
lexicographic order of the source indices (R's `combn`). For a subset S:

- X_S: the rows of X in S; W_S: the rows and columns of the **full** W
  belonging to those rows (no new inversion);
- rank |S| − 1; fit, chi² and dof as above, and weights as in section 4;
- `pattern`: one character per source, '0' if in S, '1' if dropped;
  `wt` = (m − 1) − rank; `f4rank` = rank; `feasible` = every weight in
  [0, 1]; `weights`: per source, the weight or None if dropped;
  `p_value` from `chisq_p_value`.

Then the nested path is marked, and the rows are sorted by (dof, pattern).

### Nested path (popmix's corrected rule, DEVIATIONS.md ADMIXTOOLS2-0001)

Start at the full model (the subset of all sources), mark it `best`. Then
repeatedly drop one more source: among the current model's children (rank
one lower, every source dropped by the current model also dropped, and
feasible) take the one with the lowest chi² (the first in subset order on a
tie), mark it `best`, and record on it the comparison with the current
model: `dofdiff` = its dof − current dof, `chisqdiff` = its chi² − current
chi², `p_nested` = `chisq_upper_tail(chisqdiff, dofdiff)`. Stop at rank 0 or
when no child is feasible.

### admixtools2 2.0.10's path (`compat.popdrop_sibling_path`)

Mark every model of rank m − 2 as `best`; order them by decreasing dof
(stable, keeping subset order among equals); compare each with the next in
that order (the same three fields, on the earlier of the two). Nothing else
is marked. With m = 1 nothing is marked.

## 7. Helpers the callers use

- `to_dmatrix(&ndarray::Array2<f64>) -> nalgebra::DMatrix<f64>`.
- `fudged_inverse(cov, f)`: section 1.
- `is_model_failure(&PopmixError) -> bool`: true for `SingularMatrix` and
  `Infeasible`.
- `by_p_descending(a, b)`: ordering of p-values, largest first, NaN after
  every number.
- `PopDropRow`: public struct with fields `pattern: String, wt: usize,
  dof: usize, chi_sq: f64, p_value: f64, f4rank: usize,
  weights: Vec<Option<f64>>, feasible: bool, best: bool,
  dofdiff: Option<usize>, chisqdiff: Option<f64>, p_nested: Option<f64>`.

## 8. Signatures

The module's items, with these exact signatures (callers and tests use
them; `Result` is `crate::Result`, matrices are nalgebra's):

```text
pub(crate) fn to_dmatrix(a: &Array2<f64>) -> DMatrix<f64>
pub(crate) fn is_model_failure(e: &PopmixError) -> bool
pub(crate) fn fudged_inverse(cov: &Array2<f64>, fudge: f64) -> Result<DMatrix<f64>>
pub(crate) struct RankFit { pub a: DMatrix<f64>, pub b: DMatrix<f64> }
pub(crate) fn fit_rank(xmat: &DMatrix<f64>, qinv: &DMatrix<f64>, rnk: usize, fudge: f64) -> Result<RankFit>
pub(crate) fn evaluate_fit(xmat: &DMatrix<f64>, qinv: &DMatrix<f64>, fit: &RankFit, rnk: usize) -> (f64, usize)
pub(crate) fn qpadm_weights(xmat: &DMatrix<f64>, qinv: &DMatrix<f64>, rnk: usize, fudge: f64, constrained: bool) -> Result<Vec<f64>>
pub(crate) fn weights_from_fit(fit: &RankFit, rnk: usize, constrained: bool) -> Result<Vec<f64>>
pub struct PopDropRow { ... }    // derives Debug, Clone; fields in section 7
pub(crate) fn drop_pops(xmat: &DMatrix<f64>, qinv: &DMatrix<f64>, fudge: f64, constrained: bool, compat: Admixtools2Compat) -> Result<Vec<PopDropRow>>
fn mark_nested_path(rows: &mut [PopDropRow], full_rank: usize, compat: Admixtools2Compat)   // rows in subset order
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>>   // k-subsets of 0..n, lexicographic
fn chisq_upper_tail(x: f64, df: usize) -> f64
pub(crate) fn simplex_qp(d_mat: &DMatrix<f64>, d_vec: &DVector<f64>) -> Result<DVector<f64>>
pub(crate) fn qpadm_weight_std_errors(leave_one_out: &[Array2<f64>], qinv: &DMatrix<f64>, rnk: usize, fudge: f64, constrained: bool) -> Result<Vec<f64>>
pub(crate) fn by_p_descending(a: f64, b: f64) -> std::cmp::Ordering
pub(crate) fn chisq_p_value(chisq: f64, dof: usize) -> f64
```

`Admixtools2Compat` (in `super::compat`) has a `popdrop_sibling_path: bool`.
