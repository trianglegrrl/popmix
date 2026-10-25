# Clean-room specification: qpGraph's drift lengths at fixed admixture weights

`Model::edge_lengths` in `crates/popmix-core/src/qpgraph/qpgraph_fit.rs`.
Written from Patterson et al. 2012 (Genetics, qpGraph) and Maier et al. 2023
(eLife, admixtools2), admixtools2's documentation and its black-box outputs
(`tests/golden/admixtools2/v2.0.10/*/qpgraph.json`, compared at fixed
weights to rtol 1e-10). No code.

## Problem

At fixed admixture weights each f3-statistic the graph predicts is linear in
the drift lengths l (one per drift edge): fitted f3 = D l, where D (the
"design", f3 statistics x edges) is built by the caller. With y the
observed f3 estimates and P the inverse of their covariance (both in
`self.f3`: `est` and `ppinv`), the lengths minimise

  (y − D l)ᵀ P (y − D l) + λ ‖l‖²,   λ = `diag` · tr(Dᵀ P D) / k,

where k is the number of edges and `diag` the option (`self.diag`, default
1e-4): a ridge whose size is `diag` times the mean diagonal element of
DᵀPD. Unconstrained (`self.constrained` false), l solves the normal
equations (DᵀPD + λI) l = DᵀP y; a singular system is a `SingularMatrix`
error ("qpgraph edge-length system").

Constrained (the default): the same objective subject to 0 <= l_i <=
`EDGE_UPPER` for every edge (`EDGE_UPPER` is a constant in the module).
It is a strictly convex box-constrained quadratic program; solve it exactly
with the module's `box_qp(d, q, lo, hi)` (popmix's own active-set solver,
which minimises ½ zᵀdz − qᵀz on lo <= z <= hi). For numerical balance,
solve it in scaled variables z_i = s_i l_i with s_i = sqrt of the i-th
diagonal element of (DᵀPD + λI): the matrix becomes (DᵀPD + λI)_ij /
(s_i s_j), the linear term (DᵀP y)_i / s_i, the bounds 0 <= z_i <=
`EDGE_UPPER`·s_i, and l_i = z_i / s_i. (In exact arithmetic the scaling
changes nothing; the fixed-weight goldens at rtol 1e-10 are met with it.)

## Signature

```text
fn edge_lengths(&self, design: &DMatrix<f64>) -> Result<DVector<f64>>
```
