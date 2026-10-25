# Clean-room specification: DATES's mesh of genetic distance

Two helpers of `crates/popmix-core/src/dates/ld_decay.rs`, used when
`qbin > 0` (the FFT path of the weighted LD decay): placing SNPs on a mesh
of genetic distance, and turning per-lag sums into distance bins. Written
from Narasimhan et al. 2019 (Science, supplementary methods of DATES) and
the outputs of DATES v753 run as a black box (`tests/golden/dates/v753/`,
and the Rust tests in `crates/popmix-core/tests/dates_v753_golden.rs`).
DATES's source carries no licence; the implementation is written from
this document alone. No code.

## Mesh cells (replaces `mesh_tags` and `LagSums`'s naming)

`mesh_cells(snps: &[Snp], binsize: f64, qbin: usize) -> Vec<usize>`: the
mesh cell of every SNP, in the order given. The SNPs are laid end to end
along one genome-wide axis: start at 0 at the first SNP; moving to the next
SNP adds the difference of their genetic positions (`position_morgans`)
when both are on the same chromosome, and adds exactly 5 (Morgans) — not
the position difference — when the chromosome changes. A SNP's cell is
floor(axis position / c), with c = binsize / qbin the cell width, computed
in floating point in that order.

## Lag sums to bins (replaces `LagSums` and `to_bins`)

For one chromosome, each mesh cell holds three numbers: the number of SNPs
in it, the sum of their values y, and the sum of y². For a lag d (cells
apart, d = 0, ..., D), the six lag sums are, over cell pairs (i, i + d):

| name (choose your own) | sum over cell pairs (i, i + d) |
|---|---|
| pairs | count(i) · count(i + d) |
| later | count(i) · Σy(i + d) |
| earlier | Σy(i) · count(i + d) |
| cross | Σy(i) · Σy(i + d) |
| later_sq | count(i) · Σy²(i + d) |
| earlier_sq | Σy²(i) · count(i + d) |

(They are computed by the FFT code that already exists in the module; only
the container and its conversion to bins are rewritten.) The container must
support `new(D)` (all zero, lags 0..=D), merging another (adding
element-wise), and conversion to distance bins:

`to_bins(binsize, qbin, n_bins) -> Vec<CorrBin>`: for every lag d >= 1
(lag 0, pairs within one cell, is not used) whose pair count is at least
0.5 (the FFT gives counts with rounding error; 0.5 means at least one
pair): its distance is d · (binsize / qbin), and its bin is that distance
divided by binsize, truncated to an integer, computed in floating point in
that order; lags whose bin is >= n_bins are left out. Into the bin's
`CorrBin` add: pairs to `s0`, later to `s1`, earlier to `s2`, cross to
`s12`, later_sq to `s11`, earlier_sq to `s22`. (`CorrBin`'s fields are the
pair count and the sums of y1, y2, y1², y1·y2 and y2² over pairs, y1 being
the later SNP of a pair.)

The FFT routine that fills the container writes these six arrays; rename
its uses to the new names. Results must not change: the DATES goldens
compare bins exactly (pair counts) and at 1e-6 (correlations).
