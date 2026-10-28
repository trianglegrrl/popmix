Measured 2026-10-28 on AMD EPYC 7B13 (16 threads, 62.8 GB); popmix 1.0.0 (c9f07b8), admixtools2 unavailable; median of 1 runs after a warm-up run.

| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix | Agreement | Same output at every thread count |
|---|---|---|---|---|---|---|---|
| convertf_aadr: AADR v66.p1 TGENO to PACKEDANCESTRYMAP, 23,089 x 1.23M | AdmixTools 8.0.2 convertf | 41.7 min (CV 0.0%) | 44.83 s (16; CV 0.0%) | 55.8x | 7597 / 529 MB | pass (bytes 7118189822, differing_files []) | yes |

popmix by thread count (median wall time):

- convertf_aadr: 1: 79.00 s, 4: 53.26 s, 16: 44.83 s
