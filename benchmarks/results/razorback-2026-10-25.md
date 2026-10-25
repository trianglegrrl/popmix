Measured 2026-10-25 on 13th Gen Intel(R) Core(TM) i7-1360P (16 threads, 62.4 GB); popmix 1.0.0 (5b27c90), admixtools2 2.0.10; median of 5 runs after a warm-up run.

| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix | Agreement | Same output at every thread count |
|---|---|---|---|---|---|---|---|
| convertf: packed to EIGENSTRAT text | AdmixTools 8.0.2 convertf | 78.86 s (CV 4.8%) | 3.51 s (1; CV 4.6%) | 22.5x | 1018 / 1044 MB | pass (snps 1233013, differing_snps 0) | yes |
| extract_f2: f2 blocks, 16 populations | admixtools2 2.0.10 f2_from_geno | 39.49 s (CV 1.6%) | 5.55 s (16; CV 8.8%) | 7.1x | 7885 / 759 MB | pass (blocks 710, max_abs 1.1e-16, max_rel 1.2e-14) | yes |
| qpadm: qpAdm from genotypes (2 sources, 7 right) | admixtools2 2.0.10 qpadm | 22.57 s (CV 4.5%) | 4.71 s (16; CV 7.7%) | 4.8x | 3731 / 631 MB | pass (weights_rel 5.3e-12, se_rel 1.2e-11, p_rel 9.1e-12) | yes |
| pca: PCA, 8 present-day populations, 10 PCs | EIGENSOFT 7.2.1 smartpca | 83.35 s (CV 4.8%) | 5.87 s (16; CV 12.5%) | 14.2x | 755 / 592 MB | pass (eigenvalues 10, max_rel 4.3e-12) | yes |
| kinship: KING-robust kinship, every pair of 8 populations' samples | PLINK 2 --make-king-table (16 threads) | 2.18 s (CV 36.1%) | 3.99 s (16; CV 6.1%) | 0.5x | 131 / 737 MB | pass (pairs 24090, max_abs 5.0e-13) | yes |
| qpgraph: qpGraph, 6 populations, 1 admixture, 10 starts | admixtools2 2.0.10 qpgraph | 18.24 s (CV 8.5%) | 3.65 s (16; CV 1.0%) | 5.0x | 2321 / 618 MB | pass (score 3.0e+02, score_rel 2.8e-12) | yes |
| qpdstat: 20 f4-statistics from genotypes | admixtools2 2.0.10 qpdstat | 17.19 s (CV 3.4%) | 4.04 s (16; CV 5.2%) | 4.3x | 670 / 341 MB | pass (statistics 20, est_rel 3.4e-11, se_rel 7.6e-12) | yes |

popmix by thread count (median wall time):

- convertf: 1: 3.51 s, 4: 3.79 s, 16: 3.99 s
- extract_f2: 1: 6.58 s, 4: 6.37 s, 16: 5.55 s
- qpadm: 1: 5.14 s, 4: 5.11 s, 16: 4.71 s
- pca: 1: 10.64 s, 4: 6.44 s, 16: 5.87 s
- kinship: 1: 5.75 s, 4: 4.25 s, 16: 3.99 s
- qpgraph: 1: 4.41 s, 4: 4.03 s, 16: 3.65 s
- qpdstat: 1: 4.14 s, 4: 4.33 s, 16: 4.04 s
