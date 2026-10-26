Measured 2026-10-25 on 13th Gen Intel(R) Core(TM) i7-1360P (16 threads, 62.4 GB); popmix 1.0.0 (4b093cd), admixtools2 2.0.10; median of 5 runs after a warm-up run.

| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix | Agreement | Same output at every thread count |
|---|---|---|---|---|---|---|---|
| convertf: packed to EIGENSTRAT text | AdmixTools 8.0.2 convertf | 77.96 s (CV 8.4%) | 3.22 s (1; CV 7.1%) | 24.2x | 1018 / 1044 MB | pass (snps 1233013, differing_snps 0) | yes |
| extract_f2: f2 blocks, 16 populations | admixtools2 2.0.10 f2_from_geno | 36.90 s (CV 1.0%) | 5.44 s (16; CV 10.6%) | 6.8x | 7885 / 759 MB | pass (blocks 710, max_abs 1.1e-16, max_rel 1.2e-14) | yes |
| qpadm: qpAdm from genotypes (2 sources, 7 right) | admixtools2 2.0.10 qpadm | 21.47 s (CV 2.9%) | 4.29 s (16; CV 6.1%) | 5.0x | 3731 / 630 MB | pass (weights_rel 5.3e-12, se_rel 2.0e-12, p_rel 2.0e-11) | yes |
| pca: PCA, 8 present-day populations, 10 PCs | EIGENSOFT 7.2.1 smartpca | 83.24 s (CV 8.5%) | 5.40 s (4; CV 10.7%) | 15.4x | 755 / 591 MB | pass (eigenvalues 10, max_rel 4.3e-12) | yes |
| kinship: KING-robust kinship, every pair of 8 populations' samples | PLINK 2 --make-king-table (16 threads) | 1.23 s (CV 17.4%) | 3.92 s (4; CV 3.7%) | 0.3x | 130 / 737 MB | pass (pairs 24090, max_abs 5.0e-13) | yes |
| qpgraph: qpGraph, 6 populations, 1 admixture, 10 starts | admixtools2 2.0.10 qpgraph | 16.78 s (CV 3.5%) | 3.66 s (4; CV 11.3%) | 4.6x | 2321 / 618 MB | pass (score 3.0e+02, score_rel 2.8e-12) | yes |
| qpdstat: 20 f4-statistics from genotypes | admixtools2 2.0.10 qpdstat | 16.75 s (CV 2.6%) | 3.92 s (16; CV 6.4%) | 4.3x | 670 / 341 MB | pass (statistics 20, est_rel 3.4e-11, se_rel 7.6e-12) | yes |
| admixture: ADMIXTURE K=3, 3 present-day populations, LD-pruned, seed 1 | ADMIXTURE 1.3.0 (1 thread) | 13.87 s (CV 2.7%) | 19.30 s (16; CV 7.1%) | 0.7x | 41 / 129 MB | pass (loglik -4.4e+06, reference_loglik -4.4e+06, loglik_diff 4.2e-03, seeds 3, popmix_spread 0.0e+00, reference_spread 7.9e-02) | yes |
| rolloff: ROLLOFF, French = Sardinian + Russia_Samara_EBA_Yamnaya, to 5 cM, chromosomes 21-22, jackknife over them | ADMIXTOOLS 7.0.2 rolloff (full-precision build) | 4.3 min (CV 0.9%) | 14.60 s (4; CV 3.7%) | 17.5x | 915 / 307 MB | pass (bins 50, zscore_rel 3.9e-12, regression_rel 4.1e-12, corr_rel 3.7e-12, bins_equal True) | yes |
| dates: DATES, French = Sardinian + Russia_Samara_EBA_Yamnaya, FFT mesh (qbin 10), to 50 cM, jackknife | DATES v753 (numbins fix, as the goldens) | 30.73 s (CV 1.3%) | 5.70 s (4; CV 2.7%) | 5.4x | 758 / 1125 MB | pass (bins_compared 20, pairs_equal True, max_cov_abs 4.6e-10, max_corr_abs 5.0e-07, fit_abs 3.2e-04, generations_abs 2.0e-03, se_abs 3.4e-04) | yes |
| read_kinship: READ kinship, every pair of Spain_EN samples | READv2 (73739be) | 28.26 s (CV 1.2%) | 1.89 s (16; CV 2.7%) | 14.9x | 950 / 508 MB | pass (pairs 253, classes_equal True, max_rel 4.8e-12) | yes |
| classic_qpdstat: 20 f4-statistics from genotypes | ADMIXTOOLS 7.0.2 qpDstat | 42.89 s (CV 1.5%) | 3.42 s (16; CV 9.6%) | 12.5x | 2169 / 341 MB | pass (statistics 20, est_abs 4.8e-07, z_abs 4.8e-04, snps_equal True) | yes |
| classic_qp3pop: 16 outgroup f3-statistics from genotypes | ADMIXTOOLS 7.0.2 qp3Pop | 2.4 min (CV 0.6%) | 3.34 s (16; CV 8.6%) | 43.8x | 755 / 341 MB | pass (statistics 16, est_abs 4.0e-10, se_abs 5.0e-10, z_abs 5.0e-04) | yes |
| classic_qpadm: qpAdm from genotypes (2 sources, 7 right) | ADMIXTOOLS 7.0.2 qpAdm | 36.61 s (CV 1.3%) | 3.67 s (16; CV 3.0%) | 10.0x | 1173 / 630 MB | FAIL (weights_abs 9.6e-04, se_abs 7.8e-04, chisq_abs 1.8e-01, p_rel 7.0e-02) | yes |
| classic_qpgraph: qpGraph, 6 populations, 1 admixture | ADMIXTOOLS 7.0.2 qpGraph | 90.98 s (CV 1.2%) | 3.05 s (16; CV 5.4%) | 29.8x | 760 / 618 MB | pass (score 3.1e+02, score_abs 3.9e-04) | yes |
| qpadm_rotate: qpAdm rotation, 5 populations between left and right, 4 fixed right | admixtools2 2.0.10 qpadm_rotate | 21.35 s (CV 4.8%) | 4.00 s (16; CV 1.2%) | 5.3x | 3765 / 633 MB | pass (models 30, rank_dof_feasible_equal True, chisq_rel 2.9e-11, p_rel 6.9e-10) | yes |
| find_graphs: find_graphs, 6 populations, 1 admixture, Chimp as outgroup: popmix exhaustive, admixtools2 evolutionary (20 generations, seed 1) | admixtools2 2.0.10 find_graphs | 45.24 s (CV 1.3%) | 4.40 s (16; CV 4.8%) | 10.3x | 2328 / 618 MB | pass (graphs_fitted 1965, score 3.0e+02, reference_score 3.0e+02, score_rel -2.8e-12) | yes |
| qpgraph_resample: qpGraph bootstrap, 10 resamples of the blocks, 6 populations, 1 admixture | admixtools2 2.0.10 qpgraph_resample_snps | 18.87 s (CV 4.8%) | 2.98 s (4; CV 8.7%) | 6.3x | 2328 / 616 MB | pass (replicates 10, same_blocks_drawn True, score_rel 4.1e-12) | yes |

popmix by thread count (median wall time):

- convertf: 1: 3.22 s, 4: 3.56 s, 16: 3.41 s
- extract_f2: 1: 6.27 s, 4: 5.56 s, 16: 5.44 s
- qpadm: 1: 4.66 s, 4: 4.42 s, 16: 4.29 s
- pca: 1: 10.04 s, 4: 5.40 s, 16: 7.06 s
- kinship: 1: 5.42 s, 4: 3.92 s, 16: 3.94 s
- qpgraph: 1: 4.08 s, 4: 3.66 s, 16: 3.75 s
- qpdstat: 1: 4.32 s, 4: 4.17 s, 16: 3.92 s
- admixture: 1: 51.24 s, 4: 22.43 s, 16: 19.30 s
- rolloff: 1: 24.68 s, 4: 14.60 s, 16: 16.04 s
- dates: 1: 7.92 s, 4: 5.70 s, 16: 5.72 s
- read_kinship: 1: 1.99 s, 4: 1.98 s, 16: 1.89 s
- classic_qpdstat: 1: 3.54 s, 4: 3.61 s, 16: 3.42 s
- classic_qp3pop: 1: 3.40 s, 4: 3.43 s, 16: 3.34 s
- classic_qpadm: 1: 3.87 s, 4: 3.68 s, 16: 3.67 s
- classic_qpgraph: 1: 3.11 s, 4: 3.09 s, 16: 3.05 s
- qpadm_rotate: 1: 4.30 s, 4: 4.08 s, 16: 4.00 s
- find_graphs: 1: 7.78 s, 4: 4.93 s, 16: 4.40 s
- qpgraph_resample: 1: 3.28 s, 4: 2.98 s, 16: 3.29 s
