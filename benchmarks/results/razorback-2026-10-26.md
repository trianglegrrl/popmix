Measured 2026-10-26 on 13th Gen Intel(R) Core(TM) i7-1360P (16 threads, 62.4 GB); popmix 1.0.0 (4b3962a), admixtools2 2.0.10; median of 5 runs after a warm-up run. Reference times are from the run of 2026-10-25 (the reference tools are unchanged; agreement is checked against that run's outputs); popmix was run again.

| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix | Agreement | Same output at every thread count |
|---|---|---|---|---|---|---|---|
| convertf: packed to EIGENSTRAT text | AdmixTools 8.0.2 convertf | 77.96 s (CV 8.4%) | 1.83 s (1; CV 11.3%) | 42.5x | 1018 / 287 MB | pass (snps 1233013, differing_snps 0) | yes |
| extract_f2: f2 blocks, 16 populations | admixtools2 2.0.10 f2_from_geno | 36.90 s (CV 1.0%) | 1.01 s (16; CV 8.5%) | 36.6x | 7885 / 726 MB | pass (blocks 710, max_abs 2.4e-15, differing 0) | yes |
| qpadm: qpAdm from genotypes (2 sources, 7 right) | admixtools2 2.0.10 qpadm | 21.47 s (CV 2.9%) | 1.26 s (4; CV 3.3%) | 17.0x | 3731 / 623 MB | pass (weights_rel 1.5e-15, se_rel 2.0e-12, p_rel 4.0e-12) | yes |
| pca: PCA, 8 present-day populations, 10 PCs | EIGENSOFT 7.2.1 smartpca | 83.24 s (CV 8.5%) | 4.25 s (4; CV 7.3%) | 19.6x | 755 / 633 MB | pass (eigenvalues 10, max_rel 4.3e-12) | yes |
| kinship: KING-robust kinship, every pair of 8 populations' samples | PLINK 2 --make-king-table (16 threads) | 1.23 s (CV 17.4%) | 2.12 s (16; CV 23.1%) | 0.6x | 130 / 452 MB | pass (pairs 24090, max_abs 5.0e-13) | yes |
| qpgraph: qpGraph, 6 populations, 1 admixture, 10 starts | admixtools2 2.0.10 qpgraph | 16.78 s (CV 3.5%) | 1.07 s (16; CV 18.2%) | 15.7x | 2321 / 630 MB | pass (score 3.0e+02, score_rel 1.4e-12, seeds 3, popmix_spread 4.0e-10, reference_spread 0.0e+00) | yes |
| qpdstat: 20 f4-statistics from genotypes | admixtools2 2.0.10 qpdstat | 16.75 s (CV 2.6%) | 1.98 s (16; CV 17.3%) | 8.5x | 670 / 491 MB | pass (statistics 20, est_rel 3.4e-11, se_rel 7.6e-12) | yes |
| admixture: ADMIXTURE K=3, 3 present-day populations, LD-pruned, seed 1 | ADMIXTURE 1.3.0 (1 thread) | 13.87 s (CV 2.7%) | 6.77 s (4; CV 2.2%) | 2.1x | 41 / 137 MB | pass (loglik -4.4e+06, reference_loglik -4.4e+06, loglik_diff 4.2e-03, seeds 3, popmix_spread 0.0e+00, reference_spread 7.9e-02) | yes |
| rolloff: ROLLOFF, French = Sardinian + Russia_Samara_EBA_Yamnaya, to 5 cM, chromosomes 21-22, jackknife over them | ADMIXTOOLS 7.0.2 rolloff (full-precision build) | 4.3 min (CV 0.9%) | 16.11 s (16; CV 5.4%) | 15.9x | 915 / 484 MB | pass (bins 50, zscore_rel 3.9e-12, regression_rel 4.1e-12, corr_rel 3.7e-12, bins_equal True) | yes |
| dates: DATES, French = Sardinian + Russia_Samara_EBA_Yamnaya, FFT mesh (qbin 10), to 50 cM, jackknife | DATES v753 (numbins fix, as the goldens) | 30.73 s (CV 1.3%) | 2.56 s (4; CV 14.2%) | 12.0x | 758 / 604 MB | pass (bins_compared 20, pairs_equal True, max_cov_abs 4.6e-10, max_corr_abs 5.0e-07, fit_abs 3.2e-04, generations_abs 2.0e-03, se_abs 3.4e-04) | yes |
| read_kinship: READ kinship, every pair of Spain_EN samples | READv2 (73739be) | 28.26 s (CV 1.2%) | 1.05 s (4; CV 12.7%) | 26.8x | 950 / 394 MB | pass (pairs 253, classes_equal True, max_rel 4.8e-12) | yes |
| classic_qpdstat: 20 f4-statistics from genotypes | ADMIXTOOLS 7.0.2 qpDstat | 42.89 s (CV 1.5%) | 1.84 s (16; CV 9.4%) | 23.3x | 2169 / 426 MB | pass (statistics 20, est_abs 4.8e-07, z_abs 4.8e-04, snps_equal True) | yes |
| classic_qp3pop: 16 outgroup f3-statistics from genotypes | ADMIXTOOLS 7.0.2 qp3Pop | 2.4 min (CV 0.6%) | 1.61 s (16; CV 12.6%) | 90.9x | 755 / 422 MB | pass (statistics 16, est_abs 4.0e-10, se_abs 5.0e-10, z_abs 5.0e-04) | yes |
| classic_qpadm: qpAdm from genotypes (2 sources, 7 right) | ADMIXTOOLS 7.0.2 qpAdm | 36.61 s (CV 1.3%) | 0.95 s (16; CV 8.5%) | 38.7x | 1173 / 628 MB | FAIL (weights_abs 9.6e-04, se_abs 7.8e-04, chisq_abs 1.8e-01, p_rel 7.0e-02) | yes |
| classic_qpgraph: qpGraph, 6 populations, 1 admixture | ADMIXTOOLS 7.0.2 qpGraph | 90.98 s (CV 1.2%) | 0.90 s (16; CV 9.6%) | 101.0x | 760 / 617 MB | pass (score 3.1e+02, score_abs 3.9e-04) | yes |
| qpadm_rotate: qpAdm rotation, 5 populations between left and right, 4 fixed right | admixtools2 2.0.10 qpadm_rotate | 21.35 s (CV 4.8%) | 1.47 s (16; CV 16.9%) | 14.6x | 3765 / 641 MB | pass (models 30, rank_dof_feasible_equal True, chisq_rel 7.0e-12, p_rel 9.8e-10) | yes |
| find_graphs: find_graphs, 6 populations, 1 admixture, Chimp as outgroup: popmix exhaustive, admixtools2 evolutionary (20 generations, seed 1) | admixtools2 2.0.10 find_graphs | 45.24 s (CV 1.3%) | 2.37 s (16; CV 30.2%) | 19.1x | 2328 / 630 MB | pass (graphs_fitted 1965, score 3.0e+02, reference_score 3.0e+02, score_rel -1.4e-12, seeds 3, reference_spread 8.0e-12) | yes |
| qpgraph_resample: qpGraph bootstrap, 10 resamples of the blocks, 6 populations, 1 admixture | admixtools2 2.0.10 qpgraph_resample_snps | 18.87 s (CV 4.8%) | 1.51 s (16; CV 3.0%) | 12.5x | 2328 / 626 MB | pass (replicates 10, same_blocks_drawn True, score_rel 3.4e-12) | yes |

popmix by thread count (median wall time):

- convertf: 1: 1.83 s, 4: 1.97 s, 16: 2.10 s
- extract_f2: 1: 2.23 s, 4: 1.14 s, 16: 1.01 s
- qpadm: 1: 2.36 s, 4: 1.26 s, 16: 1.45 s
- pca: 1: 6.52 s, 4: 4.25 s, 16: 5.07 s
- kinship: 1: 3.85 s, 4: 2.30 s, 16: 2.12 s
- qpgraph: 1: 1.87 s, 4: 1.15 s, 16: 1.07 s
- qpdstat: 1: 3.30 s, 4: 1.98 s, 16: 1.98 s
- admixture: 1: 16.70 s, 4: 6.77 s, 16: 7.21 s
- rolloff: 1: 25.41 s, 4: 16.74 s, 16: 16.11 s
- dates: 1: 5.16 s, 4: 2.56 s, 16: 3.03 s
- read_kinship: 1: 1.31 s, 4: 1.05 s, 16: 1.12 s
- classic_qpdstat: 1: 2.71 s, 4: 1.86 s, 16: 1.84 s
- classic_qp3pop: 1: 3.06 s, 4: 1.78 s, 16: 1.61 s
- classic_qpadm: 1: 1.77 s, 4: 0.97 s, 16: 0.95 s
- classic_qpgraph: 1: 1.56 s, 4: 0.92 s, 16: 0.90 s
- qpadm_rotate: 1: 3.03 s, 4: 1.51 s, 16: 1.47 s
- find_graphs: 1: 6.46 s, 4: 2.89 s, 16: 2.37 s
- qpgraph_resample: 1: 2.22 s, 4: 1.60 s, 16: 1.51 s
