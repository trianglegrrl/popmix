Measured 2026-10-27 on 13th Gen Intel(R) Core(TM) i7-1360P (16 threads, 62.4 GB); popmix 1.0.0 (c44ea2b), admixtools2 2.0.10; median of 5 runs after a warm-up run. Reference times are from the run of 2026-10-25 (the reference tools are unchanged; agreement is checked against that run's outputs); popmix was run again.

| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix | Agreement | Same output at every thread count |
|---|---|---|---|---|---|---|---|
| convertf: packed to EIGENSTRAT text | AdmixTools 8.0.2 convertf | 77.96 s (CV 8.4%) | 1.61 s (4; CV 15.1%) | 48.5x | 1018 / 278 MB | pass (snps 1233013, differing_snps 0) | yes |
| extract_f2: f2 blocks, 16 populations | admixtools2 2.0.10 f2_from_geno | 36.90 s (CV 1.0%) | 0.93 s (16; CV 2.6%) | 39.5x | 7885 / 341 MB | pass (blocks 710, max_abs 2.4e-15, differing 0) | yes |
| qpadm: qpAdm from genotypes (2 sources, 7 right) | admixtools2 2.0.10 qpadm | 21.47 s (CV 2.9%) | 0.95 s (16; CV 16.4%) | 22.5x | 3731 / 335 MB | pass (weights_rel 1.5e-15, se_rel 2.0e-12, p_rel 4.0e-12) | yes |
| pca: PCA, 8 present-day populations, 10 PCs | EIGENSOFT 7.2.1 smartpca | 83.24 s (CV 8.5%) | 4.29 s (16; CV 9.0%) | 19.4x | 755 / 641 MB | pass (eigenvalues 10, max_rel 4.3e-12) | yes |
| kinship: KING-robust kinship, every pair of 8 populations' samples | PLINK 2 --make-king-table (16 threads) | 1.23 s (CV 17.4%) | 1.08 s (16; CV 35.3%) | 1.1x | 130 / 299 MB | pass (pairs 24090, max_abs 5.0e-13) | yes |
| qpgraph: qpGraph, 6 populations, 1 admixture, 10 starts | admixtools2 2.0.10 qpgraph | 16.78 s (CV 3.5%) | 0.85 s (16; CV 11.2%) | 19.7x | 2321 / 326 MB | pass (score 3.0e+02, score_rel 1.4e-12, seeds 3, popmix_spread 4.0e-10, reference_spread 0.0e+00) | yes |
| qpdstat: 20 f4-statistics from genotypes | admixtools2 2.0.10 qpdstat | 16.75 s (CV 2.6%) | 1.58 s (16; CV 1.1%) | 10.6x | 670 / 413 MB | pass (statistics 20, est_rel 3.4e-11, se_rel 7.6e-12) | yes |
| admixture: ADMIXTURE K=3, 3 present-day populations, LD-pruned, seed 1 | ADMIXTURE 1.3.0 (1 thread) | 13.87 s (CV 2.7%) | 3.59 s (4; CV 14.0%) | 3.9x | 41 / 136 MB | pass (loglik -4.4e+06, reference_loglik -4.4e+06, loglik_diff 4.2e-03, seeds 3, popmix_spread 0.0e+00, reference_spread 7.9e-02) | yes |
| rolloff: ROLLOFF, French = Sardinian + Russia_Samara_EBA_Yamnaya, to 5 cM, chromosomes 21-22, jackknife over them | ADMIXTOOLS 7.0.2 rolloff (full-precision build) | 4.3 min (CV 0.9%) | 15.66 s (16; CV 6.1%) | 16.3x | 915 / 600 MB | pass (bins 50, zscore_rel 3.9e-12, regression_rel 4.1e-12, corr_rel 3.7e-12, bins_equal True) | yes |
| dates: DATES, French = Sardinian + Russia_Samara_EBA_Yamnaya, FFT mesh (qbin 10), to 50 cM, jackknife | DATES v753 (numbins fix, as the goldens) | 30.73 s (CV 1.3%) | 2.26 s (4; CV 7.7%) | 13.6x | 758 / 608 MB | pass (bins_compared 20, pairs_equal True, max_cov_abs 4.6e-10, max_corr_abs 5.0e-07, fit_abs 3.2e-04, generations_abs 2.0e-03, se_abs 3.4e-04) | yes |
| read_kinship: READ kinship, every pair of Spain_EN samples | READv2 (73739be) | 28.26 s (CV 1.2%) | 0.80 s (4; CV 3.4%) | 35.5x | 950 / 193 MB | pass (pairs 253, classes_equal True, max_rel 4.8e-12) | yes |
| classic_qpdstat: 20 f4-statistics from genotypes | ADMIXTOOLS 7.0.2 qpDstat | 42.89 s (CV 1.5%) | 1.54 s (16; CV 15.5%) | 27.8x | 2169 / 404 MB | pass (statistics 20, est_abs 4.8e-07, z_abs 4.8e-04, snps_equal True) | yes |
| classic_qp3pop: 16 outgroup f3-statistics from genotypes | ADMIXTOOLS 7.0.2 qp3Pop | 2.4 min (CV 0.6%) | 1.51 s (16; CV 5.7%) | 96.8x | 755 / 399 MB | pass (statistics 16, est_abs 4.0e-10, se_abs 5.0e-10, z_abs 5.0e-04) | yes |
| classic_qpadm: qpAdm from genotypes (2 sources, 7 right) | ADMIXTOOLS 7.0.2 qpAdm | 36.61 s (CV 1.3%) | 0.82 s (16; CV 12.0%) | 44.5x | 1173 / 321 MB | FAIL (weights_abs 9.6e-04, se_abs 7.8e-04, chisq_abs 1.8e-01, p_rel 7.0e-02) | yes |
| classic_qpgraph: qpGraph, 6 populations, 1 admixture | ADMIXTOOLS 7.0.2 qpGraph | 90.98 s (CV 1.2%) | 0.75 s (4; CV 16.7%) | 121.5x | 760 / 300 MB | pass (score 3.1e+02, score_abs 3.9e-04) | yes |
| qpadm_rotate: qpAdm rotation, 5 populations between left and right, 4 fixed right | admixtools2 2.0.10 qpadm_rotate | 21.35 s (CV 4.8%) | 1.33 s (16; CV 16.5%) | 16.1x | 3765 / 332 MB | pass (models 30, rank_dof_feasible_equal True, chisq_rel 7.0e-12, p_rel 9.8e-10) | yes |
| find_graphs: find_graphs, 6 populations, 1 admixture, Chimp as outgroup: popmix exhaustive, admixtools2 evolutionary (20 generations, seed 1) | admixtools2 2.0.10 find_graphs | 45.24 s (CV 1.3%) | 2.15 s (16; CV 14.3%) | 21.1x | 2328 / 325 MB | pass (graphs_fitted 1965, score 3.0e+02, reference_score 3.0e+02, score_rel -1.4e-12, seeds 3, reference_spread 8.0e-12) | yes |
| qpgraph_resample: qpGraph bootstrap, 10 resamples of the blocks, 6 populations, 1 admixture | admixtools2 2.0.10 qpgraph_resample_snps | 18.87 s (CV 4.8%) | 1.25 s (16; CV 13.8%) | 15.1x | 2328 / 325 MB | pass (replicates 10, same_blocks_drawn True, score_rel 3.4e-12) | yes |

popmix by thread count (median wall time):

- convertf: 1: 1.73 s, 4: 1.61 s, 16: 1.68 s
- extract_f2: 1: 2.13 s, 4: 0.96 s, 16: 0.93 s
- qpadm: 1: 2.09 s, 4: 1.17 s, 16: 0.95 s
- pca: 1: 6.46 s, 4: 4.35 s, 16: 4.29 s
- kinship: 1: 2.37 s, 4: 1.21 s, 16: 1.08 s
- qpgraph: 1: 1.68 s, 4: 0.98 s, 16: 0.85 s
- qpdstat: 1: 2.67 s, 4: 1.77 s, 16: 1.58 s
- admixture: 1: 9.29 s, 4: 3.59 s, 16: 3.99 s
- rolloff: 1: 25.71 s, 4: 15.83 s, 16: 15.66 s
- dates: 1: 4.53 s, 4: 2.26 s, 16: 3.34 s
- read_kinship: 1: 0.90 s, 4: 0.80 s, 16: 0.80 s
- classic_qpdstat: 1: 2.36 s, 4: 1.72 s, 16: 1.54 s
- classic_qp3pop: 1: 2.49 s, 4: 1.67 s, 16: 1.51 s
- classic_qpadm: 1: 1.65 s, 4: 1.02 s, 16: 0.82 s
- classic_qpgraph: 1: 1.32 s, 4: 0.75 s, 16: 0.80 s
- qpadm_rotate: 1: 2.62 s, 4: 1.50 s, 16: 1.33 s
- find_graphs: 1: 5.87 s, 4: 2.58 s, 16: 2.15 s
- qpgraph_resample: 1: 2.15 s, 4: 1.32 s, 16: 1.25 s
