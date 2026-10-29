# The popmix Python API

Use the Python API when an analysis needs more than one command can give:
loops over many models, resampling, results fed into pandas or plots, or
admixtools2 functions that have no command (`qpadm_multi`, `qpf4ratio`,
`compare_fits`, ...). Every command's analysis is a function here, with
the same defaults.

```python
import popmix as pm
pm.init_threads(4)
f2 = pm.extract_f2("family", ["French", "Yoruba", "Simulation.1"])
```

Every name in `popmix.__all__` is listed below, and every field of every
result (tests/python/test_skill.py checks both, so this file cannot fall
behind the package). `help(pm.<name>)` gives the full docstring, with every
parameter.

Most f-statistic functions take `data` or `f2_blocks`: an `F2Blocks`, or a
genotype prefix (f2 blocks are then extracted for the populations used), as
the commands' `--input` and `--f2-dir`.

## Setup, versions and errors

| Name | What |
|------|------|
| `init_threads`, `get_threads` | Set or read the thread count (4 until set). Results are identical at any count. |
| `version`, `__version__` | The popmix version. |
| `PopmixError` | Base class of every popmix error. |
| `InputError` | A bad argument or malformed file (also a `ValueError`). |
| `PopulationNotFoundError` | A population not in the data (an `InputError`); its message suggests names. |
| `MissingFileError` | A file that is not there (also a `FileNotFoundError`). |
| `FileError` | Reading or writing failed (also an `OSError`). |
| `ComputationError` | Valid input, failed computation, such as a singular system (also a `RuntimeError`). |
| `InternalError` | A bug in popmix (also a `RuntimeError`). |
| `PopmixWarning` | A `UserWarning` for records left out and results that need care; `warnings.simplefilter("error", pm.PopmixWarning)` stops on one. |

## Reading, converting and merging data

| Name | What |
|------|------|
| `read_genotypes` | Any format, detected from the files; `populations=` reads only those. Returns a `GenotypeMatrix`. |
| `read_plink`, `read_eigenstrat`, `read_packedancestrymap`, `read_vcf` | One format each. |
| `detect_format` | The format at a prefix, such as `"PACKEDPED"` or `"VCF"`. |
| `dataset_info` | Samples, SNPs and populations, without reading the genotypes. |
| `list_populations` | (population, sample count) pairs. |
| `suggest_populations` | Names close to a mistyped one, with counts. |
| `sample_coverage` | DataFrame: each sample's SNPs called, heterozygosity and ploidy (`popmix samples`). |
| `sample_ploidy` | {sample: 1 or 2}: pseudo-haploid or diploid, by the `pseudohaploid` rule. |
| `convertf` | Convert between formats, with sample and population selection and relabelling (`popmix convertf`). |
| `convertf_par` | Run an AdmixTools convertf parameter file. |
| `merge` | Merge two datasets as AdmixTools mergeit (`popmix merge`). |

`GenotypeMatrix` has `n_snps`, `n_samples`, `samples` (`Sample`: `id`,
`population`, `ploidy`), `snps` (`Snp`: `id`, `chromosome`, `position_bp`,
`position_cm`, `allele_ref`, `allele_alt`), `populations`, `to_numpy()`
(SNPs x samples, 0/1/2 copies of the counted allele, 3 for missing) and
`allele_frequencies(populations=...)`, which returns `AlleleFrequencies`
(`populations`, `n_snps`, `n_populations`, `to_numpy()`, `counts_to_numpy()`,
`get`, `filter_populations`).

## f2 blocks

| Name | What |
|------|------|
| `extract_f2` | f2 blocks from genotypes, with admixtools2's SNP filters (`maxmiss`, `transitions`, `minmaf`, ...). Returns `F2Blocks`. |
| `extract_f2_to_dir` | The same, written straight to an f2 directory (`popmix extract-f2`). |
| `f2_from_frequencies` | f2 blocks from `AlleleFrequencies`. |
| `read_f2` | Load blocks saved by `F2Blocks.save` or `extract_f2_to_dir`. |
| `read_f2_populations` | An f2 directory's populations, from its metadata alone. |
| `read_admixtools2_f2`, `write_admixtools2_f2` | Read or write an admixtools2 f2 directory. |
| `extract_f2_subset` | These populations' blocks (also `F2Blocks.subset`). |
| `loo_list`, `boo_list` | The blocks with each block left out, or bootstrap resamples (admixtools2's). |

`F2Blocks` has `populations`, `n_populations`, `n_blocks`, `block_lengths`,
`counts_per_block`, `f2(a, b)`, `f3(a, b, c)`, `f4(a, b, c, d)` and
`fst(a, b)` (each a dict with `est`, `se`, `z`, `p` and a count),
`f2_matrix`, `f2_per_block`, `fst_per_block`, `fst_block_lengths`,
`subset`, `select_blocks`, `from_block_means`, `save(path)` and
`provenance`.

## f-statistic scans

| Name | What |
|------|------|
| `qpdstat` | f4 or D (`f4mode=False`) for many combinations, as admixtools2 (`popmix qpdstat`). DataFrame. |
| `qp3pop` | f3 for many combinations; `outgroupmode=True` for outgroup f3 (`popmix qp3pop`). DataFrame. |
| `fstat_combinations` | The combinations a scan runs, in admixtools2's order. |
| `f4blockdat_from_geno`, `f3blockdat_from_geno` | Per-block f4 (or D parts) and f3 straight from genotypes, every SNP used. |
| `outgroup_f3_matrix` | Outgroup f3 for every pair (`popmix outgroup-f3`). DataFrame. |
| `dstat_matrix` | D or f4 for every pair X, Y (`popmix dstat-matrix`). DataFrame. |
| `cluster_order` | A row order that puts similar populations together, for heatmaps. |
| `qpdstat_resample_snps`, `qp3pop_resample_snps` | The scan on each jackknife or bootstrap resample of the blocks. |
| `qpdstat_resample_inds`, `qp3pop_resample_inds` | The scan with each individual left out in turn. |

## qpWave and qpAdm

| Name | What |
|------|------|
| `qpwave` | Rank tests of `left` against `right` (`popmix qpwave`). Returns `QpWaveResult`. |
| `qpadm` | `target` as a mixture of `left`, tested against `right` (`popmix qpadm`). Returns `QpAdmResult`. |
| `qpadm_p` | One model at one rank, without SEs or popdrop: fast, for screening. |
| `qpadm_rotate`, `rotate_models` | Every model of a rotation, or just the models (`popmix qpadm-rotate`). |
| `qpadm_multi`, `qpadm_models` | Several models at once; every left/right split of a population list. |
| `qpwave_pairs`, `test_cladality` | Cladality of every pair of `left`, or of all of them, with respect to `right`. |
| `qpf4ratio` | The f4-ratio share alpha of `source_a` in `test`. Returns `F4RatioResult`. |
| `qpadm_resample_snps`, `qpwave_resample_snps` | The fit on each resample of the blocks. |
| `qpadm_resample_inds`, `qpwave_resample_inds` | The fit with each individual left out in turn. |
| `calibrate_qpadm`, `calibrate_f4` | How often a model is rejected, or an f4 non-zero, over simulated replicates (`popmix calibrate-qpadm`, `calibrate-f4`). |

### `QpWaveResult`

- `left`, `right`: the populations.
- `rank_tests`: one `RankTest` per rank (`rank`, `chisq`, `dof`, `p`).
- `estimated_rank`: the lowest rank not rejected at `alpha`; `full_rank()`
  is the test of the highest rank (the last of `rank_tests`).
- `f4_matrix`, `f4_covariance`: the f4 matrix and its jackknife covariance.

### `QpAdmResult`

- `target`, `sources`, `right`: the populations.
- `weights`, `weight_se`, `weight_z`, `weights_dict()`: admixture weights,
  their SEs, weight over SE, and {source: weight}.
- `chisq`, `dof`, `p`: the model's fit; `feasible`: every weight in [0, 1].
- `constrained`: whether the weights were constrained to [0, 1].
- `popdrop`: one `PopDropRow` per subset of the sources (`pattern`, 1 for a
  source left out; `wt`, `dof`, `chisq`, `p`, `f4rank`, `weights`,
  `feasible`; and admixtools2's nested path `best`, `dofdiff`, `chisqdiff`,
  `p_nested`).
- `rankdrop`: `RankTest`s of the f4 rank.
- `n_blocks`, `n_snps`: the data used.
- `assess(alpha=0.05)`: a `QpAdmAssessment` with `status` (`"not tested"`,
  `"rejected"`, `"infeasible"` or `"not rejected"`), `alpha` and `notes`
  on what weakens a model that fits (weights near 0, simpler models that
  also fit). Read it before calling a model good: not rejected is not
  proven.

### `F4RatioResult`

`test`, `source_a`, `source_b`, `outgroup`, `reference`, `alpha` (the
share), `se`, `z`, `n_blocks`, `f4_numerator`, `f4_denominator`.

## Admixture graphs

| Name | What |
|------|------|
| `qpgraph` | Fit a graph, given as (from, to) edges (`popmix qpgraph`). Returns `QpGraphResult`. |
| `build_tree` | A tree (graph without admixture) from its nodes; returns `AdmixGraph`. |
| `find_graphs` | Search for the best-fitting graphs (`popmix find-graphs`). Returns (graphs, number fitted): `ScoredGraph`s, best first. |
| `all_graphs` | Every graph with these leaves and `nadmix` admixture nodes. |
| `graph_key` | admixtools2's `graph_hash` key, to tell graphs apart. |
| `qpgraph_out_of_sample` | A graph's score fitted on alternate jackknife blocks. |
| `qpgraph_resample_snps`, `qpgraph_resample_multi`, `compare_fits` | Fits on resamples, several graphs on the same bootstrap resamples, and whether one graph fits better than another. |
| `qpgraph_resample_inds` | The fit with each individual left out in turn. |
| `write_dot`, `parse_dot` | A graph to and from DOT. |

`AdmixGraph` has `leaves`, `node_count`, `edge_count`, `num_drift_edges`,
`get_edge_weights` and `expected_f2`.

`ScoredGraph`: `edges` (the root-to-first-population path first) and
`score`.

### `QpGraphResult`

- `pops`: the leaves; `score`: the fit (lower is better);
  `penalized_score`, `n_f3`, `n_params`, `dof`; `score_test`: the
  out-of-sample score, when computed.
- `edges`: `QpGraphEdge`s (`from_node`, `to_node`, `kind` drift or admix,
  `weight`, and `low`, `high` bounds where known).
- `f3`: `QpGraphF3`s, each f3 as estimated and as fitted (`pop1`, `pop2`,
  `pop3`, `est`, `se`, `fit`, `diff`, `z`, `p`). The worst `z` shows where
  the graph fails.
- `opt`: one `QpGraphStart` per random start (`start`, `weights`, `score`,
  `evaluations`, `convergence`, `converged`).
- `ppinv`: the inverse f3 covariance used in the score.

## Dating admixture

| Name | What |
|------|------|
| `dates` | DATES: the date of `admixed` from two `sources` (`popmix dates`). Returns `DatesResult`. |
| `rolloff` | ROLLOFF (`popmix rolloff`). Returns `RolloffResult`. |

### `DatesResult`

- `admixed`, `sources`.
- `generations`, `generations_se`: the jackknife date and SE DATES reports.
- `years_bp`, `years_bp_se`: `sample_age` + generations x
  `generation_time`, with every SE combined (`generation_time`,
  `generation_time_se`, `sample_age`, `sample_age_se` as given).
- `n_snps`, `n_chromosomes`; `source_alleles`: each source's mean alleles
  per SNP (below 10 is a warning).
- `jackknife_by_chromosome`: (chromosome, estimate, weight).
- `ld_bins`: `LdBin`s (`distance_cm`, `n_pairs`, `weighted_ld`, `se`, `z`,
  `correlation`, `regression_coef`); `to_ld_dataframe()` gives them as a
  DataFrame.
- `fit`: the `DecayFit` (`amplitude`, `decay_rate`, `decay_rate_se`,
  `intercept`, `r_squared`, `rss`).
- `pulse_check`: a `PulseCheck` (`generations`, `generations_se`,
  `amplitudes`, `amplitudes_se`, `difference_se`, `rss`,
  `multiple_pulses`), or None.

### `RolloffResult`

`bins` (DataFrame: d_cm, zscore, regression, corr, count), `decay` (the
date in generations), `affine`, `amplitude`, `decay_jackknife` (mean, SE
over chromosomes), `chromosomes`, `n_snps_used`, `jackknife_bins` and
`pulse_check`.

## Structure, ancestry and relatedness

| Name | What |
|------|------|
| `smartpca` | PCA as smartpca, with projection and Tracy-Widom tests (`popmix pca`). Returns `PcaResult`. |
| `admixture` | ADMIXTURE: unsupervised, `supervised=`, `fixed_p=` projection, `cv=` (`popmix admixture`). Returns `AdmixtureResult`. |
| `admixture_loglik` | ADMIXTURE's log-likelihood of a given Q and P. |
| `king` | KING-robust kinship of every pair, as PLINK 2 (`popmix kinship`). DataFrame. |
| `read_kinship` | READ kinship for pseudo-haploid samples, as READv2 (`popmix kinship --method read`). DataFrame. |
| `sfs` | Site-frequency spectrum projected as dadi (`popmix sfs`). Returns `Sfs`. |
| `joint_sfs` | Joint spectrum as admixtools2 `joint_sfs()` (`popmix sfs --joint`). DataFrame. |

### `PcaResult`

`eigenvalues`, `coords` (id, pop, role and PC1, PC2, ...),
`tracy_widom` (eigenvalue, twstat, p, neff), `outliers` (removed samples
with the iteration, pc and score), `no_data`, `n_pca_samples`,
`n_snps_used`, `n_snps_killed`, `shrinkage` (each component's measured
shrinkage with `holdout`) and `missing_fraction`.

### `AdmixtureResult`

`Q` (DataFrame: id, pop, Q1..QK), `q_matrix()` (Q as an array), `k`, `P`
(SNPs x K), `snps`, `loglik`, `iterations`, `converged` and `cv_error`.

### `Sfs`

`pops`, `data` (counts by allele count in each population), `mask`,
`folded`, `not_integer` (SNPs left out for half an allele),
`sample_sizes`, and the statistics `S`, `pi`, `watterson_theta` and
`tajima_d`.
