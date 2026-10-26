# popmix dates

DATES: estimate when a population was admixed from the decay, with genetic
distance, of the weighted ancestry covariance (LD) of its samples, using two
populations as proxies of the sources (Narasimhan et al. 2019, Chintalapati
et al. 2022). Works with a single admixed sample and with pseudo-haploid
ancient samples.

## Synopsis

```text
popmix dates --input PREFIX --test-pop T --source1 A --source2 B [--bin-width CM]
             [--max-distance CM] [--min-distance CM] [--qbin N] [--compat dates-v753]
             [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix, with a genetic map in the `.snp` (or `.bim`). All genotypes are read into memory. |
| `--test-pop` | text | required | The admixed population (one sample or more). |
| `--source1` | text | required | A proxy for the first source. |
| `--source2` | text | required | A proxy for the second source. |
| `--bin-width` | number | 0.1 | Distance bin width in cM (DATES `binsize` x 100). |
| `--max-distance` | number | 100 | Largest SNP-pair distance in cM (DATES `maxdis` x 100). Smaller values (for example 50) are much faster. |
| `--min-distance` | number | 0.45 | Fit only bins from this distance in cM (DATES `lovalfit`); shorter distances carry background LD. |
| `--qbin` | integer (at least 0) | 0 | DATES `qbin`: 0 visits every SNP pair (exact); N > 0 bins SNPs on an FFT mesh of N cells per bin, much faster on large data and slightly approximate (10 is typical). |
| `--compat` | `dates-v753` | none (corrected) | Reproduce DATES v753's quirks (DATES-0001 to DATES-0004) for output identical to DATES. |
| `--format` | `text`, `json`, `tsv` | `text` | `text`: the date and the first 10 bins; `json`: date, SE, years, fit (amplitude, decay_rate, r_squared), SNPs and the first 20 bins; `tsv`: every LD bin (distance_cm, n_pairs, correlation, weighted_ld). |

## Inputs and outputs

Progress goes to stderr; the result to stdout:

```text
DATES Results: Simulation.1 = mix of French and Yoruba
============================================================
Amplitude (fit): 0.031292
Estimated generations: 65.0 ± 25.8
Estimated years BP: 1884 ± 749 (29 years/gen)
Fit R²: 0.0695
SNPs used: 20881
LD bins: 1000
```

followed by the first LD bins.

## Examples

```bash
# runs
popmix dates -i family --test-pop Simulation.1 --source1 French --source2 Yoruba
popmix dates -i family --test-pop Simulation --source1 French --source2 Yoruba --qbin 10 --format json
popmix dates -i family --test-pop Simulation.1 --source1 French --source2 Yoruba \
  --max-distance 50 --format tsv > ld_decay.tsv
popmix dates -i family --test-pop Simulation.1 --source1 French --source2 Yoruba --compat dates-v753
```

```bash
# not run: needs the AADR
popmix --threads 8 dates -i /references/aadr/v62.0_1240k_public --test-pop French.DG \
  --source1 Sardinian.DG --source2 Russia_Samara_EBA_Yamnaya --qbin 10 --max-distance 50
```

## Reading the result

- **Estimated generations ± SE**: a weighted block jackknife over
  chromosomes (what DATES reports); the full-data fit's rate is
  `fit.decay_rate` in JSON. Years = generations x 29 (fixed); for ancient
  samples add their age.
- A date is credible when the decay curve (the `tsv` bins) shows a clear
  exponential above noise, the SE is well below the estimate (Z = estimate
  / SE above about 2-3) and the date is stable across `--min-distance` and
  source choices. Here the 4 Simulation.1 samples on 3 chromosomes give 65
  ± 26 against a true 50: the toy data is small.
- The amplitude scales with the admixture proportion and with how well the
  proxies represent the sources. Very old dates (hundreds of generations)
  are hard to tell from background LD.

## Caveats

- The jackknife is over chromosomes, so use data on many chromosomes.
- Corrections over DATES v753 (default) and `--compat dates-v753`:
  mesh FFT lags that wrap around and the last lag left unscaled
  ([DATES-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0001-mesh-fft-lags-wrap-around),
  [DATES-0002](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0002-last-fft-lag-left-unscaled), only with `--qbin`),
  fitting printed values
  ([DATES-0003](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0003-fitting-printed-values)),
  and pseudo-haploid calls counted as two alleles in the sources
  ([DATES-0004](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0004-every-sample-counted-as-diploid)).
- DATES's `mincount` is ignored by DATES; popmix matches (a SNP is used
  once each source has one called sample).
- The generation time (29 years) and the jackknife cannot be changed on the
  CLI; the Python API (`pm.dates(..., generation_time=...)`) takes them.

## Reference

DATES v753 (Narasimhan et al. 2019; Chintalapati, Patterson and Moorjani
2022), matched to its printed precision with `--compat dates-v753`.
