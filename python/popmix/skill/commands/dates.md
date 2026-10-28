# popmix dates

DATES: estimate when a population was admixed from the decay, with genetic
distance, of the weighted ancestry covariance (LD) of its samples, using two
populations as proxies of the sources (Narasimhan et al. 2019, Chintalapati
et al. 2022). Works with a single admixed sample and with pseudo-haploid
ancient samples.

## Synopsis

```text
popmix dates --input PREFIX --admixed T --sources A,B [--bin-width CM]
             [--max-distance CM] [--min-distance CM] [--generation-time Y] [--generation-time-se Y]
             [--sample-age Y] [--sample-age-se Y] [--no-jackknife]
             [--qbin N] [--compat dates-v753]
             [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix, with a genetic map in the `.snp` (or `.bim`). All genotypes are read into memory. |
| `--admixed` | text | required | The admixed population (one sample or more): DATES's test population. |
| `--sources` | text (comma-separated) | required | Proxies for its two sources. |
| `--bin-width` | number | 0.1 | Distance bin width in cM (DATES `binsize` x 100). |
| `--max-distance` | number | 100 | Largest SNP-pair distance in cM (DATES `maxdis` x 100). Smaller values (for example 50) are much faster. |
| `--min-distance` | number | 0.45 | Fit only bins from this distance in cM (DATES `lovalfit`); shorter distances carry background LD. |
| `--generation-time` | number | 29 | Years per generation, for the date in years. |
| `--generation-time-se` | number | 0 | Its standard error in years; enters the date's SE. |
| `--sample-age` | number | 0 | Age of the admixed samples in years BP, so the date is before present (0: before the samples lived). |
| `--sample-age-se` | number | 0 | Its standard error in years. |
| `--jackknife` / `--no-jackknife` | flag | `--jackknife` | Standard errors by a weighted jackknife over chromosomes; `--no-jackknife` is faster and gives no SE. |
| `--qbin` | integer (at least 0) | 0 | DATES `qbin`: 0 visits every SNP pair (exact); N > 0 bins SNPs on an FFT mesh of N cells per bin, much faster on large data and slightly approximate (10 is typical). |
| `--compat` | `dates-v753` | none (corrected) | Reproduce DATES v753's quirks (DATES-0001 to DATES-0004) for output identical to DATES. |
| `--format` | `text`, `json`, `tsv` | `text` | `text`: the date and the first 10 bins; `json`: admixed, sources, date, SE, years, generation_time(_se), sample_age(_se), source_alleles, fit (amplitude, decay_rate, r_squared), SNPs and the first 20 bins; `tsv`: every LD bin (distance_cm, n_pairs, correlation, weighted_ld). |

## Inputs and outputs

Progress goes to stderr; the result to stdout:

```text
DATES Results: Simulation.1 = mix of French and Yoruba
============================================================
Amplitude (fit): 0.031292
Estimated generations: 65.0 ± 25.8
Date: 1884 ± 749 years before the admixed samples lived (29 years/gen, taken as
      exact; give --sample-age for years BP)
Fit R²: 0.0695 (how well one exponential fits; it does not show a single pulse)
SNPs used: 20881
LD bins: 1000
```

followed by the first LD bins.

## Examples

```bash
# runs
popmix dates -i family --admixed Simulation.1 --sources French,Yoruba
popmix dates -i family --admixed Simulation --sources French,Yoruba --qbin 10 --format json
popmix dates -i family --admixed Simulation.1 --sources French,Yoruba \
  --max-distance 50 --format tsv > ld_decay.tsv
popmix dates -i family --admixed Simulation.1 --sources French,Yoruba --compat dates-v753
```

```bash
# not run: needs the AADR
popmix --threads 8 dates -i aadr_v66p1 --admixed French.DG \
  --sources Sardinian.DG,Russia_Samara_EBA_Yamnaya --qbin 10 --max-distance 50
```

## Reading the result

The two-pulse check fits A₁·e^(−n₁·d) + A₂·e^(−n₂·d) + K to the same
bins, with jackknife SEs over chromosomes (JSON `pulse_check`; Python
`pulse_check`, a `PulseCheck`). Two pulses are reported only when both
amplitudes and n₂ − n₁ are more than 2 SE from zero: then the single date
averages admixture at different times. Otherwise one pulse is not
rejected, which does not show a single pulse: admixture over a short span
fits one exponential too, and the test has little power with few
chromosomes. It needs the jackknife.

- **Estimated generations ± SE**: a weighted block jackknife over
  chromosomes (what DATES reports); the full-data fit's rate is
  `fit.decay_rate` in JSON.
- **Date** in years = `--sample-age` + generations x `--generation-time`
  (29 by default). Without a sample age it is years before the admixed
  samples lived, not BP. The generation time is an assumption (25 to 30
  years are common); give its uncertainty with `--generation-time-se`, and
  the date's SE includes it and the sample age's.
- **Fit R²** says how well one exponential fits the decay, not that there
  was a single pulse: continuous or repeated admixture also gives a decay
  that one exponential fits, with a date between the events.
- A date is credible when the decay curve (the `tsv` bins) shows a clear
  exponential above noise, the SE is well below the estimate (Z = estimate
  / SE above about 2-3) and the date is stable across `--min-distance` and
  source choices. Here the 4 Simulation.1 samples on 3 chromosomes give 65
  ± 26 against a true 50: the toy data is small.
- The amplitude scales with the admixture proportion and with how well the
  proxies represent the sources. Very old dates (hundreds of generations)
  are hard to tell from background LD.

## Caveats

- The jackknife is over chromosomes, so use data on many chromosomes;
  popmix warns below 10.
- Corrections over DATES v753 (default) and `--compat dates-v753`:
  mesh FFT lags that wrap around and the last lag left unscaled
  ([DATES-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0001-mesh-fft-lags-wrap-around),
  [DATES-0002](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0002-last-fft-lag-left-unscaled), only with `--qbin`),
  fitting printed values
  ([DATES-0003](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0003-fitting-printed-values)),
  and pseudo-haploid calls counted as two alleles in the sources
  ([DATES-0004](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#dates-0004-every-sample-counted-as-diploid)).
- DATES's `mincount` is ignored by DATES; popmix matches (a SNP is used
  once each source has one called sample). So a source of few samples, or
  of pseudo-haploid samples with much missing data, gives noisy weights:
  popmix warns when a source averages fewer than 10 observed alleles per
  SNP used (`source_alleles` in JSON).

## Reference

DATES v753 (Narasimhan et al. 2019; Chintalapati, Patterson and Moorjani
2022), matched to its printed precision with `--compat dates-v753`.
