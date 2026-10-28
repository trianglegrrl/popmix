# popmix rolloff

ROLLOFF (ADMIXTOOLS 7.0.2): date admixture from the decay with genetic
distance of the correlation between the admixed population's LD and the
two references' allele-frequency differences (Moorjani et al. 2011).
[dates](dates.md) is its successor and is preferred for new work; ROLLOFF
needs several admixed samples.

## Synopsis

```text
popmix rolloff --input PREFIX --admixed T --sources A,B [--binsize M] [--maxdis M]
               [--runmode 0|1] [--chithresh X] [--chromosomes 1,2,...] [--lval CM] [--hval CM]
               [-o PREFIX]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix with a genetic map. |
| `--admixed` | text | required | The admixed population. |
| `--sources` | text (comma-separated) | required | The two source (reference) populations. |
| `--binsize` | number | 0.0005 | Distance bin in Morgans (0.05 cM). |
| `--maxdis` | number | 0.05 | Largest SNP-pair distance in Morgans (5 cM). |
| `--runmode` | `0`, `1` | `1` | Weights: `1` pA - pB (recommended, the one tested against ROLLOFF); `0` divided by the square root of the mean frequency. |
| `--chithresh` | number | 6 | Leave out SNP pairs in LD in a reference (chi-square at least this). |
| `--chromosomes` | text (comma-separated) | none (all autosomes) | Use only these chromosomes. |
| `--lval` | number | 0.5 | Fit the decay from this distance in cM. |
| `--hval` | number | 100 | Fit the decay up to this distance in cM. |
| `--output-prefix`, `-o` | text (path) | none | Also write `PREFIX.bins.tsv`: one row per distance bin (d_cm, zscore, regression, corr, count). |

## Inputs and outputs

Prints one line:

```text
decay 500.000 generations; jackknife 807.373 +/- 356.316 (3 chromosomes, 2018 SNPs)
```

and with `-o` the bins table. Bins are labeled by their upper edge, as
ROLLOFF does.

## Examples

```bash
# runs
popmix rolloff -i family --admixed Simulation.1 --sources French,Yoruba -o rolloff_sim
popmix rolloff -i family --admixed Simulation.1 --sources French,Yoruba --chromosomes 20,21 --maxdis 0.02
```

```bash
# not run: needs the AADR
popmix rolloff -i aadr_v66p1 --admixed ASW.DG --sources CEU.DG,YRI.DG -o asw
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

`decay` is the full-data fit C + A exp(-m d) (m in generations) between 0.5
and 100 cM; the jackknife over chromosomes gives the estimate and SE to
report. The decay is bounded to (0, 500]: a value at 500, as in the toy
example above (4 admixed samples), means the curve does not determine a
date. Check the `corr` column of the bins for a smooth decay.

## Caveats

- ROLLOFF needs many admixed samples; with one or a few, use
  [dates](dates.md).
- The fit is kept within (0, 500] generations, where ROLLOFF's nls is
  unbounded; on noisy curves the two then differ, and where the curve
  determines a date they agree
  ([ROLLOFF matched behaviour](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#rolloff-admixtools-702)).
- The jackknife weights each chromosome present by its SNPs; ROLLOFF's
  expfit.sh assumes chromosomes 1-22.
- `zdipcorrmode`, `mincount` and `minparentcount` are fixed on the CLI;
  the Python API (`pm.rolloff`) takes them.

## Reference

ROLLOFF from ADMIXTOOLS 7.0.2 with rexpfit.r's fit and expfit.sh's weighted
jackknife: curves to 1e-9.
