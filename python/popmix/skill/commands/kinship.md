# popmix kinship

Kinship for every pair of samples: KING-robust (Manichaikul et al. 2010)
for diploid genotype calls, as PLINK 2 `--make-king-table`, or READ-style
(Monroy Kuhn et al. 2018, READv2) for pseudo-haploid ancient samples.

## Synopsis

```text
popmix kinship --input PREFIX [--populations A,B,...] [--method king|read]
               [--norm median|mean|max] [--norm-value X] [--window-size BP]
               [--compat readv2-2.01] [--format tsv|json] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix. |
| `--populations` | text (comma-separated) | none (all populations except `Ignore`) | Populations whose samples are compared, pooled as one group (every pair across them). |
| `--method` | `king`, `read` | `king` | `king`: KING-robust, needs diploid calls (autosomes); `read`: READ-style, for pseudo-haploid samples, normalized within the group (every SNP in the data). |
| `--norm` | `median`, `mean`, `max` | `median` | With `--method read`: normalize each pair's P0 by this statistic over all pairs. `median` assumes most pairs are unrelated; use `max` for a small family such as a trio. |
| `--norm-value` | number | none | With `--method read`: normalize by this P0 instead (for example the median from a larger sample of the same population). |
| `--window-size` | integer | 5000000 | With `--method read`: jackknife window in base pairs (READv2 `--window_size`). |
| `--compat` | `readv2-2.01` | none (corrected) | With `--method read`: reproduce READv2 2.01 exactly (READV2-0001 to READV2-0003). |
| `--format` | `tsv`, `json` | `tsv` | Table format: tab-separated with a header, or JSON records. |
| `--output`, `-o` | text (path) | none (stdout) | Write the table here instead of stdout. |

## Inputs and outputs

A TSV table, one row per pair.

- `king`: `id1`, `id2`, `nsnp`, `hethet`, `ibs0`, `het1_hom2`, `het2_hom1`
  (counts) and `kinship`; pairs in PLINK 2's order ((2nd, 1st), (3rd, 1st),
  (3rd, 2nd), ...).
- `read`: `id1`, `id2`, `relationship`, `z_upper`, `z_lower`,
  `normalized_p0`, `normalized_p0_se`, `p0`, `p0_se`, `first_degree_type`,
  `first_degree_windows_outside`, `overlap_snps`, `effective_snps`,
  `kinship`, `kinship_ci_low`, `kinship_ci_high`; the normalization value
  is printed on stderr.

## Examples

```bash
# runs
popmix kinship -i family --populations Simulation.1
popmix kinship -i family --populations French --method read -o french_read.tsv
popmix kinship -i family --populations French --method read --norm max --window-size 10000000 \
  --compat readv2-2.01 -o french_read_compat.tsv
```

```bash
# not run: needs ancient data
popmix kinship -i aadr_v66p1 --populations France_Gurgy_MN --method read -o gurgy.tsv
```

## Reading the result

- KING `kinship`: about 0.5 for duplicates or monozygotic twins, 0.25
  first degree, 0.125 second, 0.0625 third, 0 unrelated (negative values
  are unrelated pairs, often from different populations).
- READ `normalized_p0`: about 1 unrelated, 0.875 second degree, 0.75 first
  degree, 0.5 the same individual; `kinship` is 1 minus it.
  `relationship` is READv2's call; `first_degree_type` tells
  parent-offspring from siblings for first-degree pairs; `z_upper` and
  `z_lower` are the distances, in SEs, to the neighbouring classes.

## Caveats

- KING-robust warns when pairs share fewer than 10,000 SNPs: it gives no
  standard error. [samples](samples.md) lists each sample's calls and
  heterozygosity.
- KING needs heterozygous calls: a pseudo-haploid sample gives kinship minus
  infinity, as in PLINK 2, and popmix warns; use `--method read`.
- READ normalizes within the group given, so compare samples of one
  population, with enough unrelated pairs for the median (or give
  `--norm-value`): its thresholds are fixed and take the median to be an
  unrelated pair's, so popmix warns below 10 pairs. A population with much
  background relatedness or inbreeding shifts every class.
- Coverage: pairs with few overlapping SNPs (`overlap_snps`) have wide
  intervals and READ may not call third degree; low-coverage genomes
  called pseudo-haploid compare one random allele each, which READ is
  built for, while KING-robust on low-coverage diploid calls is biased by
  genotyping error (hets missed or false) toward unrelated or related.
- Contamination adds alleles from another individual: in READ it raises
  P0, so relatives look more distant; in KING it adds false heterozygotes.
  Check contamination estimates before calling a pair unrelated.
- `--norm`, `--norm-value`, `--window-size` and `--compat` apply only with
  `--method read`.
- PLINK 2 labels its het-hom columns the wrong way round; popmix's
  `het1_hom2` is heterozygous in id1, homozygous in id2
  ([PLINK2-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#plink2-0001-het-hom-counts-labeled-the-wrong-way-round)).
- READ's corrections over READv2: disjoint windows, every window left out
  in the jackknife, and the kinship coefficient's own interval
  ([READV2-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#readv2-0001-jackknife-windows-overlap-by-one-snp)
  to [READV2-0003](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#readv2-0003-the-kinship-interval-is-the-normalized-p0s)).

## Reference

PLINK 2 `--make-king-table` (v2.0.0-a.7) and READv2 2.01: 1e-10.
