# popmix qpwave

qpWave: how many independent streams of ancestry, as seen from the right
populations, are needed to relate the left populations? Rank r of the f4
matrix f4(L0, Li; R0, Rj) means at least r + 1 streams. The rank is a lower
bound: other right populations, or more data, can show more streams, so not
rejecting rank 0 does not show that the left populations form a clade.

## Synopsis

```text
popmix qpwave --left L1,L2,... --right R1,R2,... (--input PREFIX | --f2-dir DIR) [options] [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--left` | text (comma-separated) | required | Left populations, at least two; the first is the basis L0. |
| `--right` | text (comma-separated) | required | Right populations, at least two; the first is the basis R0. |
| `--alpha` | number (0 to 1) | 0.05 | Significance level of the rank tests: the estimated rank is the lowest rank with p above it. |
| `--fudge` | number | 1e-4 | Regularisation of the f4 covariance: this fraction of its trace is added to its diagonal (admixtools2 `fudge`). |
| `--input`, `-i` | text (path) | none | Genotype prefix. f2 blocks are extracted for exactly this command's populations, on the SNPs with data in all of them (with `--maxmiss 0`). Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Read f2 blocks saved by [extract-f2](extract-f2.md) instead of genotypes; every population must be in them. The SNP set is the one extract-f2 chose. |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--maxmiss` | number (0 to 1) | 0 | With `--input`: keep SNPs missing in at most this fraction of the populations. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--transversions-only` | flag | off | With `--input`: use only transversions (drop A<->G and C<->T SNPs), for ancient DNA without UDG treatment, as admixtools2 `transitions = FALSE`. A SNP without two known bases is dropped. |
| `--transitions-only` | flag | off | With `--input`: use only transitions, as admixtools2 `transversions = FALSE`. |
| `--minmaf` | number (0-0.5) | 0 | With `--input`: drop SNPs whose minor allele frequency, the populations' alleles pooled, is below this (admixtools2 `minmaf`). |
| `--maxmaf` | number (0-0.5) | 0.5 | With `--input`: drop SNPs whose pooled minor allele frequency is above this (admixtools2 `maxmaf`). |
| `--keep-snps` | text (path) | none | With `--input`: use only the SNP IDs in this file (one per line, or a `.snp`/`.bim` file). As admixtools2 `keepsnps`, it overrides `--maxmiss`, the MAF bounds, `--transversions-only` and the autosome restriction (chromosome X SNPs listed are used); polymorphism and `--exclude-snps` still apply. |
| `--exclude-snps` | text (path) | none | With `--input`: never use the SNP IDs in this file (popmix's own; admixtools2 has none). |
| `--keep-inds` | text (path) | none | With `--input`: use only the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--exclude-inds` | text (path) | none | With `--input`: never use the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--relabel` | text (path) | none | With `--input`: give samples new populations: `ID POPULATION` per line, as admixtools2 `inds`/`pops`. Population names given to the command are those after relabelling. |
| `--relabel-pops` | text (path) | none | With `--input`: rename populations: `OLD NEW` per line, applied after `--relabel` (merge two populations by giving them one name). |
| `--format` | `text`, `json`, `tsv` | `text` | `json`: left, right, alpha, estimated_rank, rank_tests (rank, chisq, dof, p). `tsv`: one row per rank test. |

## Inputs and outputs

Genotypes or an f2 directory; the result on stdout:

```text
qpWave: Simulation.1, French_1, Yoruba_1
Outgroups: French_2, Yoruba_2, French_3, Yoruba_3

Rank tests (* rejected at alpha 0.05):
  Rank   chisq    dof          p
     0  1240.0343     6  0.000000 *
     1    1.2836     2  0.526345

Estimated rank: 1 (the lowest rank not rejected)
=> At least 2 streams of ancestry, as seen from these right populations: a
   lower bound, as other right populations or more data can show more streams
   (the rank can be at most 2 here).
```

## Examples

```bash
# runs
popmix qpwave -i family6 --left Simulation.1,French_1,Yoruba_1 --right French_2,Yoruba_2,French_3,Yoruba_3
popmix qpwave -i family6 --left French_1,French_2 --right Yoruba_1,Yoruba_2,Yoruba_3 --format json
```

```bash
# not run: needs the AADR
popmix qpwave -i aadr_v66p1 --left French.DG,Sardinian.DG \
  --right Mbuti.DG,Han.DG,Karitiana.DG,Papuan.DG
```

## Reading the result

Each row tests "rank at most r" against full rank; `*` marks p at or below
`--alpha` (rejected). The estimated rank is the lowest rank not rejected:
at least that many plus one streams. It is a lower bound: the f4 matrix's
rank is at most min(left, right) - 1, right populations related alike to
the left ones hide streams, and a test without power reads as "not
rejected". So rank 0 not rejected means these right populations do not
distinguish the left ones, not that they form a clade. A qpAdm model of
a target with n sources is the qpWave test of rank n - 1 on the target plus
the sources (the last row above equals the qpAdm chi-squared).

## Caveats

- chi-squared treats the jackknife covariance of the f4 statistics as
  known. With g blocks and p = (left - 1)(right - 1) statistics it is
  inflated about (g - 1)/(g - p - 2)-fold (Hotelling's T²), so p is too
  small; popmix warns past 10% and when g <= p + 2, and below 20 blocks.
- The ranks tested go from 0 to (number of left populations - 2). If every
  test is rejected the estimated rank is one above the highest rank tested:
  full rank, as many streams as left populations.
- The same SNP-set note as [qpadm](qpadm.md) applies to `--input`.

## Reference

admixtools2 `qpwave()` (2.0.10), `drop_ranks`; Reich et al. 2012.
