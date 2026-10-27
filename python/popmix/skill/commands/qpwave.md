# popmix qpwave

qpWave: how many independent streams of ancestry, as seen from the right
populations, are needed to relate the left populations? Rank r of the f4
matrix f4(L0, Li; R0, Rj) means r + 1 streams; rank 0 means the left
populations form a clade with respect to the right ones.

## Synopsis

```text
popmix qpwave --left L1,L2,... --right R1,R2,... (--input PREFIX | --f2-dir DIR) [options] [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--left` | text (comma-separated) | required | Left populations, at least two; the first is the basis L0. |
| `--right` | text (comma-separated) | required | Right populations, at least two; the first is the basis R0. |
| `--fudge` | number | 1e-4 | Regularisation of the f4 covariance: this fraction of its trace is added to its diagonal (admixtools2 `fudge`). |
| `--input`, `-i` | text (path) | none | Genotype prefix. f2 blocks are extracted for exactly this command's populations, on the SNPs with data in all of them (with `--maxmiss 0`). Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Read f2 blocks saved by [extract-f2](extract-f2.md) instead of genotypes; every population must be in them. The SNP set is the one extract-f2 chose. |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--maxmiss` | number (0 to 1) | 0 | With `--input`: keep SNPs missing in at most this fraction of the populations. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `text`, `json`, `tsv` | `text` | `json`: left, right, estimated_rank, rank_tests (rank, chisq, dof, p), is_clade. `tsv`: one row per rank test. |

## Inputs and outputs

Genotypes or an f2 directory; the result on stdout:

```text
qpWave: Simulation.1, French_1, Yoruba_1
Outgroups: French_2, Yoruba_2, French_3, Yoruba_3

Rank tests:
  Rank  Chi-sq     df    P-value
     0  1240.0343     6  0.000000 *
     1    1.2836     2  0.526345

Estimated rank: 1
=> 2 independent ancestry streams
```

## Examples

```bash
# runs
popmix qpwave -i family6 --left Simulation.1,French_1,Yoruba_1 --right French_2,Yoruba_2,French_3,Yoruba_3
popmix qpwave -i family6 --left French_1,French_2 --right Yoruba_1,Yoruba_2,Yoruba_3 --format json
```

```bash
# not run: needs the AADR
popmix qpwave -i /references/aadr/v62.0_1240k_public --left French.DG,Sardinian.DG \
  --right Mbuti.DG,Han.DG,Karitiana.DG,Papuan.DG
```

## Reading the result

Each row tests "rank at most r" against full rank; `*` marks p < 0.05
(rejected). The estimated rank is the lowest rank not rejected at 0.05.
`=> Populations form a clade` when rank 0 is not rejected. A qpAdm model of
a target with n sources is the qpWave test of rank n - 1 on the target plus
the sources (the last row above equals the qpAdm chi-squared).

## Caveats

- The ranks tested go from 0 to (number of left populations - 2). If every
  test is rejected the estimated rank is one above the highest rank tested:
  full rank, as many streams as left populations.
- The same SNP-set note as [qpadm](qpadm.md) applies to `--input`.
- The regularization (fudge 1e-4) is fixed on the CLI.

## Reference

admixtools2 `qpwave()` (2.0.10), `drop_ranks`; Reich et al. 2012.
