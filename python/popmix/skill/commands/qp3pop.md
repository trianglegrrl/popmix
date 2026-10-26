# popmix qp3pop

f3-statistics f3(pop1; pop2, pop3) for many combinations, as admixtools2
`qp3pop()`: admixture tests (pop1 the target) or, with `--outgroup-mode`,
outgroup f3 (pop1 the outgroup).

## Synopsis

```text
popmix qp3pop --pop1 A[,...] [--pop2 B,... --pop3 C,...] (--input PREFIX | --f2-dir DIR)
              [--outgroup-mode] [--no-bias-correction] [options] [--format tsv|json] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--pop1` | text (comma-separated) | required | Targets. Alone: every distinct f3 on these populations (each 3-subset in its three rotations). |
| `--pop2` | text (comma-separated) | none | Second positions; give with `--pop3`: every combination of one population from each list, repeats dropped. |
| `--pop3` | text (comma-separated) | none | Third positions. |
| `--outgroup-mode` | flag | off | With `--input`: do not divide f3 by the target's heterozygosity (admixtools2 `outgroupmode = TRUE`), for outgroup f3. |
| `--no-bias-correction` | flag | off | With `--input`: do not subtract the target's sampling bias (admixtools2 `apply_corr = FALSE`). |
| `--input`, `-i` | text (path) | none | Genotype prefix. Each statistic is computed from the genotypes on the SNPs with data in its own populations, as admixtools2 qpdstat/qp3pop on genotypes. Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Use f2 blocks saved by [extract-f2](extract-f2.md) instead (f4 and f3 only, on extract-f2's SNP set). |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `tsv`, `json` | `tsv` | Table format: tab-separated with a header, or JSON. |
| `--output`, `-o` | text (path) | none (stdout) | Write the table to this file instead of stdout. |

## Inputs and outputs

A table: `pop1`, `pop2`, `pop3`, `est`, `se`, `z`, `p` and, from genotypes,
`n` (SNPs used).

## Examples

```bash
# runs
popmix qp3pop -i family --pop1 Simulation.1,Simulation --pop2 French --pop3 Yoruba
popmix qp3pop -i family --pop1 Yoruba,French,Simulation.1 --format json
popmix qp3pop -i family --pop1 Yoruba --pop2 French --pop3 Simulation,Simulation.1 --outgroup-mode
```

```bash
# not run: needs the AADR
popmix qp3pop -i /references/aadr/v62.0_1240k_public --pop1 Uygur.DG,Hazara.DG \
  --pop2 Han.DG --pop3 French.DG,Sardinian.DG
```

## Reading the result

For admixture tests, z < -3 means pop1 is admixed between populations
related to pop2 and pop3. From genotypes the estimate is normalized by the
target's heterozygosity (admixtools2's default), so it is larger in
magnitude than the f3 of [f3](f3.md) (-0.0275 against -0.0088 for
Simulation.1 here) with a similar z. In outgroup mode, larger values mean
more drift shared by pop2 and pop3 since they split from pop1.

## Caveats

- `--outgroup-mode` and `--no-bias-correction` apply only with `--input`;
  f3 from `--f2-dir` is never normalized, and the flags are ignored there.
- A target that is a single pseudo-haploid sample has no bias correction
  (0/0): its f3 is NaN with `n` 0 (admixtools2 returns NA).
- From f2 blocks the combinations keep their order; admixtools2 sorts them
  ([interface differences](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#interface-differences)).

## Reference

admixtools2 `qp3pop()` (2.0.10), rtol 1e-10; classic ADMIXTOOLS qp3Pop.
