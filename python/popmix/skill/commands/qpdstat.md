# popmix qpdstat

f4- or D-statistics f(pop1, pop2; pop3, pop4) for many population
combinations at once, as admixtools2 `qpdstat()`. Use it to scan for gene
flow (which populations share extra drift with a test population).

## Synopsis

```text
popmix qpdstat --pop1 A[,...] [--pop2 B,... --pop3 C,... --pop4 D,...] (--input PREFIX | --f2-dir DIR)
               [--d-mode] [options] [--format tsv|json] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--pop1` | text (comma-separated) | required | First positions. Alone: every distinct statistic on these populations (each 4-subset in its three splits). |
| `--pop2` | text (comma-separated) | none | Second positions. `--pop2`, `--pop3` and `--pop4` go together: every combination of one population from each of the four lists (pop1 varying slowest), repeats dropped. |
| `--pop3` | text (comma-separated) | none | Third positions. |
| `--pop4` | text (comma-separated) | none | Fourth positions. |
| `--d-mode` | flag | off | D-statistics, (a-b)(c-d) / ((a+b-2ab)(c+d-2cd)) summed over SNPs, instead of f4. Needs `--input`: with `--f2-dir` it is an error. |
| `--input`, `-i` | text (path) | none | Genotype prefix. Each statistic is computed from the genotypes on the SNPs with data in its own populations, as admixtools2 qpdstat/qp3pop on genotypes. Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Use f2 blocks saved by [extract-f2](extract-f2.md) instead (f4 and f3 only, on extract-f2's SNP set). |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `tsv`, `json` | `tsv` | Table format: tab-separated with a header, or JSON. |
| `--output`, `-o` | text (path) | none (stdout) | Write the table to this file instead of stdout. |

## Inputs and outputs

A table, one row per statistic: `pop1`, `pop2`, `pop3`, `pop4`, `est`,
`se`, `z`, `p` and, from genotypes, `n` (SNPs used).

```text
pop1	pop2	pop3	pop4	est	se	z	p	n
Yoruba	French	Simulation	French	0.0464031209408	0.00252332505691	18.389672315	1.58937684565e-75	21162
```

## Examples

```bash
# runs
popmix qpdstat -i family --pop1 Yoruba --pop2 French --pop3 Simulation,Simulation.1 --pop4 French
popmix qpdstat -i family --pop1 Yoruba,French,Simulation,Simulation.1 --d-mode -o dstats.tsv
popmix extract-f2 -i family -o f2_family --populations French,Yoruba,Simulation,Simulation.1
popmix qpdstat --f2-dir f2_family --pop1 Yoruba,French,Simulation,Simulation.1 --format json
```

```bash
# not run: needs the AADR
popmix qpdstat -i aadr_v66p1 --d-mode \
  --pop1 Mbuti.DG --pop2 French.DG,Sardinian.DG,Basque.DG --pop3 Han.DG --pop4 Karitiana.DG
```

## Reading the result

`z` is est / se from the block jackknife and `p` its two-sided normal
p-value; |z| > 3 is the usual threshold, and with many statistics correct
for the number of tests. Sign: positive when pop1-pop2 differences
correlate with pop3-pop4 differences (pop2 closer to pop4, or pop1 to pop3).

## Caveats

- From genotypes each statistic uses the SNPs with data in its own four
  populations (as admixtools2 qpdstat on genotypes); from `--f2-dir` all use
  extract-f2's SNP set, and `n` is not reported.
- D needs genotypes; admixtools2 silently returns f4 when given f2 blocks
  and asked for D, popmix stops with an error
  ([interface differences](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#interface-differences)).
- A statistic that cannot be estimated reports NaN and `n` 0, where
  admixtools2 returns NA everywhere.
- `--pop2` to `--pop4` must all be given, or none.

## Reference

admixtools2 `qpdstat()` (2.0.10): combinations and order as admixtools2,
estimates to rtol 1e-10, SEs 1e-6. Classic ADMIXTOOLS qpDstat.
