# popmix dstat-matrix

D- (or f4-) statistics for every pair X, Y of the populations, with one pair
fixed: D(X, Y; pop3, pop4) or D(pop1, pop2; X, Y). Shows at a glance which
populations share more with one side of a fixed pair.

## Synopsis

```text
popmix dstat-matrix --populations X1,X2,... (--pop3 C --pop4 D | --pop1 A --pop2 B)
                    (--input PREFIX | --f2-dir DIR --f4-mode) [--value ...] [--order ...] [--heatmap FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--populations` | text (comma-separated) | required | X and Y: the matrix rows and columns. |
| `--pop1` | text | none | Fixed pop1: entries are D(pop1, pop2; X, Y). Give with `--pop2`. |
| `--pop2` | text | none | Fixed pop2. |
| `--pop3` | text | none | Fixed pop3: entries are D(X, Y; pop3, pop4). Give with `--pop4`. |
| `--pop4` | text | none | Fixed pop4. |
| `--f4-mode` | flag | off | f4 instead of D; needed with `--f2-dir`. |
| `--value` | `est`, `se`, `z`, `p` | `est` | Which statistic fills the matrix: estimate, standard error, Z-score or p-value. |
| `--order` | `input`, `cluster` | `input` | Row and column order: as given in `--populations`, or clustered (average linkage on the rows) so similar populations are adjacent. |
| `--heatmap` | text (path) | none | Also draw the matrix as a heatmap to this image file (needs `popmix[plot]`). |
| `--input`, `-i` | text (path) | none | Genotype prefix. Each statistic is computed from the genotypes on the SNPs with data in its own populations, as admixtools2 qpdstat/qp3pop on genotypes. Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Use f2 blocks saved by [extract-f2](extract-f2.md) instead (f4 and f3 only, on extract-f2's SNP set). |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `tsv`, `json` | `tsv` | Table format: tab-separated with a header, or JSON. |
| `--output`, `-o` | text (path) | none (stdout) | Write the table to this file instead of stdout. |

## Inputs and outputs

An antisymmetric matrix: entry (X, Y) is the statistic with X in the first
free position and Y in the second; (Y, X) is its negative (for `est` and
`z`; `se` and `p` are symmetric); the diagonal is empty (NaN).

## Examples

```bash
# runs
popmix dstat-matrix -i family --populations French,Simulation,Simulation.1 --pop3 Yoruba --pop4 French
popmix dstat-matrix -i family --populations French,Simulation,Simulation.1 --pop1 Yoruba --pop2 French \
  --value z --format json
popmix extract-f2 -i family -o f2_family --populations French,Yoruba,Simulation,Simulation.1
popmix dstat-matrix --f2-dir f2_family --f4-mode --populations French,Simulation,Simulation.1 \
  --pop3 Yoruba --pop4 French -o f4.tsv
```

```bash
# not run: needs the AADR
popmix dstat-matrix -i aadr_v66p1 --pop1 Mbuti.DG --pop2 Han.DG \
  --populations French.DG,Sardinian.DG,Karitiana.DG,Papuan.DG --value z --heatmap d.png
```

## Reading the result

Row X, column Y with fixed (pop3, pop4): D(X, Y; pop3, pop4) > 0 when X is
closer to pop3 than Y is, or Y closer to pop4 (here French, with pop4
French, gives negative entries in its row). With fixed (pop1, pop2):
D(pop1, pop2; X, Y) > 0 when pop1 is closer to X, or pop2 to Y. Swapping
the fixed pair flips every entry; |z| > 3 is the usual threshold.

## Caveats

- Give exactly one fixed pair: `--pop1` and `--pop2`, or `--pop3` and
  `--pop4`.
- D needs genotypes (`--input`); with `--f2-dir` add `--f4-mode`, or the
  command stops with an error.

## Reference

admixtools2 `qpdstat()` (2.0.10) arranged as a matrix; Green et al. 2010,
Patterson et al. 2012.
