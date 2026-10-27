# popmix outgroup-f3

A matrix of outgroup f3(outgroup; A, B) for every pair of the populations:
the drift A and B share since they split from the outgroup. The standard
way to see which populations are closest to each other; with `--heatmap`,
also a picture.

## Synopsis

```text
popmix outgroup-f3 --outgroup O --populations A,B,... (--input PREFIX | --f2-dir DIR)
                   [--value est|se|z|p] [--order input|cluster] [--heatmap FILE] [options]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--outgroup` | text | required | The outgroup, the f3 target: a population equally distant from all the others (Mbuti or Yoruba for non-Africans). |
| `--populations` | text (comma-separated) | required | The matrix rows and columns. |
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

A symmetric matrix with the populations as row and column labels (TSV with
an empty first header cell, or JSON keyed by row then column). The diagonal
is f3(outgroup; A, A), the drift from the outgroup to A.

## Examples

```bash
# runs
popmix outgroup-f3 -i family --outgroup Yoruba --populations French,Simulation,Simulation.1
popmix outgroup-f3 -i family --outgroup Yoruba --populations French,Simulation,Simulation.1 \
  --value z --order cluster --format json -o og.json
```

```bash
# runs
popmix outgroup-f3 -i family --outgroup Yoruba --populations French,Simulation,Simulation.1 \
  --heatmap og.png
```

```bash
# not run: needs the AADR
popmix outgroup-f3 -i aadr_v66p1 --outgroup Mbuti.DG \
  --populations French.DG,Sardinian.DG,Han.DG,Papuan.DG,Karitiana.DG --order cluster --heatmap og.pdf
```

## Reading the result

Larger values: more shared drift, closer relationship. Compare off-diagonal
entries within a row (which B is closest to A); the diagonal is on another
scale. With `--value z` the entries are est / se.

## Caveats

- From genotypes this is `qp3pop --outgroup-mode` (not normalized, bias
  corrected); from `--f2-dir`, f3 from the blocks. Each pair uses its own
  SNPs from genotypes.
- `--heatmap` needs matplotlib (`pip install popmix[plot]`).

## Reference

admixtools2 `qp3pop(outgroupmode = TRUE)` arranged as a matrix (2.0.10).
