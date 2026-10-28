# popmix find-graphs

Search for the admixture graphs with given leaves and number of admixture
nodes that best fit the data; each candidate is fitted as
[qpgraph](qpgraph.md) fits it and ranked by score.

## Synopsis

```text
popmix find-graphs --populations A,B,... (--input PREFIX | --f2-dir DIR) [--nadmix N]
                   [--outgroup O] [--method stepwise|exhaustive] [--beam-width N] [--top N]
                   [--output-graph FILE] [--plot FILE] [--numstart N] [--seed S] [--format text|json]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--populations` | text (comma-separated) | required | The graphs' leaves; the first is the f3 base of every fit. |
| `--nadmix` | integer (at least 0) | 0 | Number of admixture nodes. |
| `--outgroup` | text | none | Only graphs where this population is a direct child of the root. |
| `--method` | `stepwise`, `exhaustive` | `stepwise` | `stepwise`: add populations, then admixture edges, one at a time, keeping the `--beam-width` best graphs at each step (fast, may miss the optimum); `exhaustive`: fit every graph, as admixtools2 `generate_all_graphs()` (in practice up to about 5 leaves with 1-2 admixtures). |
| `--plot` | text (path) | none | Also draw the best graph to this image file (needs `popmix[plot]`). |
| `--beam-width` | integer (at least 1) | 20 | Graphs kept after each stepwise step. |
| `--top` | integer (at least 1) | 5 | Graphs reported. |
| `--output-graph` | path (file) | none | Write the best graph here, as `popmix qpgraph --graph` reads it (with its score in a `#` comment). |
| `--numstart` | integer (at least 1) | 10 | Random starts for each graph's admixture weights. |
| `--seed` | integer | 1 | Seed for the random starts. |
| `--out-of-sample` / `--no-out-of-sample` | flag | `--out-of-sample` | Also report each graph shown out of sample (fitted on alternate jackknife blocks, scored on the others; two more fits per graph). |
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
| `--format` | `text`, `json` | `text` | `json`: n_evaluated and graphs (score, edges), best first. |

## Inputs and outputs

Genotypes or f2 blocks holding every population. Prints the number of
graphs fitted and the `--top` best, each with its score and edges; internal
nodes are named n1, n2, ...

## Examples

```bash
# runs
popmix find-graphs -i family --populations Yoruba,French,Simulation,Simulation.1 --nadmix 1 \
  --outgroup Yoruba --beam-width 3 --top 2 --output-graph best.graph
popmix qpgraph -i family --graph best.graph
popmix find-graphs -i family --populations Yoruba,French,Simulation,Simulation.1 \
  --method exhaustive --top 3 --format json
```

```bash
# not run: needs the AADR
popmix find-graphs -i aadr_v66p1 \
  --populations Mbuti.DG,Han.DG,French.DG,Karitiana.DG,Papuan.DG --nadmix 2 --outgroup Mbuti.DG \
  --method exhaustive --plot best.png --output-graph best.graph
```

## Reading the result

Lower score is better. Several graphs within a few score units fit about
equally well; do not over-read the top one: the best of many fits is
optimistic. The out-of-sample score of each graph shown (fitted on
alternate jackknife blocks, scored on the others) is a fairer comparison,
as is the penalized score across numbers of admixtures (`popmix qpgraph`). Compare graphs with resampling
in the Python API (`pm.compare_fits`), and check the best one with
[qpgraph](qpgraph.md) (worst residual |z|).

## Caveats

- The number of graphs grows very fast with leaves and admixtures; use
  `stepwise` beyond small cases, and `--outgroup` to cut the space.
- The search is popmix's own; admixtools2's evolutionary `find_graphs()`
  is not reproduced, but `exhaustive` gives the same graphs, order and
  scores as `generate_all_graphs()`.

## Reference

admixtools2 `generate_all_graphs()` and `qpgraph()` (2.0.10) for the fits;
the stepwise search is popmix's.
