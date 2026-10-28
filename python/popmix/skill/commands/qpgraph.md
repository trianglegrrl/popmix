# popmix qpgraph

Fit an admixture graph to the data, as admixtools2 `qpgraph()`: drift
lengths and admixture weights are chosen to best predict f3(base; i, j) for
every pair of the graph's leaves, and the fit is scored.

## Synopsis

```text
popmix qpgraph --graph FILE (--input PREFIX | --f2-dir DIR) [--admix-weights W1,...]
               [--numstart N] [--seed S] [--diag X] [--diag-f3 X] [--unconstrained]
               [--plot FILE] [--plot-fit FILE] [--dot FILE] [--compat admixtools2-2.0.10] [options]
               [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--graph` | path (file) | required | The graph: one edge per line, `from to`, whitespace-separated; blank lines and text after `#` are ignored. Nodes without outgoing edges are the leaves, which must be population names; the root is the node without incoming edges (exactly one); a node with two incoming edges is an admixture node. The first leaf in the file is the f3 base population. |
| `--admix-weights` | text (comma-separated numbers) | none (optimized) | Fit at these weights instead of optimizing them: one per admixture node, in order of first appearance, each the weight of the node's incoming edge whose source node appears first in the file. |
| `--numstart` | integer (at least 1) | 10 | Random starts for the admixture-weight optimization; the best is reported. |
| `--seed` | integer | 1 | Seed for the random starts, as R's `set.seed`: with the same seed the starts are admixtools2's. |
| `--diag` | number | 1e-4 | Regularization of the edge-length fit (admixtools2 `diag`). |
| `--diag-f3` | number | 1e-5 | Regularization of the f3 covariance (admixtools2 `diag_f3`). |
| `--unconstrained` | flag | off | Allow negative drift lengths and weights outside [0, 1] (admixtools2 `constrained = FALSE`). |
| `--plot` | text (path) | none | Also draw the fitted graph to this image file (needs `popmix[plot]`). |
| `--plot-fit` | text (path) | none | Also plot the residual z of each f3 to this image file (needs `popmix[plot]`). |
| `--out-of-sample` / `--no-out-of-sample` | flag | `--out-of-sample` | Also report the out-of-sample score: the graph fitted on alternate jackknife blocks and scored on the others, both ways, averaged (two more fits; `pm.qpgraph_out_of_sample`). |
| `--dot` | text (path) | none | Also write the fitted graph in Graphviz DOT (drift x 1000 and weights in % as labels; names quoted). |
| `--compat` | `admixtools2-2.0.10` | none (corrected) | Reproduce admixtools2: stop the weight fit where its L-BFGS-B stalls before a zero-bound kink (ADMIXTOOLS2-0009), and with `--dot` write its DOT text exactly, which strips "." and "-" from names on edges only (ADMIXTOOLS2-0004). |
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
| `--format` | `text`, `json`, `tsv` | `text` | `text`: score, worst residual and the edges; `json`: pops, score, edges (from, to, type, weight, low, high), f3 (pop1, pop2, pop3, est, se, fit, diff, z, p) and opt (each start's weights, score, evaluations, convergence); `tsv`: the edge table. |

## Inputs and outputs

A graph file (the format [find-graphs](find-graphs.md) `--output-graph`
writes) and genotypes or f2 blocks holding every leaf. Output:

```text
qpGraph: 3 populations, 6 edges
f3 base population: Yoruba
Score: 0.0714
Worst f3 residual: f3(Yoruba; Simulation.1, Simulation.1) z = -0.21

  from             to               type   weight
  R                Yoruba           edge   0.000000
  R                n                edge   0.056653
  n                French           edge   0.000000
  n                adm              admix  0.189144
  R                adm              admix  0.810856
  adm              Simulation.1     edge   0.000000
```

## Examples

```bash
# runs
cat > admixed.graph <<'GRAPH'
# Simulation.1 = French-like + Yoruba-like
R    Yoruba
R    n
n    French
n    adm
R    adm
adm  Simulation.1
GRAPH
popmix qpgraph -i family --graph admixed.graph
popmix qpgraph -i family --graph admixed.graph --admix-weights 0.2 --format tsv
popmix qpgraph -i family --graph admixed.graph --numstart 3 --seed 7 --dot admixed.dot --format json > fit.json
```

```bash
# runs
printf 'R Yoruba\nR a\na French\na b\nb Simulation\nb Simulation.1\n' > tree.graph
popmix qpgraph -i family --graph tree.graph --plot tree.png --plot-fit tree_fit.png
```

```bash
# not run: needs the AADR
popmix extract-f2 -i aadr_v66p1 -o f2_graph \
  --populations Mbuti.DG,Han.DG,French.DG,Karitiana.DG,Papuan.DG
popmix qpgraph --f2-dir f2_graph --graph my.graph --numstart 50 --plot my.pdf
```

## Reading the result

- **Score**: the weighted squared error of the fitted against the observed
  f3 (admixtools2's score); lower is better, and it is comparable between
  graphs on the same populations and f2 blocks. It falls with every
  parameter added, so a lower score alone does not favour a graph with
  more admixtures.
- **Parameters** (drift lengths and admixture weights) against the f3
  statistics fitted: with 0 dof or fewer the graph can fit anything. The
  **penalized score** adds twice the parameters (AIC-style), and the
  **out-of-sample score** fits on alternate jackknife blocks and scores on
  the others: neither falls just because a graph has more parameters.
- A warning that starts reaching the best score disagree on the weights
  means the f3 statistics do not determine them: report the range, not
  one value.
- **Worst f3 residual**: |z| over about 3 means the graph misses that
  relation; `--plot-fit` shows all of them.
- Edges: `edge` weights are drift lengths (f2 units); `admix` weights are
  admixture proportions, the two into one node summing to 1. `low` and
  `high` (JSON, TSV) give the range of an admixture weight over the random
  starts; they are empty for drift edges and with `--admix-weights`.
- Drift lengths of 0 are common: edges the data cannot resolve, or
  constrained at the 0 bound.

## Caveats

- The f3 base is the first leaf in the graph file; reorder the lines to
  change it. A graph with two roots, or a leaf that is not a population, is
  an error.
- Admixture weights use the same random starts and L-BFGS-B as admixtools2
  and match its score to 3e-12; on flat likelihood ridges a start can end
  in a neighbouring optimum
  ([qpGraph matched behaviour](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#matched-not-corrected-1)).
  Use more `--numstart` for graphs with several admixtures.
- `--dot` keeps names as they are, quoted; admixtools2 `write_dot()`
  removes "." and "-"
  ([ADMIXTOOLS2-0004](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0004-write_dot-renames-nodes-on-edges-but-not-where-it-declares-them)).
- Resampling (jackknife, bootstrap, comparing two graphs) is in the Python
  API (`pm.qpgraph_resample_snps`, `pm.compare_fits`).

## Reference

admixtools2 `qpgraph()` (2.0.10), every start's score to 3e-12; Patterson et
al. 2012, Maier et al. 2023.
