# popmix admixture

Maximum-likelihood ancestry proportions of each sample from K ancestral
populations, as ADMIXTURE (Alexander et al. 2009): unsupervised,
supervised (reference populations fixed) or projection onto an earlier
run's allele frequencies, with cross-validation to choose K.

## Synopsis

```text
popmix admixture --input PREFIX -k K --output-prefix PREFIX [--populations A,B,...]
                 [--supervised R1,...,RK | --fixed-p FILE.P] [--cv FOLDS]
                 [--chromosomes 1,2,...] [--seed S] [--tolerance X] [--plot FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix. |
| `--populations` | text (comma-separated) | none (all populations except `Ignore`) | Populations whose samples are used. |
| `-k` | integer | required | Number of ancestral populations (clusters). |
| `--supervised` | text (comma-separated) | none (unsupervised) | K reference populations: their samples keep Q fixed at their own cluster (ADMIXTURE `--supervised`); they must be among the samples used. |
| `--fixed-p` | path (file) | none | Projection (ADMIXTURE `-P`): a `.P` file (SNPs x K) from an earlier run on the same SNPs; only Q is estimated. |
| `--cv` | integer | none | Folds of cross-validation (ADMIXTURE `--cv`); prints the CV error. |
| `--chromosomes` | text (comma-separated) | none (every chromosome in the data) | Use only SNPs on these chromosomes. |
| `--seed` | integer | 1 | Seed for the random start (and the folds). |
| `--tolerance` | number | 1e-4 | Stop when the log-likelihood improves by less than this (ADMIXTURE `-C`). |
| `--output-prefix`, `-o` | text (path) | required | Write PREFIX.K.Q.tsv (id, pop, Q1..QK), PREFIX.K.Q and PREFIX.K.P (ADMIXTURE's formats). |
| `--plot` | text (path) | none | Also draw a bar plot of Q to this image file (needs `popmix[plot]`). |

## Inputs and outputs

| File | Content |
|------|---------|
| `PREFIX.K.Q.tsv` | One row per sample: `id`, `pop`, `Q1` ... `QK`. |
| `PREFIX.K.Q` | ADMIXTURE's `.Q`: the proportions only, 6 decimals. |
| `PREFIX.K.P` | ADMIXTURE's `.P`: SNPs x K frequencies of the allele the genotype file does not count. |

stdout: `Loglikelihood: -274582.780506 (8 iterations)`, `not converged` if
the iteration limit was hit, and `CV error (K=2): ...` with `--cv`.

## Examples

```bash
# runs
popmix admixture -i family -k 2 --chromosomes 22 -o fam
popmix admixture -i family -k 2 --chromosomes 22 --supervised French,Yoruba -o fam_sup
popmix admixture -i family -k 2 --chromosomes 22 --populations Simulation,Simulation.1 \
  --fixed-p fam.2.P -o fam_proj
```

```bash
# runs
popmix admixture -i family -k 2 --chromosomes 22 --populations French,Yoruba --cv 2 --seed 3 \
  --tolerance 1e-3 -o fam_cv --plot fam_cv.png
```

```bash
# not run: needs large data
for k in 2 3 4 5; do popmix --threads 8 admixture -i panel -k $k --cv 5 -o panel; done
```

## Reading the result

Each row of Q sums to 1. Unsupervised cluster numbers are arbitrary (they
change with the seed); name them by the populations that load on them.
Supervised clusters are numbered in the order their reference population's
first sample appears in the `.ind` (here Yoruba before French, so Q1 is
Yoruba), not in the order of `--supervised`. Choose K by the lowest CV
error. A higher log-likelihood is a better fit for the same K and data.

## Caveats

- Thin linked SNPs first (ADMIXTURE assumes independent markers), and
  remove close relatives ([kinship](kinship.md)).
- Pseudo-haploid samples are read as called (homozygous), as ADMIXTURE
  does.
- `--fixed-p` needs the P of a run on exactly the same SNPs: SNPs with no
  genotype in the samples used are left out, so a different
  `--populations` set can change the SNPs and the projection then fails.
- Different random starts can end at different local maxima for larger K;
  run several seeds and keep the best log-likelihood
  ([ADMIXTURE optimizer differences](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#differences-in-the-optimizer-not-deviations)).
- CV folds come from popmix's random numbers, not ADMIXTURE's, so CV errors
  agree with ADMIXTURE's only as closely as two random splits.
- Q and P are kept within [1e-5, 1 - 1e-5], so supervised references show
  0.99999.

## Reference

ADMIXTURE 1.3.0: the same likelihood, Q to 1e-4 where the maximum is well
defined, projection to its printed precision.
