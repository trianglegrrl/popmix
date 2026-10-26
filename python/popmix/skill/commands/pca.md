# popmix pca

Principal component analysis as EIGENSOFT smartpca: the PCA is computed on
one set of populations and other samples (typically ancient, with missing
data) are projected onto it, with outlier removal and Tracy-Widom tests of
the eigenvalues.

## Synopsis

```text
popmix pca --input PREFIX --output-prefix PREFIX [--pca-pops A,B,...] [--project-pops C,D,...]
           [--components N] [--outlier-iterations N] [--outlier-components N] [--outlier-sigma X]
           [--lsqproject] [--compat eigensoft-7.2.1] [--plot FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix. |
| `--pca-pops` | text (comma-separated) | none (every population not in `--project-pops`, except `Ignore`) | Populations the PCA is computed on (smartpca `poplistname`). |
| `--project-pops` | text (comma-separated) | none | Populations projected onto the PCA; they do not influence the axes. |
| `--components` | integer | 10 | Eigenvectors reported (smartpca `numoutevec`). |
| `--outlier-iterations` | integer | 5 | Rounds of outlier removal, 0 for none (smartpca `numoutlieriter`). |
| `--outlier-components` | integer | 10 | Eigenvectors examined for outliers (smartpca `numoutlierevec`). |
| `--outlier-sigma` | number | 6 | A PCA sample more than this many standard deviations out on one examined eigenvector is removed (smartpca `outliersigmathresh`). |
| `--lsqproject` | flag | off | Project by least squares on each sample's observed genotypes (smartpca `lsqproject: YES`); use it for samples with missing data, such as ancient ones. |
| `--compat` | `eigensoft-7.2.1` | none (corrected) | Reproduce smartpca exactly, including its Tracy-Widom p-values beyond the table (EIGENSOFT-0001). |
| `--output-prefix`, `-o` | text (path) | required | Write PREFIX.coords.tsv, PREFIX.eval, PREFIX.tracy_widom.tsv and PREFIX.outliers.tsv. |
| `--plot` | text (path) | none | Also plot PC1 against PC2 to this image file (needs `popmix[plot]`). |

## Inputs and outputs

| File | Content |
|------|---------|
| `PREFIX.coords.tsv` | One row per sample (PCA and projected, in `.ind` order): `id`, `pop`, `role` (`pca` or `project`), `PC1` ... |
| `PREFIX.eval` | Every eigenvalue, one per line, decreasing. |
| `PREFIX.tracy_widom.tsv` | `eigenvalue`, `twstat`, `p`, `neff` (effective number of markers) per eigenvalue. |
| `PREFIX.outliers.tsv` | Removed samples: `id`, `pop`, `iteration`, `pc` (from 1), `score` (standard deviations). |

stdout: `49 PCA samples, 5 projected, 21071 SNPs, 0 outliers removed`.

## Examples

```bash
# runs
popmix pca -i family --pca-pops French,Yoruba --project-pops Simulation,Simulation.1 --components 4 -o fam_pca
popmix pca -i family --pca-pops French,Yoruba --project-pops Simulation,Simulation.1 --lsqproject \
  --outlier-iterations 0 --compat eigensoft-7.2.1 -o fam_pca_lsq --plot fam_pca.png
```

```bash
# not run: needs the AADR
popmix --threads 8 pca -i /references/aadr/v62.0_1240k_public \
  --pca-pops French.DG,Sardinian.DG,Basque.DG,Russian.DG,Han.DG \
  --project-pops Russia_Samara_EBA_Yamnaya,Turkey_N --lsqproject -o west_eurasia
```

## Reading the result

Plot PC1 against PC2 from `coords.tsv`, colored by `pop`. The share of
variance of PC k is its eigenvalue over the sum in `.eval`. A Tracy-Widom
`p` below 0.01 marks a significant axis (structure, not noise). Outliers
are removed only from the PCA samples, before the final PCA; check
`outliers.tsv` to see whom it dropped, and set `--outlier-iterations 0` to
keep everyone.

## Caveats

- Eigenvector signs are arbitrary (they can differ from smartpca's).
- Pseudo-haploid samples are read as called (0/2), as smartpca does; that
  is why ancient samples should be projected, not used to build the axes.
- Only autosomal SNPs polymorphic in the PCA samples are used.
- The default corrects smartpca's Tracy-Widom p-values for statistics
  beyond its table (they were too large by sqrt(TW))
  ([EIGENSOFT-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#eigensoft-0001-tracy-widom-density-reported-as-the-p-value-beyond-the-table));
  `--compat eigensoft-7.2.1` reproduces them.

## Reference

EIGENSOFT 7.2.1 smartpca with `poplistname`, `lsqproject`, outlier removal
and `nostatslim` 10: eigenvalues to 1e-9 (Patterson, Price and Reich 2006).
