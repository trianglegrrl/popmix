# popmix qpadm-rotate

Fit every qpAdm model of a rotation (admixtools2 `qpadm_rotate`): each
subset of the `--leftright` populations in turn is the sources, and the
rest of them, after `--rightfix`, the right populations. The target is the
same in every model. One TSV row per model.

## Synopsis

```text
popmix qpadm-rotate --target T --leftright P1,P2,... [--rightfix R1,R2,...] (--input PREFIX | --f2-dir DIR) [options] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--target` | text | required | The population modelled in every model. |
| `--leftright` | text (comma-separated) | required | Populations distributed between sources and right. Every non-empty subset small enough that the right set stays larger is a model's sources. |
| `--rightfix` | text (comma-separated) | none | Populations on the right in every model, first. |
| `--constrained`/`--unconstrained` | flag | `--unconstrained` | Constrain weights to [0, 1]. |
| `--alpha` | number (0 to 1) | 0.05 | Significance level for each model's `status`. |
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
| `--output`, `-o` | text (path) | none (stdout) | Write the table here. |

## Inputs and outputs

Genotypes or an f2 directory holding every population; a TSV table, one
row per model: `model`, `target`, `left` (the sources), `right`,
`n_sources`, `chisq`, `dof`, `p`, `feasible`, `status` (as
[qpadm](qpadm.md)'s verdict: `not tested`, `rejected`, `infeasible`, `not
rejected`), `weights` and `se` (comma-separated, in source order),
`n_blocks`, `n_snps`. Numbers have 12 significant digits.

```text
model	target	left	right	n_sources	chisq	dof	p	feasible	status	weights	se	n_blocks	n_snps
1	Simulation.1	French_1	French_3,Yoruba_3,Yoruba_2,Yoruba_1,French_2	1	1136.2114752	4	1.0715340757e-244	True	rejected	1	0	50	18721
4	Simulation.1	French_1,Yoruba_1	French_3,Yoruba_3,Yoruba_2,French_2	2	0.440087915487	2	0.802483521822	True	not rejected	0.212194016284,0.787805983716	0.0113003963273,0.0113003963273	50	18721
```

## Examples

```bash
# runs
popmix qpadm-rotate -i family6 --target Simulation.1 --leftright French_1,Yoruba_1,French_2 \
  --rightfix French_3,Yoruba_3,Yoruba_2
```

```bash
# not run: needs the AADR
popmix qpadm-rotate -i aadr_v66p1 --target Russia_Samara_EBA_Yamnaya.AG \
  --leftright Russia_EHG.AG,Iran_GanjDareh_N.AG,Turkey_Barcin_N.AG \
  --rightfix Mbuti.DG,Han.DG,Papuan.DG,Karitiana.DG -o rotation.tsv
```

## Reading the result

Models not rejected with every weight in [0, 1] (`status` `not rejected`)
are candidates, not answers: several can fit, and a model with fewer
sources that fits is simpler. The p-values are not corrected for testing
every model; popmix warns with the number of models tested, which a report
should give with any model chosen (Harney et al. 2021).

## Caveats

- Each model is a full qpAdm fit, so the number of models grows quickly
  with `--leftright`.
- The right set changes between models; a model's p depends on it.

## Reference

admixtools2 `qpadm_rotate` (Maier et al. 2023); the right set is
`--rightfix` followed by the rest of `--leftright`, as admixtools2 orders
it. Each fit is [qpadm](qpadm.md)'s, checked against admixtools2 in
`tests/golden/test_qpadm_models_admixtools2.py`.
