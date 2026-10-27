# popmix qpadm

qpAdm: model a target population as a mixture of source ("left")
populations, using reference ("right", outgroup) populations to tell the
sources apart. Reports the admixture weights and a chi-squared test of the
model.

## Synopsis

```text
popmix qpadm --target T --left S1,S2,... --right R1,R2,... (--input PREFIX | --f2-dir DIR)
             [--constrained | --unconstrained] [options] [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--target` | text | required | The population modeled. |
| `--left` | text (comma-separated) | required | Source populations, at least one. |
| `--right` | text (comma-separated) | required | Right (reference) populations, at least two; the first is the basis (R0) of the f4 matrix f4(target or source, ...; R0, Rj). |
| `--constrained` / `--unconstrained` | flag | `--unconstrained` | `--constrained` keeps every weight in [0, 1] (quadratic programming), as admixtools2 `constrained = TRUE`; the default allows any weights summing to 1. |
| `--input`, `-i` | text (path) | none | Genotype prefix. f2 blocks are extracted for exactly this command's populations, on the SNPs with data in all of them (with `--maxmiss 0`). Give this or `--f2-dir`. |
| `--f2-dir` | path (directory) | none | Read f2 blocks saved by [extract-f2](extract-f2.md) instead of genotypes; every population must be in them. The SNP set is the one extract-f2 chose. |
| `--block-size` | number | 0.05 | With `--input`: jackknife block length in Morgans. |
| `--maxmiss` | number (0 to 1) | 0 | With `--input`: keep SNPs missing in at most this fraction of the populations. |
| `--all-chromosomes` | flag | off | With `--input`: use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | With `--input`: which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `text`, `json`, `tsv` | `text` | `json`: target, sources, weights and se (by source), chisq, dof, p (null when untestable), feasible. `tsv`: one row per source with weight, chisq, dof, p. |

## Inputs and outputs

Genotypes or an f2 directory holding every population; the result on stdout:

```text
qpAdm: Simulation.1 = mix of French_1, Yoruba_1
Outgroups: French_2, Yoruba_2, French_3, Yoruba_3

Admixture weights:
  French_1: 0.1940 (19.4%)
  Yoruba_1: 0.8060 (80.6%)

Model fit: chi-sq=1.2836, df=2, p=0.5263
Feasible: True
=> Good model fit (p > 0.05)
```

## Examples

```bash
# runs
popmix qpadm -i family6 --target Simulation.1 --left French_1,Yoruba_1 \
  --right French_2,Yoruba_2,French_3,Yoruba_3
popmix qpadm -i family6 --target Simulation.1 --left French_1,Yoruba_1 \
  --right French_2,Yoruba_2,French_3,Yoruba_3 --format json
popmix qpadm -i family6 --target Simulation.1 --left French_1,Yoruba_1,Simulation \
  --right French_2,Yoruba_2,French_3,Yoruba_3 --constrained --format tsv
```

```bash
# not run: needs the AADR
popmix extract-f2 -i /references/aadr/v62.0_1240k_public -o f2_eu --maxmiss 0.1 \
  --populations Spain_EN,Turkey_N,Serbia_IronGates_Mesolithic,Mbuti.DG,Onge.DG,Papuan.DG,Han.DG,Karitiana.DG,Ami.DG
popmix qpadm --f2-dir f2_eu --target Spain_EN --left Turkey_N,Serbia_IronGates_Mesolithic \
  --right Mbuti.DG,Onge.DG,Papuan.DG,Han.DG,Karitiana.DG,Ami.DG
```

## Reading the result

- **Weights** sum to 1. **Feasible: True** when all lie in [0, 1]; an
  infeasible model (a negative weight) is rejected however good its p.
- **p > 0.05** (`=> Good model fit`): the target's f4 relations to the
  right populations are explained by the sources; the model is not
  rejected. p <= 0.05 rejects it. Report p with the weights.
- **dof** is the number of right populations minus the number of sources.
  With dof = 0 (as many right as left) there is nothing to test: p is `NA`
  in text and tsv, `null` in JSON, and the output says the model is not
  tested.
- Weight standard errors (block jackknife) are only in the JSON output
  (`se`).

## Caveats

- The right populations must be differently related to the sources;
  otherwise weights are unstable (large SEs) and the test has no power.
- `--input` extracts f2 for exactly the model's populations, on the SNPs
  they all share, which is admixtools2's `qpadm(f2_from_geno(prefix,
  pops))`. admixtools2's `qpadm(prefix, ...)` uses each f4's own SNPs
  instead, which differs slightly
  ([interface differences](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#interface-differences)).
  For many models on one SNP set, use [extract-f2](extract-f2.md) once and
  `--f2-dir`.
- Constrained fits with three or more sources work, where admixtools2
  stops with an error
  ([ADMIXTOOLS2-0002](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0002-constrained-qpadm-stops-with-three-or-more-sources)).
- The CLI reports the full model only; rank drops, popdrop (every subset of
  the sources), rotations and model sets are in the Python API
  (`pm.qpadm_p`, `pm.qpadm_rotate`, `pm.qpadm_multi`, `pm.qpadm_models`),
  with `compat="admixtools2-2.0.10"` for popdrop
  ([ADMIXTOOLS2-0001](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0001-popdrops-nested-path-compares-siblings)).
- Classic ADMIXTOOLS qpAdm uses a different SNP rule and covariance
  estimate, so its numbers differ slightly
  ([ADMIXTOOLS 7.0.2 (qpAdm)](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools-702-qpadm)).
- The regularization (fudge 1e-4) is fixed on the CLI.

## Reference

admixtools2 `qpadm()` (2.0.10): weights, SEs, chi-squared to 1e-6, p to
1e-4. Haak et al. 2015, Harney et al. 2021.
