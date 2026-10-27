# popmix info

Summarize a genotype dataset: its format, the number of samples, SNPs and
populations, and the population names.

## Synopsis

```text
popmix info PREFIX [--populations]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `PREFIX` | text (path) | required | Genotype prefix without extension (EIGENSTRAT, PACKEDANCESTRYMAP, TGENO, PLINK, PED, ANCESTRYMAP), or the full path of a `.vcf` / `.vcf.gz`. |
| `--populations`, `-p` | flag | off | List every population with its sample count. Without it, up to 10 populations are listed by name (the first 5 and a count when there are more). |

## Inputs and outputs

Reads the dataset's sample and SNP files (not the genotypes). Prints to
stdout:

```text
Dataset: family
Format: packedancestrymap
Samples: 54
SNPs: 21162
Populations: 4
  - French
  ...
```

## Examples

```bash
# runs
popmix info family
popmix info family --populations
popmix info small_tgeno
popmix info small_vcf.vcf
```

```bash
# not run: needs the AADR
popmix info aadr_v66p1 -p
```

## Reading the result

`Format` is the detected format (`packedancestrymap`, `tgeno`, `eigenstrat`,
`plink`, `ped`, `ancestrymap`, `vcf`). For a VCF every sample is its own
population, so `Populations` equals `Samples`.

## Caveats

- The order of the short population list depends on the format; use
  [populations](populations.md) for a sorted list with counts.
- `info` does not read genotypes, so it does not validate them; a damaged
  `.geno` shows up only when a command reads it.

## Reference

No single reference tool; it plays the role of looking at the `.ind` and
`.snp` files (or `wc -l`) before a convertf or admixtools2 run.
