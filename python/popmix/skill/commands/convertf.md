# popmix convertf

Convert a genotype dataset from one file format to another (VCF to PLINK,
PLINK to EIGENSTRAT, anything to PACKEDANCESTRYMAP, ...).

## Synopsis

```text
popmix convertf --input PREFIX --output PREFIX --output-format FORMAT [--input-format FORMAT]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Input genotype prefix, or the full path of a `.vcf` / `.vcf.gz`. |
| `--output`, `-o` | text (path) | required | Output prefix; the files are written as PREFIX plus the format's extensions, overwriting any that exist. |
| `--input-format`, `-if` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno`, `vcf` | none (auto-detected) | Input format, when detection from the files is not wanted. |
| `--output-format`, `-of` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno` | required | Output format. `plink` and `packedped` are the same (`.bed/.bim/.fam`). VCF cannot be written. |

## Inputs and outputs

| Format | Files written |
|--------|---------------|
| `eigenstrat` | `.geno` (text), `.snp`, `.ind` |
| `packedancestrymap` | `.geno` (packed, one record per SNP), `.snp`, `.ind` |
| `tgeno` | `.geno` (packed, one record per sample), `.snp`, `.ind` |
| `plink` / `packedped` | `.bed`, `.bim`, `.fam` (population as family ID; genetic positions in Morgans) |
| `ped` | `.ped`, `.map`, `.pedind` |
| `ancestrymap` | `.ancestrymapgeno`, `.snp`, `.ind` |

Prints the detected input format and `Converted N SNPs x M samples`.

## Examples

```bash
# runs
popmix convertf -i small_eigenstrat -o small_plink -of plink
popmix convertf -i small_vcf.vcf -o from_vcf -of eigenstrat
popmix convertf -i family -o family_text -of eigenstrat
popmix convertf -i small_tgeno -if tgeno -o small_packed -of packedancestrymap
popmix info small_plink
```

```bash
# not run: needs large data
# convert once to PACKEDANCESTRYMAP, which extract-f2 streams instead of loading
popmix convertf -i my_samples.vcf.gz -o my_samples -of packedancestrymap
```

## Reading the result

`Converted N SNPs x M samples` counts what was written. From a VCF, records
that are not biallelic SNPs (indels, multiallelic, symbolic alleles) are
skipped, so N can be smaller than the number of records.

## Caveats

- VCF input: each sample becomes its own population (edit the `.ind` or
  `.fam` afterwards to group samples); haploid calls count as homozygous;
  FILTER is ignored, so filter the VCF first
  ([PLINK 2: reading VCF](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#plink-2-reading-vcf)).
- A `.snp` line needs six columns and single-character alleles; a longer
  allele is read as its first base, as convertf does; a genetic position
  that is not a number is an error, where convertf reads it as 0
  ([convertf: reading .snp and .bim](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#convertf-admixtools-802-reading-snp-and-bim)).
- PLINK chromosome codes X, Y, XY and MT are read; they are not autosomes.
- PACKEDANCESTRYMAP is the format to keep large data in: `extract-f2` and
  the f-statistic commands stream it, and read other formats into memory.
- Sex is carried between formats: M/F/U in `.ind`, 1/2/0 in `.fam` and `.ped`;
  a VCF has none, so its samples are U.

## Reference

AdmixTools `convertf` (8.0.2; outputs checked file for file) and, for VCF
input, `plink2 --vcf ... --snps-only just-acgt --max-alleles 2`.
