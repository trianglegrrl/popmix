# popmix convertf

Convert a genotype dataset from one file format to another (VCF to PLINK,
PLINK to EIGENSTRAT, anything to PACKEDANCESTRYMAP, ...).

## Synopsis

```text
popmix convertf --input PREFIX --output PREFIX --output-format FORMAT [--input-format FORMAT]
popmix convertf --par FILE.par
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--par`, `-p` | text (path) | none | Run a convertf parameter file instead of the options (see below); give it alone. |
| `--input`, `-i` | text (path) | required (unless `--par`) | Input genotype prefix, or the full path of a `.vcf` / `.vcf.gz`. |
| `--output`, `-o` | text (path) | required (unless `--par`) | Output prefix; the files are written as PREFIX plus the format's extensions, overwriting any that exist. |
| `--input-format`, `-if` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno`, `vcf` | none (auto-detected) | Input format, when detection from the files is not wanted. |
| `--output-format`, `-of` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno`, `vcf` | required (unless `--par`) | Output format. `plink` and `packedped` are the same (`.bed/.bim/.fam`). `vcf` writes `OUTPUT.vcf` (or `OUTPUT` itself if it ends in `.vcf`, or `.vcf.gz`, gzip-compressed), as PLINK 2's `--export vcf`: REF is the second allele, ALT the counted first one, and males' homozygous chrX calls are haploid. |
| `--keep-inds` | text (path) | none | Use only the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--exclude-inds` | text (path) | none | Never use the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--relabel` | text (path) | none | Give samples new populations: `ID POPULATION` per line, as admixtools2 `inds`/`pops`. Population names given to the command are those after relabelling. |
| `--relabel-pops` | text (path) | none | Rename populations: `OLD NEW` per line, applied after `--relabel` (merge two populations by giving them one name). |
| `--populations` | text (comma-separated) | none (all) | Keep only samples of these populations, by their labels after `--relabel` and `--relabel-pops` (convertf `poplistname`). |
| `--keep-snps` | text (path) | none | Keep only the SNP IDs in this file (one per line, or a `.snp`/`.bim` file). |
| `--exclude-snps` | text (path) | none | Leave out the SNP IDs in this file. |
| `--chromosomes` | text (comma-separated) | none (all) | Keep only SNPs on these chromosomes (`1`-`22`, `X`, `Y`, `MT`). |
| `--region` | text (`CHROM:START-END`) | none | Keep only SNPs in this region, base pairs, both ends included. |

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

## Parameter files

`popmix convertf -p file.par` runs a convertf parameter file, so a pipeline
can swap popmix in. Keys: `genotypename`, `snpname`, `indivname` (inputs;
their names need not share a stem); `outputformat` (`ANCESTRYMAP`,
`EIGENSTRAT`, `PED`, `PACKEDPED`, `PACKEDANCESTRYMAP` (the default),
`TRANSPOSE_PACKED`); `genotypeoutname`, `snpoutname`, `indivoutname` (or
`genooutfilename`, `snpoutfilename`, `indoutfilename`); `poplistname` (a
file of populations to keep); `badsnpname` (SNPs to leave out; ones the data
lacks are ignored); `chrom`, with `lopos` and `hipos` for a region. A value
may start with a name defined earlier and `/` (`DIR: /data`,
`genotypename: DIR/x.geno`); `#` starts a comment. Any other key is an
error, not silently ignored. Each key is checked against AdmixTools 8
convertf on the same file.

```text
DIR: /data/aadr
genotypename: DIR/v66.p1.geno
snpname: DIR/v66.p1.snp
indivname: DIR/v66.p1.ind
outputformat: EIGENSTRAT
poplistname: pops.txt
genotypeoutname: subset.geno
snpoutname: subset.snp
indivoutname: subset.ind
```

## Reading the result

`Converted N SNPs x M samples` counts what was written. From a VCF, records
that are not biallelic SNPs (indels, multiallelic, symbolic alleles) are
skipped, so N can be smaller than the number of records.

Selection runs in this order: samples kept or left out by ID, given new
populations by ID (`--relabel`), populations renamed (`--relabel-pops`),
then samples chosen by population (`--populations`); SNPs are kept by ID,
chromosome and region together. Samples and SNPs keep their file order. An
ID or population named that is not in the data is an error, so a typo
never passes silently.

## Caveats

- VCF input: each sample becomes its own population (group them with
  `--relabel`, a file of `ID POPULATION` lines); haploid calls count as homozygous;
  FILTER is ignored, so filter the VCF first
  ([PLINK 2: reading VCF](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#plink-2-reading-vcf)).
- A `.snp` line needs six columns and single-character alleles; a longer
  allele is read as its first base, as convertf does; a genetic position
  that is not a number is an error, where convertf reads it as 0
  ([convertf: reading .snp and .bim](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#convertf-admixtools-802-reading-snp-and-bim)).
- PLINK chromosome codes X, Y, XY and MT are read; they are not autosomes.
- PACKEDANCESTRYMAP is the format to keep large data in: `extract-f2` and
  the f-statistic commands stream it (and PLINK `.bed` and text
  EIGENSTRAT); from TGENO they read the selected populations' samples, and
  other formats whole.
- `convertf` itself reads packed, PLINK, text EIGENSTRAT, TGENO and VCF
  input a part at a time, for every output format and subset: the whole of
  AADR v66.p1 (TGENO, 23,089 samples x 1.23M SNPs) converts to
  PACKEDANCESTRYMAP holding about 0.6 GB. PED and ANCESTRYMAP input are
  read whole.
- A PLINK sample is named by its IID and its population is the FID, where
  convertf names it `FID:IID`
  ([convertf: reading .snp and .bim](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#convertf-admixtools-802-reading-snp-and-bim)).
- Sex is carried between formats: M/F/U in `.ind`, 1/2/0 in `.fam` and `.ped`;
  a VCF has none, so its samples are U.

## Reference

AdmixTools `convertf` (8.0.2; outputs checked file for file) and, for VCF
input, `plink2 --vcf ... --snps-only just-acgt --max-alleles 2`.
