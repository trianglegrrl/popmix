#!/usr/bin/env bash
#
# Reference for VCF output (TASKS.md P1): AdmixTools 8 convertf
# (build_convertf.sh) writes tests/golden/formats/vcf_out/source.* (het,
# missing, a monomorphic SNP with an unknown allele, chromosome X) as
# PACKEDPED, and PLINK 2 exports that as VCF into vcf_out/reference.vcf.
# Neither step is popmix's.
#
# Usage: tests/golden/formats/generate_vcf_out.sh

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)/vcf_out"
convertf="${POPMIX_CONVERTF:-$HOME/.local/share/popmix-references/admixtools-8.0.2/bin/convertf}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

cat > "$work/par" <<PAR
genotypename: $here/source.geno
snpname: $here/source.snp
indivname: $here/source.ind
outputformat: PACKEDPED
genotypeoutname: $work/bed.bed
snpoutname: $work/bed.bim
indivoutname: $work/bed.fam
PAR
(cd "$work" && "$convertf" -p par > convertf.log) || { cat "$work/convertf.log"; exit 1; }
plink2 --bfile "$work/bed" --export vcf id-paste=iid --out "$work/reference" > /dev/null
grep -v "^##fileDate\|^##source" "$work/reference.vcf" > "$here/reference.vcf"
echo "wrote $here/reference.vcf"
