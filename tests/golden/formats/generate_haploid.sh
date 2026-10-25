#!/usr/bin/env bash
#
# Reference for reading a VCF with haploid calls and a record that fails its
# FILTER (tests/fixtures/small_vcf_haploid.vcf): PLINK 2 reads a haploid call
# as homozygous and keeps every record whatever its FILTER, then AdmixTools 8
# convertf (build_convertf.sh) writes them as EIGENSTRAT. The result goes to
# tests/golden/formats/haploid/reference.{geno,snp,ind}.
#
# Usage: tests/golden/formats/generate_haploid.sh

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
convertf="${POPMIX_CONVERTF:-$HOME/.local/share/popmix-references/admixtools-8.0.2/bin/convertf}"
out="$here/haploid"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

mkdir -p "$out"
plink2 --vcf "$repo/tests/fixtures/small_vcf_haploid.vcf" \
  --make-bed --allow-extra-chr --out "$work/bed" > /dev/null
cat > "$work/par" <<EOF
genotypename: $work/bed.bed
snpname: $work/bed.bim
indivname: $work/bed.fam
outputformat: EIGENSTRAT
genooutfilename: $out/reference.geno
snpoutfilename: $out/reference.snp
indoutfilename: $out/reference.ind
EOF
# convertf writes a profile (gmon.out) where it runs
(cd "$work" && "$convertf" -p par > /dev/null)
echo "wrote $out/reference.{geno,snp,ind}"
