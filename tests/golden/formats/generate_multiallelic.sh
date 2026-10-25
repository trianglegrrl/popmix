#!/usr/bin/env bash
#
# Reference for reading a VCF with records that are not biallelic SNPs
# (tests/fixtures/small_vcf_multiallelic.vcf: a multiallelic site, indels,
# a site with no ALT allele, a symbolic allele): PLINK 2 keeps the
# biallelic SNPs (--snps-only just-acgt --max-alleles 2), then AdmixTools 8
# convertf (build_convertf.sh) writes them as EIGENSTRAT. The result goes to
# tests/golden/formats/multiallelic/reference.{geno,snp,ind}.
#
# Usage: tests/golden/formats/generate_multiallelic.sh

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
convertf="${POPMIX_CONVERTF:-$HOME/.local/share/popmix-references/admixtools-8.0.2/bin/convertf}"
out="$here/multiallelic"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

mkdir -p "$out"
plink2 --vcf "$repo/tests/fixtures/small_vcf_multiallelic.vcf" --snps-only just-acgt --max-alleles 2 \
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
