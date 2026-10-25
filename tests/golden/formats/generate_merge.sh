#!/usr/bin/env bash
#
# Reference for merging two datasets: AdmixTools 8 mergeit (build_convertf.sh)
# merges tests/golden/formats/merge/{a,b}.* (a SNP on the other strand,
# with and without the alleles swapped, an A/T SNP, a third allele, a
# monomorphic SNP on either side, two SNPs at one position) into
# tests/golden/formats/merge/reference.{geno,snp,ind}.
#
# Usage: tests/golden/formats/generate_merge.sh

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)/merge"
mergeit="${POPMIX_MERGEIT:-$HOME/.local/share/popmix-references/admixtools-8.0.2/bin/mergeit}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

cat > "$work/par" <<PAR
geno1: $here/a.geno
snp1: $here/a.snp
ind1: $here/a.ind
geno2: $here/b.geno
snp2: $here/b.snp
ind2: $here/b.ind
outputformat: EIGENSTRAT
genooutfilename: $here/reference.geno
snpoutfilename: $here/reference.snp
indoutfilename: $here/reference.ind
PAR
# mergeit writes a profile (gmon.out) where it runs
(cd "$work" && "$mergeit" -p par > mergeit.log) || { cat "$work/mergeit.log"; exit 1; }
grep -A6 "^Histogram" "$work/mergeit.log"
echo "wrote $here/reference.{geno,snp,ind}"
