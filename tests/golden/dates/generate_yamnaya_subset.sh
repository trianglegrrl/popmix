#!/usr/bin/env bash
#
# Make tests/golden/dates/yamnaya_chr21_22.{geno,snp,ind}: chromosomes 21
# and 22 of 74 AADR v66.p1 samples, EIGENSTRAT, written by convertf.
#
# Usage: CONVERTF=/path/to/convertf tests/golden/dates/generate_yamnaya_subset.sh
#
# French (31, the test population) = Sardinian (31) + Russia_Samara_EBA_Yamnaya
# (11 pseudo-haploid samples and the diploid genome I6729.DG). A source that
# mixes pseudo-haploid and diploid samples is what DATES-0004 is about: DATES
# counts every sample as diploid. DATES v753 references are then made with
# generate_dates_reference.sh yamnaya_chr21_22_q10.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AADR="/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB"
: "${CONVERTF:=$HOME/.local/share/popmix-references/admixtools-8.0.2/bin/convertf}"
OUT="$HERE/yamnaya_chr21_22"

YAMNAYA="I6729.DG I0231.AG.TW I0443.AG I0438.AG I0245.AG I1450.AG I6727.AG I6731.AG
I8291.AG I7479.AG I6728.AG I7677.AG"

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# every French and Sardinian sample, the listed Yamnaya ones; the rest Ignore
echo $YAMNAYA | tr ' ' '\n' > "$work/keep"
awk 'NR == FNR { keep[$1] = 1; next }
     { pop = ($3 == "French" || $3 == "Sardinian" || ($1 in keep)) ? $3 : "Ignore"; print $1, $2, pop }' \
    "$work/keep" "$AADR.ind" > "$work/sel.ind"

cat > "$work/par" <<EOF
genotypename: $AADR.geno
snpname: $AADR.snp
indivname: $work/sel.ind
outputformat: EIGENSTRAT
genooutfilename: $work/out.geno
snpoutfilename: $work/out.snp
indoutfilename: $work/out.ind
EOF
# convertf writes a profile (gmon.out) where it runs
(cd "$work" && "$CONVERTF" -p par > convertf.log)

grep -v ' Ignore$' "$work/out.ind" > "$OUT.ind"
awk '{ print ($2 == 21 || $2 == 22) }' "$work/out.snp" > "$work/rows"
paste -d '\t' "$work/rows" "$work/out.snp" | awk -F '\t' '$1 == 1 { print $2 }' | awk '{ print $1, $2, $3, $4, $5, $6 }' > "$OUT.snp"
paste -d '\t' "$work/rows" "$work/out.geno" | awk -F '\t' '$1 == 1 { print $2 }' > "$OUT.geno"
wc -l "$OUT.ind" "$OUT.snp"
