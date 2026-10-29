#!/usr/bin/env bash
#
# Generate DATES v753 golden master reference data.
#
# Usage: DATES_BIN=/path/to/dates/bin tests/golden/dates/generate_dates_reference.sh [config...]
#
# DATES_BIN must contain the `dates` binary and its helper scripts
# (dates_jackknife, dates_wtjack, dates_expfit, dates_plot, calc_error,
# grabpars) plus `dowtjack`, which upstream only builds into src/.
# build_dates.sh builds all of them, patched, in
# ~/.cache/popmix-references/dates-v753/bin.
#
# The `dates` binary must be built with dates_numbins_fix.patch applied.
# Upstream v753 reads the uninitialized local `numbins` when sizing two
# buffers, so the exact (qbin: 0) mode crashes or not depending on stack
# contents. The patch only changes those buffer sizes; mesh-mode output is
# byte-identical to the unpatched binary.
#
# Each config writes tests/golden/dates/v753/<config>/ containing the
# parameter file, the LD table (.out), the full-data fit log (.flog) and
# curve (.fit), the per-chromosome jackknife estimates (.jin) and the final
# estimate (.jout).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="$HERE/v753"
# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache; an
# AADR config's directory gets a copy of its manifest.json (dataset.json),
# which records the release files' checksums.
AADR_DIR="${POPMIX_REFERENCE_CACHE:-$HOME/.cache/popmix-references}/aadr_v66p1"
AADR_PREFIX="$AADR_DIR/aadr_v66p1"

: "${DATES_BIN:?Set DATES_BIN to a directory containing patched dates and its helpers}"
export PATH="$DATES_BIN:$PATH"

# name | genotype prefix | source1 | source2 | test pop | maxdis (Morgans) | qbin
CONFIGS=(
    "chr20_22_q0|$HERE/family_chr20_22|French|Yoruba|Simulation.1|1.0|0"
    "chr20_22_q10|$HERE/family_chr20_22|French|Yoruba|Simulation.1|1.0|10"
    "simulation_q0|$HERE/example_data/family_packed|French|Yoruba|Simulation|1.0|0"
    "simulation_q10|$HERE/example_data/family_packed|French|Yoruba|Simulation|1.0|10"
    "aadr_french_q10|$AADR_PREFIX|Sardinian.DG|Russia_Samara_EBA_Yamnaya.AG|French.DG|0.5|10"
    # a source of mixed ploidy (DATES-0004); generate_yamnaya_subset.sh makes the input
    "yamnaya_chr21_22_q10|$HERE/yamnaya_chr21_22|Sardinian|Russia_Samara_EBA_Yamnaya|French|0.5|10"
)

run_config() {
    local name prefix src1 src2 test maxdis qbin
    IFS='|' read -r name prefix src1 src2 test maxdis qbin <<< "$1"

    local work
    work="$(mktemp -d)"
    trap 'rm -rf "$work"' RETURN

    echo "$src1 $src2 $test $test" > "$work/admixlist"
    cat > "$work/par" <<EOF
genotypename: $prefix.geno
snpname: $prefix.snp
indivname: $prefix.ind
admixlist: $work/admixlist
binsize: 0.001
maxdis: $maxdis
seed: 77
runmode: 1
chithresh: 0.0
mincount: 4
qbin: $qbin
runfit: YES
afffit: YES
lovalfit: 0.45
EOF

    echo "== $name"
    (cd "$work" && dates -p par > dates.log 2>&1)

    local res="$work/$test"
    for f in "$res/$test.out" "$res/$test.fit" "$res/$test.jin" "$res/$test.jout" \
             "$res/expfit_$test.flog"; do
        [[ -s "$f" ]] || { echo "missing DATES output: $f" >&2; tail -20 "$work/dates.log" >&2; return 1; }
    done

    local dest="$OUT_ROOT/$name"
    rm -rf "$dest"
    mkdir -p "$dest"
    # inputs in tests/golden/dates are named relative to the config directory
    sed -e "s|$work/admixlist|admixlist|" -e "s|: $HERE/|: ../../|" "$work/par" > "$dest/par"
    cp "$work/admixlist" "$dest/admixlist"
    for ext in out fit jin jout; do
        cp "$res/$test.$ext" "$dest/dates.$ext"
    done
    cp "$res/expfit_$test.flog" "$dest/dates.flog"
    if [[ "$prefix" == "$AADR_PREFIX" ]]; then
        cp "$AADR_DIR/manifest.json" "$dest/dataset.json"
    fi
    tail -3 "$dest/dates.jout"
}

selected=("$@")
for cfg in "${CONFIGS[@]}"; do
    name="${cfg%%|*}"
    if [[ ${#selected[@]} -eq 0 ]] || [[ " ${selected[*]} " == *" $name "* ]]; then
        run_config "$cfg"
    fi
done
