#!/usr/bin/env bash
#
# Build classic ADMIXTOOLS 7.0.2 (DReichLab/AdmixTools tag v7.0.2) as
# released: qpDstat, qp3Pop, qpAdm, qpWave, qpGraph, qpF4ratio and
# qpfstats, whose outputs popmix's parsers read (generate_classic.py), and
# convertf, which makes the small genotype files they run on.
#
# Usage: tests/golden/admixtools_classic/build_classic.sh [install_dir]
#   default install_dir: ~/.local/share/popmix-references/admixtools-7.0.2/bin
#
# Needs git, gcc, GSL and OpenBLAS.

set -euo pipefail

dest="${1:-$HOME/.local/share/popmix-references/admixtools-7.0.2/bin}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone --quiet --depth 1 --branch v7.0.2 https://github.com/DReichLab/AdmixTools.git "$work/src"
openblas_inc="$(ls -d /usr/include/x86_64-linux-gnu/openblas* /usr/include/openblas 2>/dev/null | head -1 || true)"
progs="qpDstat qp3Pop qpAdm qpWave qpGraph qpF4ratio qpfstats convertf"
# conda environments set compiler flags for their own toolchain
env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
  make --quiet -C "$work/src/src" $progs CC=/usr/bin/gcc CFLAGS="${openblas_inc:+-I$openblas_inc}" > /dev/null
mkdir -p "$dest"
for p in $progs; do cp "$work/src/src/$p" "$dest/"; done
echo "installed $progs in $dest"
