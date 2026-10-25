#!/usr/bin/env bash
#
# Build convertf and mergeit from AdmixTools 8.0.2 (DReichLab/AdmixTools tag v8.0.2),
# the first release line that reads and writes TGENO (transpose_packed),
# the format of AADR v66.p1. generate_formats.py and generate_merge.sh use them.
#
# Usage: tests/golden/formats/build_convertf.sh [install_dir]
#   default install_dir: ~/.local/share/popmix-references/admixtools-8.0.2/bin
#
# Needs git, gcc, GSL and OpenBLAS.

set -euo pipefail

dest="${1:-$HOME/.local/share/popmix-references/admixtools-8.0.2/bin}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone --quiet --depth 1 --branch v8.0.2 https://github.com/DReichLab/AdmixTools.git "$work/src"
openblas_inc="$(ls -d /usr/include/x86_64-linux-gnu/openblas*/ /usr/include/openblas/ 2>/dev/null | head -1 || true)"
# conda environments set compiler flags for their own toolchain
env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
  make --quiet -C "$work/src/src" convertf mergeit CC=/usr/bin/gcc CFLAGS="${openblas_inc:+-I$openblas_inc}" > /dev/null
mkdir -p "$dest"
cp "$work/src/src/convertf" "$work/src/src/mergeit" "$dest/"
echo "installed convertf and mergeit in $dest"
