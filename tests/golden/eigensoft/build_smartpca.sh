#!/usr/bin/env bash
#
# Build the reference smartpca: EIGENSOFT 7.2.1 (the version Debian/Ubuntu
# package as eigensoft 7.2.1+dfsg) with smartpca-hiprec.patch, which prints
# eigenvalues, eigenvector entries, Tracy-Widom statistics and outlier scores
# with %.17g instead of 3-6 decimals. Nothing else changes; generate_pca.py
# checks that the stock binary (/usr/lib/eigensoft/smartpca) agrees with it to
# its printed precision.
#
# Usage: tests/golden/eigensoft/build_smartpca.sh [install_dir]
#   default install_dir: ~/.local/share/popmix-references/eigensoft-7.2.1-hiprec
#
# Needs git, gcc, GSL, OpenBLAS and the LAPACKE runtime library. Only
# fastmode (not used by the references) calls LAPACKE; lapacke-shim.h
# declares the four routines it uses, for systems without the LAPACKE header.

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
dest="${1:-$HOME/.local/share/popmix-references/eigensoft-7.2.1-hiprec}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone --quiet --depth 1 --branch v7.2.1 https://github.com/DReichLab/EIG.git "$work/EIG"
cd "$work/EIG"
git apply "$here/smartpca-hiprec.patch"
mkdir -p "$work/shim"
cp "$here/lapacke-shim.h" "$work/shim/lapacke.h"

lapacke="$(ls /usr/lib/x86_64-linux-gnu/liblapacke.so* 2>/dev/null | head -1 || true)"
lapacke="${lapacke:--llapacke}"
openblas_inc="$(ls -d /usr/include/x86_64-linux-gnu/openblas* /usr/include/openblas 2>/dev/null | head -1 || true)"

# Conda environments set CFLAGS/LDFLAGS for their own toolchain; use the
# system compiler and libraries.
cd src
env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
  make eigensrc/smartpca CC=/usr/bin/gcc LDFLAGS= \
  CFLAGS="-O2 -I$work/shim -I$work/EIG/include ${openblas_inc:+-I$openblas_inc}" \
  LDLIBS="-lgsl -lopenblas -lm -lpthread $lapacke"

mkdir -p "$dest"
cp eigensrc/smartpca "$dest/smartpca"
echo "installed $dest/smartpca"
