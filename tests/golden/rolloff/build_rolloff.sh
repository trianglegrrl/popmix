#!/usr/bin/env bash
#
# Build the reference ROLLOFF: rolloff and dowtjack from ADMIXTOOLS 7.0.2
# (DReichLab/AdmixTools tag v7.0.2), as released ("stock") and with
# rolloff-hiprec.patch, which prints rolloff's output file and dowtjack's
# estimate with %.17g instead of 3-6 decimals. Nothing else changes;
# generate_rolloff.py checks the two agree to the stock build's precision.
#
# Usage: tests/golden/rolloff/build_rolloff.sh [install_dir]
#   default install_dir: ~/.local/share/popmix-references/admixtools-7.0.2
#
# Needs git, gcc, GSL and OpenBLAS. The fits use ADMIXTOOLS' rexpfit.r,
# which needs R with the DEoptim package.

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
dest="${1:-$HOME/.local/share/popmix-references/admixtools-7.0.2}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone --quiet --depth 1 --branch v7.0.2 https://github.com/DReichLab/AdmixTools.git "$work/stock"
cp -r "$work/stock" "$work/hiprec"
(cd "$work/hiprec" && git apply "$here/rolloff-hiprec.patch")

openblas_inc="$(ls -d /usr/include/x86_64-linux-gnu/openblas* /usr/include/openblas 2>/dev/null | head -1 || true)"
for b in stock hiprec; do
  # conda environments set compiler flags for their own toolchain
  env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
    make --quiet -C "$work/$b/src" rolloff dowtjack CC=/usr/bin/gcc \
    CFLAGS="${openblas_inc:+-I$openblas_inc}" > /dev/null
done
mkdir -p "$dest"
cp "$work/stock/src/rolloff" "$dest/rolloff-stock"
cp "$work/hiprec/src/rolloff" "$dest/rolloff-hiprec"
cp "$work/hiprec/src/dowtjack" "$dest/dowtjack-hiprec"
cp "$work/stock/src/rexpfit.r" "$dest/rexpfit.r"
echo "installed rolloff-stock, rolloff-hiprec, dowtjack-hiprec and rexpfit.r in $dest"
