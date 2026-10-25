#!/usr/bin/env bash
#
# Build the reference ALDER: version 1.03 from its authors' site, built
# twice, as released ("stock") and with alder-hiprec.patch, which prints the
# weighted LD curve, the fits, the LD-correlation table and the admixture
# test with %.17g instead of 2-8 decimals. Nothing else changes;
# generate_alder.py checks that the two agree to the stock build's printed
# precision.
#
# ALDER is licensed for academic and non-profit use only; it is downloaded
# and built here as a reference and is not distributed with popmix.
#
# Usage: tests/golden/alder/build_alder.sh [install_dir]
#   default install_dir: ~/.local/share/popmix-references/alder-1.03
#
# Needs curl, g++ with OpenMP, FFTW 3 and LAPACK.

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
dest="${1:-$HOME/.local/share/popmix-references/alder-1.03}"
url="https://cb.csail.mit.edu/alder/alder_v1.03.tar.gz"
sha256="aae68848dc444946aefe171322e0d98b1a0cc0c82f3669007ad9608d331b7241"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

curl -sSL -o "$work/alder.tar.gz" "$url"
echo "$sha256  $work/alder.tar.gz" | sha256sum -c -

build() {
  local dir="$1"
  mkdir -p "$dir"
  tar -xzf "$work/alder.tar.gz" -C "$dir"
  rm -f "$dir/alder/alder"  # the archive ships a prebuilt binary
  if [ "${2:-}" = patch ]; then
    (cd "$dir/alder" && patch -p1 --quiet < "$here/alder-hiprec.patch")
  fi
  # conda environments set compiler flags for their own toolchain
  env -u CFLAGS -u LDFLAGS -u CXXFLAGS -u CPPFLAGS \
    make --quiet -C "$dir/alder" CXX=/usr/bin/g++ CC=/usr/bin/gcc > /dev/null
}

build "$work/stock"
build "$work/hiprec" patch
mkdir -p "$dest"
cp "$work/stock/alder/alder" "$dest/alder-stock"
cp "$work/hiprec/alder/alder" "$dest/alder-hiprec"
echo "installed $dest/alder-stock and $dest/alder-hiprec"
