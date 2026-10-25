#!/usr/bin/env bash
#
# Build the DATES v753 that generate_dates_reference.sh needs: upstream
# priyamoorjani/DATES with dates_numbins_fix.patch applied, its helper
# scripts and dowtjack, in one directory (DATES_BIN).
#
# Usage: tests/golden/dates/build_dates.sh [install_dir]
#   default install_dir: ~/.cache/popmix-references/dates-v753/bin
#
# Needs curl, gcc, GSL and OpenBLAS.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="${1:-$HOME/.cache/popmix-references/dates-v753/bin}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

curl -sSL https://github.com/priyamoorjani/DATES/archive/refs/heads/master.tar.gz | tar -xz -C "$work"
src="$work/DATES-master/src"
patch --quiet -d "$src" -p1 < "$here/dates_numbins_fix.patch"
# the repository ships objects and binaries: build everything from source
find "$src" \( -name '*.o' -o -name '*.a' \) -delete
rm -f "$src"/{dates,dates_expfit,dowtjack,grabpars,doexpfit,expfit,expfit2,simpjack2,nicklib}
openblas_inc="$(ls -d /usr/include/x86_64-linux-gnu/openblas*/ /usr/include/openblas/ 2>/dev/null | head -1 || true)"
# conda environments set compiler flags for their own toolchain
env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
  make --quiet -C "$src/nicksrc" CC=/usr/bin/gcc > /dev/null 2>&1
env -u CFLAGS -u LDFLAGS -u CPPFLAGS -u LDLIBS \
  make --quiet -C "$src" dates dates_expfit dowtjack grabpars CC=/usr/bin/gcc \
  ${openblas_inc:+CFLAGS="-I$openblas_inc"} > /dev/null
mkdir -p "$dest"
cp "$src"/{dates,dates_expfit,dowtjack,grabpars} "$src"/perlsrc/{dates_jackknife,dates_wtjack,dates_plot,calc_error} "$dest/"
echo "installed patched DATES v753 in $dest"
