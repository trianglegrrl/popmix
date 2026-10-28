#!/usr/bin/env bash
# Install the wheel just built in dist/ into a fresh virtual environment of
# the Python on PATH, then run it (check-installed.sh). release.yml runs this
# on the oldest and newest Python.
set -euo pipefail

python -m venv check-env
if [ -d check-env/Scripts ]; then bin=check-env/Scripts; else bin=check-env/bin; fi
# the built wheel, then the dependencies it declares, from PyPI
"$bin/python" -m pip install popmix --no-index --find-links dist --no-deps
"$bin/python" -m pip install popmix --find-links dist
PATH="$PWD/$bin:$PATH" bash "$(dirname "$0")/check-installed.sh"
rm -rf check-env
