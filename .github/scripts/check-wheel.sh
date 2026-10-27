#!/usr/bin/env bash
# Install the wheel just built in dist/ into a fresh virtual environment of
# the Python on PATH, then run it: the CLI, the bundled agent skill and the
# compiled module. release.yml runs this on the oldest and newest Python.
set -euo pipefail

python -m venv check-env
if [ -d check-env/Scripts ]; then bin=check-env/Scripts; else bin=check-env/bin; fi
# the built wheel, then the dependencies it declares, from PyPI
"$bin/python" -m pip install popmix --no-index --find-links dist --no-deps
"$bin/python" -m pip install popmix --find-links dist
export PATH="$PWD/$bin:$PATH"
python --version
popmix --version
rm -rf skill-check
popmix install-skill --dest skill-check
test -f skill-check/commands/qpadm.md
python -c "import popmix, popmix._core; print(popmix.__version__)"
rm -rf check-env skill-check
