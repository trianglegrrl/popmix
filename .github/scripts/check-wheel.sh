#!/usr/bin/env bash
# Install the wheel just built in dist/ into a fresh virtual environment of
# the Python on PATH, then run it (check-installed.sh). release.yml runs this
# on the oldest and newest Python.
set -euo pipefail

# the licences every wheel must carry (NOTICE, THIRDPARTY.yml: TASKS.md C8)
python - <<'PY'
import glob, sys, zipfile
for wheel in glob.glob("dist/*.whl"):
    names = zipfile.ZipFile(wheel).namelist()
    for required in ("licenses/LICENSE", "licenses/NOTICE", "licenses/THIRDPARTY.yml"):
        if not any(n.endswith(".dist-info/" + required) for n in names):
            sys.exit(f"{wheel}: no {required}")
PY

python -m venv check-env
if [ -d check-env/Scripts ]; then bin=check-env/Scripts; else bin=check-env/bin; fi
# the built wheel, then the dependencies it declares, from PyPI
"$bin/python" -m pip install popmix --no-index --find-links dist --no-deps
"$bin/python" -m pip install popmix --find-links dist
PATH="$PWD/$bin:$PATH" bash "$(dirname "$0")/check-installed.sh"
rm -rf check-env
