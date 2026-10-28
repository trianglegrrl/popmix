#!/usr/bin/env bash
# Run the popmix installed for the Python on PATH: the CLI, the bundled agent
# skill and the compiled module. With EXPECTED_VERSION set, check that the
# installed distribution carries that version. check-wheel.sh runs this on a
# wheel just built; release.yml on one installed back from TestPyPI.
set -euo pipefail

python --version
popmix --version
rm -rf skill-check
popmix install-skill --dest skill-check
test -f skill-check/commands/qpadm.md
python -c "import popmix, popmix._core; print(popmix.__version__)"
if [ -n "${EXPECTED_VERSION:-}" ]; then
  python -c "import importlib.metadata as m, os, sys; v = m.version('popmix'); sys.exit(None if v == os.environ['EXPECTED_VERSION'] else f'installed {v}, expected ' + os.environ['EXPECTED_VERSION'])"
fi
rm -rf skill-check
