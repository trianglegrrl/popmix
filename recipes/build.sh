#!/bin/bash
set -ex

# THIRDPARTY.yml, the licences of the Rust crates compiled in, ships in the
# sdist (.github/scripts/thirdparty-licenses.py)
${PYTHON} -m pip install . -vv --no-deps --no-build-isolation
