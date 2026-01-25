#!/bin/bash
set -ex

# Build with maturin
${PYTHON} -m pip install . -vv --no-deps --no-build-isolation
