#!/bin/bash
set -ex

# licenses of the Rust crates bundled into the extension
cargo-bundle-licenses --format yaml --output THIRDPARTY.yml

${PYTHON} -m pip install . -vv --no-deps --no-build-isolation
