#!/usr/bin/env bash
set -euo pipefail

presets=(
  debug
  debug-fp64
  release
  release-fp64
  release-intel
  release-amd
  release-fp64-intel
  release-fp64-amd
)

rm -rf build

for preset in "${presets[@]}"; do
  echo
  echo "=== configuring: ${preset} ==="
  cmake --preset "${preset}"

  echo
  echo "=== building: ${preset} ==="
  cmake --build --preset "${preset}"
done
