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

for preset in "${presets[@]}"; do
  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}"
done
