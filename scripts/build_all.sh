#!/usr/bin/env bash
set -euo pipefail

presets=(
  debug
  debug-fp64
  release
  release-fp64
  release-native
  release-fp64-native
)

for preset in "${presets[@]}"; do
  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}"
done
