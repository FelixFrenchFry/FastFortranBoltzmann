#!/usr/bin/env bash
set -euo pipefail

presets=(
  release-fp32-native
  release-fp64-native
)

rm -rf build

for preset in "${presets[@]}"; do
  echo
  echo "----- [ configuring preset: ${preset} ] -----"
  cmake --preset "${preset}"

  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}"
done
