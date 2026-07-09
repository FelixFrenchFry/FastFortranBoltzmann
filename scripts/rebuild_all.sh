#!/usr/bin/env bash
set -euo pipefail

presets=(
  debug-fp32
  debug-fp64
  release-fp32
  release-fp64
  release-fp32-native
  release-fp64-native
)

targets=(
  FFB
  FFB_DXY
)

rm -rf build

for preset in "${presets[@]}"; do
  echo
  echo "----- [ configuring preset: ${preset} ] -----"
  cmake --preset "${preset}"

  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}" --target "${targets[@]}"
done
