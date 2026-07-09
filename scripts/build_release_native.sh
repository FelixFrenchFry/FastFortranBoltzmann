#!/usr/bin/env bash
set -euo pipefail

presets=(
  release-fp32-native
  release-fp64-native
)

targets=(
  FFB
  FFB_DXY
)

for preset in "${presets[@]}"; do
  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}" --target "${targets[@]}"
done
