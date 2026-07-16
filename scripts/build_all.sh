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
  FFB_CHECKPOINT_0
  FFB_CHECKPOINT_1
  FFB_CHECKPOINT_2
)

for preset in "${presets[@]}"; do
  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}" --target "${targets[@]}"
done
