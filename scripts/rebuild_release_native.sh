#!/usr/bin/env bash
set -euo pipefail

presets=(
  release-fp32-native
  release-fp64-native
)

targets=(
  FFB
  FFB_CHECKPOINT_0
  FFB_CHECKPOINT_1
  FFB_CHECKPOINT_2
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
