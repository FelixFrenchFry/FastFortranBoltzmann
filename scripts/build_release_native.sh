#!/usr/bin/env bash
set -euo pipefail

presets=(
  release-native
)

for preset in "${presets[@]}"; do
  echo
  echo "----- [ building preset: ${preset} ] -----"
  cmake --build --preset "${preset}"
done
