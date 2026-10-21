![Fortran 2018](https://img.shields.io/badge/Fortran-2018-734f96) ![Intel ifx](https://img.shields.io/badge/compiler-Intel%20ifx-0071c5)

# FastFortranBoltzmann 🧮

Distributed D2Q9 LBM fluid solver using Fortran coarrays for inter-node communication on CPU clusters.

The solver supports shear-wave decay, Couette flow, Poiseuille flow, and lid-driven cavity fluid simulation scenarios. The project investigates kernel fusion, memory layouts, domain decomposition, SIMD vectorization, FP precision, Fortran coarrays, MPI pinning and other aspects using a Fortran 2018 implementation. Comprehensive performance benchmarks for solver runs on up to 960 CPU cores across 80 cluster nodes in parallel are provided.

## How to build

Requires Linux, Intel's `ifx` Fortran compiler, Intel MPI, CMake 3.25 or newer, and Ninja. From the repository root, build the `FFB` solver in release mode with double precision and optimizations for your CPU architecture with:

```bash
# e.g. activate a standard Intel oneAPI environment
source /opt/intel/oneapi/setvars.sh

cmake --preset release-fp64-native
cmake --build --preset release-fp64-native --target FFB
```

Choose a build preset from [CMakePresets.json](CMakePresets.json) using the table below, and use the same preset in both CMake commands:

| Build | Single precision | Double precision |
|---|---|---|
| Debug (`-O0`) | `debug-fp32` | `debug-fp64` |
| Release (`-O3`) | `release-fp32` | `release-fp64` |
| Release, native CPU (`-O3 -march=native`) | `release-fp32-native` | `release-fp64-native` |

## How to run

TODO
