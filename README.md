![Fortran 2018](https://img.shields.io/badge/Fortran-2018-734f96) ![Intel ifx](https://img.shields.io/badge/compiler-Intel%20ifx-0071c5)

# FastFortranBoltzmann 🧮

Distributed D2Q9 LBM fluid solver using Fortran coarrays for inter-node communication on CPU clusters.

The solver supports shear-wave decay, Couette flow, Poiseuille flow, and lid-driven cavity fluid simulation scenarios. The project investigates kernel fusion, memory layouts, domain decomposition, SIMD vectorization, FP precision, Fortran coarrays, MPI pinning and other aspects using a Fortran 2018 implementation. Comprehensive performance benchmarks for solver runs on up to 960 CPU cores across 80 cluster nodes in parallel are provided.

## How to build

Requires Linux, Intel's `ifx` Fortran compiler, Intel MPI, CMake 3.25 or newer, and Ninja. For example, build the `FFB` solver in release mode with double precision and optimizations for your CPU architecture with:

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

For example, run one simulation on 8 coarray images in a 2 × 4 decomposition with MPI processes pinned to cores. This prints the hardware and simulation setup, live progress, and a final timing and throughput summary, and writes run metadata to `output/run_000/config.json`:

```bash
I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core I_MPI_PIN_ORDER=spread I_X=2 I_Y=4 FOR_COARRAY_NUM_IMAGES=8 ./build/release-fp64-native/bin/FFB
```

| Variable | Meaning | Range | Required |
|---|---|---|---|
| `FOR_COARRAY_NUM_IMAGES` | Total number of coarray images | Positive integer equal to `I_X × I_Y` | Yes |
| `I_X`, `I_Y` | Number of subdomains along the x-axis and y-axis | Positive divisors of `N_X` and `N_Y` | Yes |
| `I_MPI_PIN` | Enables pinning each image to fixed processor cores | `1` (default), `0` | No |
| `I_MPI_PIN_DOMAIN` | Hardware unit whose cores each image is pinned to | `core`, `socket`, `numa`, `node`, `cache1`, or [other forms](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-18/interoperability-with-openmp-api.html) | No |
| `I_MPI_PIN_ORDER` | Order in which images are assigned to pinning domains | `compact` (default), `scatter`, `spread`, `bunch` | No |

## How to benchmark

For example, benchmark one configuration with [benchmark_run.py](scripts/benchmark_run.py), or compare several configurations in a predefined list with [benchmark_sweep.py](scripts/benchmark_sweep.py). Both print the median, best, worst, mean, stddev throughput in MLUPS, and a timing breakdown for each configuration:

```bash
# benchmark one configuration, computing statistics over 5 runs
python3 scripts/benchmark_run.py --exe build/release-fp64-native/bin/FFB --runs 5 --images 8 --ix 2 --iy 4 --pin core_spread

# benchmark several configurations in a "custom" config file, computing statistics over 5 runs each
python3 scripts/benchmark_sweep.py --exe build/release-fp64-native/bin/FFB --runs 5 --case custom --skip 0 --pin core_spread
```

| Option | Meaning | Range | Required |
|---|---|---|---|
| `--images`, `--ix`, `--iy` | Image count and subdomains along the x-axis and y-axis for `benchmark_run.py` | Same as `FOR_COARRAY_NUM_IMAGES`, `I_X`, `I_Y` | Yes, for `benchmark_run.py` |
| `--exe` | Path to the solver executable | `build/release-fp64-native/bin/FFB` (default) or any other built solver | No |
| `--runs` | Run repetitions per configuration | Positive integer, `5` (default) | No |
| `--pin` | Pinning preset that sets the `I_MPI_PIN*` variables | `none`, `core_scatter`, `core_spread` (default), `numa_scatter`, `numa_spread` | No |
| `--case` | Sweep of image grids for `benchmark_sweep.py` | `custom` (default) or any sweep listed by `--help` | No |
| `--skip` | Number of image grids `benchmark_sweep.py` skips at the start | Non-negative integer below the sweep length, `0` (default) | No |
| `--timeout` | Time limit per solver run in seconds | Positive integer, `900` (default) | No |
| `--per-host` | Generates a balanced hostfile inside a Slurm allocation | Flag, off (default) | No |
