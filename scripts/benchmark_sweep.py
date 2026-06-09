#!/usr/bin/env python3
"""Run a sim config multiple times and summarize performance metrics"""

import argparse
import os
import re
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path



NUMBER = r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
STEP_RE = re.compile(rf"step time:\s*({NUMBER})\s*ms", re.IGNORECASE)
MLUPS_RE = re.compile(rf"MLUPS:\s*({NUMBER})", re.IGNORECASE)
EXECUTION_TIME_RE = re.compile(rf"^([a-z ]+?)\s*\|\s*({NUMBER})\s*\|\s*({NUMBER})\s*%", re.IGNORECASE | re.MULTILINE)
LAUNCHED_RE = re.compile(r"^\[[0-9:]+\]\s+launched", re.MULTILINE)
SIM_SIZE_RE = re.compile(r"^\s*integer\(int32\), parameter :: (N_[XY])\s*=\s*([0-9]+)", re.MULTILINE)
CMAKE_CACHE_RE = re.compile(r"^([^#/\n][^:=\n]*):[^=\n]*=(.*)$", re.MULTILINE)
HEADER_WIDTH = 75
MAX_RUNS = 999
TIMING_CATEGORIES = (
    "kernel compute",
    "halo exchange",
    "other",
    "total",
)
PINNING_ENV_NAMES = (
    "I_MPI_PIN",
    "I_MPI_PIN_DOMAIN",
    "I_MPI_PIN_ORDER",
    "I_MPI_PIN_PROCESSOR_LIST",
)
PINNING_PRESETS = {
    "none": {
        "I_MPI_PIN": "0",
    },
    "core_scatter": {
        "I_MPI_PIN": "1",
        "I_MPI_PIN_DOMAIN": "core",
        "I_MPI_PIN_ORDER": "scatter",
    },
    "core_spread": {
        "I_MPI_PIN": "1",
        "I_MPI_PIN_DOMAIN": "core",
        "I_MPI_PIN_ORDER": "spread",
    },
    "numa_scatter": {
        "I_MPI_PIN": "1",
        "I_MPI_PIN_DOMAIN": "numa",
        "I_MPI_PIN_ORDER": "scatter",
    },
    "numa_spread": {
        "I_MPI_PIN": "1",
        "I_MPI_PIN_DOMAIN": "numa",
        "I_MPI_PIN_ORDER": "spread",
    },
}

# settings
DEFAULT_EXE = "build/release-fp32/bin/FFB"
DEFAULT_RUNS = 5
DEFAULT_PIN = "numa_spread"
SKIP = 0

# domain decompositions (images, ix, iy)
CASE_CUSTOM = [
    (4, 1, 4),
    (6, 1, 6),
    (9, 1, 9),
    (12, 1, 12),
    (16, 1, 16),
    (20, 1, 20),
    (25, 1, 25),
    (30, 1, 30),
    (36, 1, 36),
    (40, 1, 40),
    (48, 1, 48),
    (64, 1, 64),
    (72, 1, 72),
    (80, 1, 80),
    (96, 1, 96),
    (120, 1, 120),
    (144, 1, 144),
    (160, 1, 160),
    (192, 1, 192)
]

CASE_SQUARES = [
    # 1D horizontal slice (ix=1)
    (4, 1, 4),
    (9, 1, 9),
    (16, 1, 16),
    (25, 1, 25),
    (36, 1, 36),
    (64, 1, 64),
    (100, 1, 100),
    (144, 1, 144),
    (225, 1, 225)
    # 2D square (ix=iy)
    (4, 2, 2),
    (9, 3, 3),
    (16, 4, 4),
    (25, 5, 5),
    (36, 6, 6),
    (64, 8, 8),
    (100, 10, 10),
    (144, 12, 12),
    (225, 15, 15)
    # 1D vertical slice (iy=1)
    (4, 4, 1),
    (9, 9, 1),
    (16, 16, 1),
    (25, 25, 1),
    (36, 36, 1),
    (64, 64, 1),
    (100, 100, 1),
    (144, 144, 1),
    (225, 225, 1)
]

CASE_SMALL = [
    (4, 1, 4),
    (6, 1, 6),
    (9, 1, 9),
    (12, 1, 12),
    (16, 1, 16),
    (20, 1, 20),
    (25, 1, 25),
    (36, 1, 36),
    (48, 1, 48),
    (64, 1, 64)
]

CASE_SMALL_FULL = [
    (1, 1, 1),
    (2, 1, 2),
    (3, 1, 3),
    (4, 1, 4),
    (5, 1, 5),
    (6, 1, 6),
    (8, 1, 8),
    (9, 1, 9),
    (10, 1, 10),
    (12, 1, 12),
    (15, 1, 15),
    (16, 1, 16),
    (18, 1, 18),
    (20, 1, 20),
    (24, 1, 24),
    (25, 1, 25),
    (30, 1, 30),
    (32, 1, 32),
    (36, 1, 36),
    (40, 1, 40),
    (45, 1, 45),
    (48, 1, 48),
    (50, 1, 50),
    (60, 1, 60),
    (64, 1, 64)
]

CASE_MEDIUM = [
    (64, 1, 64),
    (80, 1, 80),
    (96, 1, 96),
    (100, 1, 100),
    (120, 1, 120),
    (144, 1, 144),
    (160, 1, 160),
    (192, 1, 192),
    (225, 1, 225),
    (240, 1, 240),
    (288, 1, 288),
    (320, 1, 320),
    (360, 1, 360),
    (400, 1, 400),
    (480, 1, 480),
    (576, 1, 576),
    (600, 1, 600)
]

CASE_MEDIUM_FULL = [
    (64, 1, 64),
    (72, 1, 72),
    (75, 1, 75),
    (80, 1, 80),
    (90, 1, 90),
    (96, 1, 96),
    (100, 1, 100),
    (120, 1, 120),
    (144, 1, 144),
    (150, 1, 150),
    (160, 1, 160),
    (180, 1, 180),
    (192, 1, 192),
    (200, 1, 200),
    (225, 1, 225),
    (240, 1, 240),
    (288, 1, 288),
    (300, 1, 300),
    (320, 1, 320),
    (360, 1, 360),
    (400, 1, 400),
    (450, 1, 450),
    (480, 1, 480),
    (576, 1, 576),
    (600, 1, 600)
]

CASE_LARGE = [
    (600, 1, 600),
    (720, 1, 720),
    (800, 1, 800),
    (900, 1, 900),
    (960, 1, 960),
    (1200, 1, 1200),
    (1440, 1, 1440),
    (1600, 1, 1600),
    (1800, 1, 1800),
    (2400, 1, 2400),
    (2880, 1, 2880),
    (3600, 1, 3600),
    (4800, 1, 4800),
    (7200, 1, 7200),
    (14400, 1, 14400)
]

CASE_LARGE_FULL = [
    (600, 1, 600),
    (720, 1, 720),
    (800, 1, 800),
    (900, 1, 900),
    (960, 1, 960),
    (1200, 1, 1200),
    (1440, 1, 1440),
    (1600, 1, 1600),
    (1800, 1, 1800),
    (2400, 1, 2400),
    (2880, 1, 2880),
    (3600, 1, 3600),
    (4800, 1, 4800),
    (7200, 1, 7200),
    (14400, 1, 14400)
]

CASE_ALL = [
    (1, 1, 1),
    (2, 1, 2),
    (3, 1, 3),
    (4, 1, 4),
    (5, 1, 5),
    (6, 1, 6),
    (8, 1, 8),
    (9, 1, 9),
    (10, 1, 10),
    (12, 1, 12),
    (15, 1, 15),
    (16, 1, 16),
    (18, 1, 18),
    (20, 1, 20),
    (24, 1, 24),
    (25, 1, 25),
    (30, 1, 30),
    (32, 1, 32),
    (36, 1, 36),
    (40, 1, 40),
    (45, 1, 45),
    (48, 1, 48),
    (50, 1, 50),
    (60, 1, 60),
    (64, 1, 64),
    (72, 1, 72),
    (75, 1, 75),
    (80, 1, 80),
    (90, 1, 90),
    (96, 1, 96),
    (100, 1, 100),
    (120, 1, 120),
    (144, 1, 144),
    (150, 1, 150),
    (160, 1, 160),
    (180, 1, 180),
    (192, 1, 192),
    (200, 1, 200),
    (225, 1, 225),
    (240, 1, 240),
    (288, 1, 288),
    (300, 1, 300),
    (320, 1, 320),
    (360, 1, 360),
    (400, 1, 400),
    (450, 1, 450),
    (480, 1, 480),
    (576, 1, 576),
    (600, 1, 600),
    (720, 1, 720),
    (800, 1, 800),
    (900, 1, 900),
    (960, 1, 960),
    (1200, 1, 1200),
    (1440, 1, 1440),
    (1600, 1, 1600),
    (1800, 1, 1800),
    (2400, 1, 2400),
    (2880, 1, 2880),
    (3600, 1, 3600),
    (4800, 1, 4800),
    (7200, 1, 7200),
    (14400, 1, 14400)
]

DOMAIN_DECOMP_CASE_SETS = {
    "custom": CASE_CUSTOM,
    "squares": CASE_SQUARES,
    "small": CASE_SMALL,
    "small_full": CASE_SMALL_FULL,
    "medium": CASE_MEDIUM,
    "medium_full": CASE_MEDIUM_FULL,
    "large": CASE_LARGE,
    "large_full": CASE_LARGE_FULL,
    "all": CASE_ALL,
}


def print_header(title):
    print(f"--- [ {title} ] " + "-" * (HEADER_WIDTH - len(title) - 9))


def print_param(name, value):
    print(f"{name:<25} = {value}")


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", default=DEFAULT_EXE)
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    parser.add_argument("--case", choices=DOMAIN_DECOMP_CASE_SETS.keys(), default="custom")
    parser.add_argument("--skip", type=int, default=SKIP)
    parser.add_argument("--pin", choices=PINNING_PRESETS.keys(), default=DEFAULT_PIN)
    return parser.parse_args()


def read_source_sim_size():
    settings_path = Path(__file__).resolve().parents[1] / "source" / "settings.f90"
    settings = settings_path.read_text()
    values = {name: int(value) for name, value in SIM_SIZE_RE.findall(settings)}

    if "N_X" not in values or "N_Y" not in values:
        raise RuntimeError("could not parse N_X/N_Y from source/settings.f90")

    return values["N_X"], values["N_Y"]


def read_cmake_cache(path):
    cache = path.read_text()
    return {name: value.strip() for name, value in CMAKE_CACHE_RE.findall(cache)}


def read_cmake_settings_definitions(cache_values):
    definitions = cache_values.get("FFB_SETTINGS_DEFINITIONS", "")
    if not definitions:
        return {}

    values = {}
    for definition in definitions.split(";"):
        definition = definition.strip()
        if not definition:
            continue

        if "=" not in definition:
            values[definition] = ""
            continue

        name, value = definition.split("=", 1)
        values[name] = value

    return values


def read_sim_size(exe):
    exe_path = Path(exe).resolve()
    build_dir = exe_path.parent.parent if exe_path.parent.name == "bin" else exe_path.parent
    cache_path = build_dir / "CMakeCache.txt"

    if cache_path.exists():
        cache_values = read_cmake_cache(cache_path)
        cmake_settings = read_cmake_settings_definitions(cache_values)

        if "FFB_USE_CMAKE_SETTINGS" in cmake_settings:
            if "FFB_N_X" not in cmake_settings or "FFB_N_Y" not in cmake_settings:
                raise RuntimeError("could not parse FFB_N_X/FFB_N_Y from FFB_SETTINGS_DEFINITIONS")
            return int(cmake_settings["FFB_N_X"]), int(cmake_settings["FFB_N_Y"])

    return read_source_sim_size()


def parse_last(pattern, text, name):
    matches = pattern.findall(text)
    if not matches:
        raise RuntimeError(f"could not parse {name}")
    return float(matches[-1])


def parse_execution_times(output):
    values = {}
    for name, seconds, share in EXECUTION_TIME_RE.findall(output):
        values[name.strip().lower()] = {
            "seconds": float(seconds),
            "share": float(share),
        }

    missing = [name for name in TIMING_CATEGORIES if name not in values]
    if missing:
        raise RuntimeError("could not parse execution time table")

    return values


def apply_pinning_preset(env, pin):
    for name in PINNING_ENV_NAMES:
        env.pop(name, None)

    env.update(PINNING_PRESETS[pin])


def print_pinning_settings(pin):
    env = " ".join(f"{name}={value}" for name, value in PINNING_PRESETS[pin].items())

    print_param("mpi pinning preset", pin)
    print_param("mpi pinning env", env)


def run_once(exe, images, ix, iy, run_num, pin):
    env = os.environ.copy()
    apply_pinning_preset(env, pin)
    if "FOR_COARRAY_CONFIG_FILE" not in env:
        env["FOR_COARRAY_NUM_IMAGES"] = str(images)
    env["I_X"] = str(ix)
    env["I_Y"] = str(iy)

    completed = subprocess.run(
        [exe],
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    output = completed.stdout + "\n" + completed.stderr
    if completed.returncode != 0:
        print(output)
        raise RuntimeError(f"run {run_num} failed with exit code {completed.returncode}")

    return (
        parse_last(STEP_RE, output, "step time"),
        parse_last(MLUPS_RE, output, "MLUPS"),
        parse_execution_times(output),
        output,
    )


def print_static_app_output(output):
    match = LAUNCHED_RE.search(output)
    if not match:
        return

    print(output[:match.start()].strip())
    print()


def get_stats(values, higher_is_better):
    best = max(values) if higher_is_better else min(values)
    worst = min(values) if higher_is_better else max(values)
    mean = statistics.mean(values)
    stddev = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "median": statistics.median(values),
        "best": best,
        "worst": worst,
        "mean": mean,
        "stddev_percent": 100.0 * stddev / mean if mean != 0.0 else 0.0,
    }


def print_execution_time_stats(execution_times):
    print()
    print(f"{'execution time median':<40} | {'total [sec]':>14} | {'share [%]':>15}")
    print("-" * HEADER_WIDTH)

    for category in TIMING_CATEGORIES:
        seconds = [values[category]["seconds"] for values in execution_times]
        shares = [values[category]["share"] for values in execution_times]

        print(f"{category:<40} | {statistics.median(seconds):14.3f} | {statistics.median(shares):13.3f} %")


def validate_case(case_num, n_x, n_y, images, ix, iy):
    if images <= 0:
        sys.exit(f"error: case {case_num}: total images must be positive")
    if ix <= 0 or iy <= 0:
        sys.exit(f"error: case {case_num}: I_X and I_Y must be positive")
    if images != ix * iy:
        sys.exit(f"error: case {case_num}: total images must match I_X * I_Y")
    if n_x % ix != 0:
        sys.exit(f"error: case {case_num}: N_X must be divisible by I_X")
    if n_y % iy != 0:
        sys.exit(f"error: case {case_num}: N_Y must be divisible by I_Y")


def run_case(exe, runs, n_x, n_y, case_num, n_cases, images, ix, iy, pin):
    validate_case(case_num, n_x, n_y, images, ix, iy)

    mlups_values = []
    execution_times = []

    print()
    print_header("benchmark script settings")
    print_param("case", f"{case_num} / {n_cases}")
    print_param("executable", exe)
    print_param("runs", runs)
    print_param("images", images)
    print_param("image grid", f"{ix} x {iy}")
    print_param("sim size", f"{n_x} x {n_y}")
    print_pinning_settings(pin)
    print()

    runs_started_at = timestamp()

    for run_num in range(1, runs + 1):
        step_ms, mlups, execution_time, output = run_once(exe, images, ix, iy, run_num, pin)

        if run_num == 1:
            print_static_app_output(output)
            print_header(f"benchmark runs started at {runs_started_at}")

        print(f"{run_num:03d} | avg step: {step_ms:.3f} ms | MLUPS: {mlups:.3f}")

        mlups_values.append(mlups)
        execution_times.append(execution_time)

    print()
    mlups_stats = get_stats(mlups_values, higher_is_better=True)

    print_header("MLUPS metrics")
    print_param("median", f"{mlups_stats['median']:.3f}")
    print_param("best", f"{mlups_stats['best']:.3f}")
    print_param("worst", f"{mlups_stats['worst']:.3f}")
    print_param("mean", f"{mlups_stats['mean']:.3f}")
    print_param("stddev", f"{mlups_stats['stddev_percent']:.3f} %")

    print_execution_time_stats(execution_times)


def main():
    args = parse_args()

    if args.runs <= 0:
        sys.exit("error: --runs must be positive")
    if args.runs > MAX_RUNS:
        sys.exit(f"error: --runs must be <= {MAX_RUNS}")
    if args.skip < 0:
        sys.exit("error: --skip must not be negative")
    domain_decomp_cases = DOMAIN_DECOMP_CASE_SETS[args.case]
    if not domain_decomp_cases:
        sys.exit(f"error: case list '{args.case}' must not be empty")
    if args.skip >= len(domain_decomp_cases):
        sys.exit("error: --skip must be smaller than the selected case list length")
    if not os.path.exists(args.exe):
        sys.exit(f"error: executable not found: {args.exe}")

    n_x, n_y = read_sim_size(args.exe)

    for case_num, (images, ix, iy) in enumerate(domain_decomp_cases[args.skip:], start=args.skip + 1):
        run_case(
            args.exe, args.runs, n_x, n_y, case_num, len(domain_decomp_cases), images, ix, iy, args.pin)


if __name__ == "__main__":
    main()
