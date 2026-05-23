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
}

# settings
DEFAULT_EXE = "build/release/bin/FFB"
DEFAULT_RUNS = 5
DEFAULT_PIN = "core_scatter"


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
    parser.add_argument("--images", type=int, required=True)
    parser.add_argument("--ix", type=int, required=True)
    parser.add_argument("--iy", type=int, required=True)
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


def main():
    args = parse_args()

    if args.runs <= 0:
        sys.exit("error: --runs must be positive")
    if args.runs > MAX_RUNS:
        sys.exit(f"error: --runs must be <= {MAX_RUNS}")
    if args.images <= 0:
        sys.exit("error: --images must be positive")
    if args.ix <= 0 or args.iy <= 0:
        sys.exit("error: --ix and --iy must be positive")
    if args.images != args.ix * args.iy:
        sys.exit("error: --images must match --ix * --iy")
    if not os.path.exists(args.exe):
        sys.exit(f"error: executable not found: {args.exe}")
    n_x, n_y = read_sim_size(args.exe)
    if n_x % args.ix != 0:
        sys.exit("error: N_X must be divisible by --ix")
    if n_y % args.iy != 0:
        sys.exit("error: N_Y must be divisible by --iy")

    step_times = []
    mlups_values = []
    execution_times = []

    print()
    print_header("benchmark script settings")
    print_param("executable", args.exe)
    print_param("runs", args.runs)
    print_param("images", args.images)
    print_param("image grid", f"{args.ix} x {args.iy}")
    print_param("sim size", f"{n_x} x {n_y}")
    print_pinning_settings(args.pin)
    print()

    runs_started_at = timestamp()

    for run_num in range(1, args.runs + 1):
        step_ms, mlups, execution_time, output = run_once(
            args.exe, args.images, args.ix, args.iy, run_num, args.pin)

        if run_num == 1:
            print_static_app_output(output)
            print_header(f"benchmark runs started at {runs_started_at}")

        print(f"{run_num:03d} | avg step: {step_ms:.3f} ms | MLUPS: {mlups:.3f}")

        step_times.append(step_ms)
        mlups_values.append(mlups)
        execution_times.append(execution_time)

    print()
    step_stats = get_stats(step_times, higher_is_better=False)
    mlups_stats = get_stats(mlups_values, higher_is_better=True)

    print_header("avg step metrics")
    print_param("median", f"{step_stats['median']:.3f} ms")
    print_param("best", f"{step_stats['best']:.3f} ms")
    print_param("worst", f"{step_stats['worst']:.3f} ms")
    print_param("mean", f"{step_stats['mean']:.3f} ms")
    print_param("stddev", f"{step_stats['stddev_percent']:.3f} %")

    print()
    print_header("MLUPS metrics")
    print_param("median", f"{mlups_stats['median']:.3f}")
    print_param("best", f"{mlups_stats['best']:.3f}")
    print_param("worst", f"{mlups_stats['worst']:.3f}")
    print_param("mean", f"{mlups_stats['mean']:.3f}")
    print_param("stddev", f"{mlups_stats['stddev_percent']:.3f} %")

    print_execution_time_stats(execution_times)


if __name__ == "__main__":
    main()
