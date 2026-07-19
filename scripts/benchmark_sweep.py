#!/usr/bin/env python3
"""Run a sim config multiple times and summarize performance metrics"""

import argparse
import os
import re
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path



NUMBER = r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
STEP_RE = re.compile(rf"step time\s*[:=]\s*({NUMBER})\s*ms", re.IGNORECASE)
MLUPS_RE = re.compile(rf"MLUPS\s*[:=]\s*({NUMBER})", re.IGNORECASE)
TIMING_SPREAD_RE = re.compile(
    rf"^(kernel compute|halo sync|halo transfer|other|total)\s*"
    rf"\|\s*({NUMBER})\s*\((\d+)\)\s*"
    rf"\|\s*({NUMBER})\s*\((\d+)\)\s*$",
    re.IGNORECASE | re.MULTILINE,
)
LAUNCHED_RE = re.compile(r"^\[[0-9:]+\]\s+launched", re.MULTILINE)
SIM_SIZE_RE = re.compile(r"^\s*integer\(int32\), parameter :: (N_[XY])\s*=\s*([0-9]+)", re.MULTILINE)
CMAKE_CACHE_RE = re.compile(r"^([^#/\n][^:=\n]*):[^=\n]*=(.*)$", re.MULTILINE)
HEADER_WIDTH = 80
MAX_RUNS = 999
TIMING_CATEGORIES = (
    "kernel compute",
    "halo sync",
    "halo transfer",
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
DEFAULT_PIN = "core_spread"
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
    (192, 1, 192),
    (225, 1, 225),
    (288, 1, 288),
    (360, 1, 360),
    (400, 1, 400),
    (480, 1, 480)
]

CASE_CUSTOM_FULL = [
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
    (400, 1, 400)
]

CASE_SQUARES_100 = [
    (1, 1, 1),
    (4, 2, 2),
    (9, 3, 3),
    (16, 4, 4),
    (25, 5, 5),
    (36, 6, 6),
    (49, 7, 7),
    (64, 8, 8),
    (81, 9, 9),
    (100, 10, 10),
]

CASE_SQUARES_SLICES_2520 = [
    # 1D horizontal slice (ix=1)
    (4, 1, 4),
    (9, 1, 9),
    (36, 1, 36),
    # 2D square (ix=iy)
    (4, 2, 2),
    (9, 3, 3),
    (36, 6, 6),
    # 1D vertical slice (iy=1)
    (4, 4, 1),
    (9, 9, 1),
    (36, 36, 1)
]

CASE_SQUARES_SLICES_3600 = [
    # 1D horizontal slice (ix=1)
    (4, 1, 4),
    (9, 1, 9),
    (16, 1, 16),
    (25, 1, 25),
    (36, 1, 36),
    (100, 1, 100),
    (144, 1, 144),
    (225, 1, 225),
    (400, 1, 400),
    # 2D square (ix=iy)
    (4, 2, 2),
    (9, 3, 3),
    (16, 4, 4),
    (25, 5, 5),
    (36, 6, 6),
    (100, 10, 10),
    (144, 12, 12),
    (225, 15, 15),
    (400, 20, 20),
    # 1D vertical slice (iy=1)
    (4, 4, 1),
    (9, 9, 1),
    (16, 16, 1),
    (25, 25, 1),
    (36, 36, 1),
    (100, 100, 1),
    (144, 144, 1),
    (225, 225, 1),
    (400, 400, 1)
]

CASE_SQUARES_SLICES_14400 = [
    # 1D horizontal slice (ix=1)
    (4, 1, 4),
    (9, 1, 9),
    (16, 1, 16),
    (25, 1, 25),
    (36, 1, 36),
    (64, 1, 64),
    (100, 1, 100),
    (144, 1, 144),
    (225, 1, 225),
    (400, 1, 400),
    # 2D square (ix=iy)
    (4, 2, 2),
    (9, 3, 3),
    (16, 4, 4),
    (25, 5, 5),
    (36, 6, 6),
    (64, 8, 8),
    (100, 10, 10),
    (144, 12, 12),
    (225, 15, 15),
    (400, 20, 20),
    # 1D vertical slice (iy=1)
    (4, 4, 1),
    (9, 9, 1),
    (16, 16, 1),
    (25, 25, 1),
    (36, 36, 1),
    (64, 64, 1),
    (100, 100, 1),
    (144, 144, 1),
    (225, 225, 1),
    (400, 400, 1)
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

# all 1D vertical slices (iy=1)
CASE_ALL_X = [
    (1, 1, 1),
    (2, 2, 1),
    (3, 3, 1),
    (4, 4, 1),
    (5, 5, 1),
    (6, 6, 1),
    (8, 8, 1),
    (9, 9, 1),
    (10, 10, 1),
    (12, 12, 1),
    (15, 15, 1),
    (16, 16, 1),
    (18, 18, 1),
    (20, 20, 1),
    (24, 24, 1),
    (25, 25, 1),
    (30, 30, 1),
    (32, 32, 1),
    (36, 36, 1),
    (40, 40, 1),
    (45, 45, 1),
    (48, 48, 1),
    (50, 50, 1),
    (60, 60, 1),
    (64, 64, 1),
    (72, 72, 1),
    (75, 75, 1),
    (80, 80, 1),
    (90, 90, 1),
    (96, 96, 1),
    (100, 100, 1),
    (120, 120, 1),
    (144, 144, 1),
    (150, 150, 1),
    (160, 160, 1),
    (180, 180, 1),
    (192, 192, 1),
    (200, 200, 1),
    (225, 225, 1),
    (240, 240, 1),
    (288, 288, 1),
    (300, 300, 1),
    (320, 320, 1),
    (360, 360, 1),
    (400, 400, 1),
    (450, 450, 1),
    (480, 480, 1),
    (576, 576, 1),
    (600, 600, 1),
    (720, 720, 1),
    (800, 800, 1),
    (900, 900, 1),
    (960, 960, 1),
    (1200, 1200, 1),
    (1440, 1440, 1),
    (1600, 1600, 1),
    (1800, 1800, 1),
    (2400, 2400, 1),
    (2880, 2880, 1),
    (3600, 3600, 1),
    (4800, 4800, 1),
    (7200, 7200, 1),
    (14400, 14400, 1)
]

# all 1D horizontal slices (ix=1)
CASE_ALL_Y = [
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
    "custom_full": CASE_CUSTOM_FULL,
    "squares_100": CASE_SQUARES_100,
    "squares_slices_2520": CASE_SQUARES_SLICES_2520,
    "squares_slices_3600": CASE_SQUARES_SLICES_3600,
    "squares_slices_14400": CASE_SQUARES_SLICES_14400,
    "small": CASE_SMALL,
    "small_full": CASE_SMALL_FULL,
    "medium": CASE_MEDIUM,
    "medium_full": CASE_MEDIUM_FULL,
    "large": CASE_LARGE,
    "large_full": CASE_LARGE_FULL,
    "all_x": CASE_ALL_X,
    "all_y": CASE_ALL_Y,
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
    parser.add_argument("--per-host", action="store_true")
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


def parse_timing_spread(output):
    values = {}
    for name, best_seconds, best_image_id, worst_seconds, worst_image_id in TIMING_SPREAD_RE.findall(output):
        values[name.strip().lower()] = {
            "best_seconds": float(best_seconds),
            "best_image_id": int(best_image_id),
            "worst_seconds": float(worst_seconds),
            "worst_image_id": int(worst_image_id),
        }

    missing = [name for name in TIMING_CATEGORIES if name not in values]
    if missing:
        raise RuntimeError("could not parse timing spread table")

    return values


def apply_pinning_preset(env, pin):
    for name in PINNING_ENV_NAMES:
        env.pop(name, None)

    env.update(PINNING_PRESETS[pin])


def print_pinning_settings(pin):
    env = " ".join(f"{name}={value}" for name, value in PINNING_PRESETS[pin].items())

    print_param("mpi pinning preset", pin)
    print_param("mpi pinning env", env)


def run_once(exe, images, ix, iy, run_num, pin, per_host=False):
    env = os.environ.copy()
    apply_pinning_preset(env, pin)
    if "FOR_COARRAY_CONFIG_FILE" not in env:
        env["FOR_COARRAY_NUM_IMAGES"] = str(images)
    env["I_X"] = str(ix)
    env["I_Y"] = str(iy)

    if per_host:
        # clear Slurm-level process layout controls to avoid conflicts
        env.pop("I_MPI_HYDRA_BOOTSTRAP", None)
        env.pop("SLURM_DISTRIBUTION", None)

        # force Intel MPI to fail fast and clean up resources on exit/failure
        env["I_MPI_FORCE_CLEANUP"] = "yes"
        env["I_MPI_FAIL_FAST"] = "yes"

        # generate dynamic hostfile if running in Slurm
        slurm_nodelist = env.get("SLURM_JOB_NODELIST")
        if slurm_nodelist:
            try:
                res = subprocess.run(
                    ["scontrol", "show", "hostnames", slurm_nodelist],
                    capture_output=True,
                    text=True,
                    check=True
                )
                hosts = res.stdout.strip().splitlines()
            except Exception:
                hosts = []

            if hosts:
                num_nodes = len(hosts)
                base = images // num_nodes
                rem = images % num_nodes

                lines = []
                for i, host in enumerate(hosts):
                    slots = base + (1 if i < rem else 0)
                    if slots > 0:
                        lines.append(f"{host}:{slots}")

                import tempfile
                tf = tempfile.NamedTemporaryFile(
                    mode="w", delete=False, prefix=f"hostfile_{images}_", dir="."
                )
                try:
                    tf.write("\n".join(lines) + "\n")
                    tf.close()
                    env["I_MPI_HYDRA_HOST_FILE"] = tf.name
                except Exception as e:
                    print(f"Warning: Failed to write dynamic hostfile: {e}")

    completed = subprocess.run(
        [exe],
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    output = completed.stdout + "\n" + completed.stderr

    # clean up temporary hostfile
    temp_hostfile_path = env.get("I_MPI_HYDRA_HOST_FILE")
    if temp_hostfile_path and os.path.exists(temp_hostfile_path):
        try:
            os.unlink(temp_hostfile_path)
        except Exception:
            pass

    if completed.returncode != 0:
        print(output)
        raise RuntimeError(f"run {run_num} failed with exit code {completed.returncode}")

    return (
        parse_last(STEP_RE, output, "step time"),
        parse_last(MLUPS_RE, output, "MLUPS"),
        parse_timing_spread(output),
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


def get_median_timing_measurement(timing_spreads, category, seconds_key, image_key):
    measurements = sorted(
        (values[category][seconds_key], values[category][image_key])
        for values in timing_spreads
    )

    return measurements[len(measurements) // 2]


def format_timing_measurement(seconds, image_id, width):
    return f"{seconds:.3f} ({image_id})".rjust(width)


def print_timing_spread_medians(timing_spreads):
    print()
    print("image execution time spread  |      best [sec] (image) |     worst [sec] (image)")
    print("-" * HEADER_WIDTH)

    for category in TIMING_CATEGORIES:
        best_seconds, best_image_id = get_median_timing_measurement(
            timing_spreads, category, "best_seconds", "best_image_id")
        worst_seconds, worst_image_id = get_median_timing_measurement(
            timing_spreads, category, "worst_seconds", "worst_image_id")

        print(
            f"{category:<28} | "
            f"{format_timing_measurement(best_seconds, best_image_id, 23)} | "
            f"{format_timing_measurement(worst_seconds, worst_image_id, 23)}"
        )


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


def run_case(exe, runs, n_x, n_y, case_num, n_cases, images, ix, iy, pin, per_host=False):
    validate_case(case_num, n_x, n_y, images, ix, iy)

    mlups_values = []
    timing_spreads = []

    print()
    print_header("benchmark script settings")
    print_param("case", f"{case_num} / {n_cases}")
    print_param("executable", exe)
    print_param("runs", runs)
    print_param("images", images)
    print_param("image grid", f"{ix} x {iy}")
    print_param("sim size", f"{n_x} x {n_y}")
    print_pinning_settings(pin)
    print_param("mpi per host option", "enabled" if per_host else "disabled")
    print()

    runs_started_at = timestamp()

    for run_num in range(1, runs + 1):
        if run_num > 1:
            time.sleep(3.0) # extra time to clean up resources and cool down

        step_ms, mlups, timing_spread, output = run_once(exe, images, ix, iy, run_num, pin, per_host)
        total_seconds = timing_spread["total"]["worst_seconds"]

        if run_num == 1:
            print_static_app_output(output)
            print_header(f"benchmark runs started at {runs_started_at}")

        print(
            f"{run_num:03d} | avg step time: {step_ms:.3f} ms | "
            f"total time: {total_seconds:.3f} sec | MLUPS: {int(mlups)}"
        )

        mlups_values.append(mlups)
        timing_spreads.append(timing_spread)

    print()
    mlups_stats = get_stats(mlups_values, higher_is_better=True)

    print_header("MLUPS metrics")
    print_param("median", f"{int(mlups_stats['median'])}")
    print_param("best", f"{int(mlups_stats['best'])}")
    print_param("worst", f"{int(mlups_stats['worst'])}")
    print_param("mean", f"{int(mlups_stats['mean'])}")
    print_param("stddev", f"{mlups_stats['stddev_percent']:.3f} %")

    print_timing_spread_medians(timing_spreads)


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
            args.exe, args.runs, n_x, n_y, case_num, len(domain_decomp_cases), images, ix, iy, args.pin, args.per_host)


if __name__ == "__main__":
    main()
