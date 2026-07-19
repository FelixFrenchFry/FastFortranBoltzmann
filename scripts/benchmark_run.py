#!/usr/bin/env python3
"""Run a sim config multiple times and summarize performance metrics"""

import argparse
import os
import re
import signal
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
DEFAULT_TIMEOUT = 15 * 60
TIMEOUT_TERMINATION_GRACE = 30


class RunTimeoutError(RuntimeError):
    def __init__(self, run_num, timeout, output):
        self.output = output
        super().__init__(f"run {run_num} timed out after {timeout} seconds")


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
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--images", type=int, required=True)
    parser.add_argument("--ix", type=int, required=True)
    parser.add_argument("--iy", type=int, required=True)
    parser.add_argument("--pin", choices=PINNING_PRESETS.keys(), default=DEFAULT_PIN)
    parser.add_argument("--per-host", action="store_true")
    return parser.parse_args()


def terminate_process(process, force=False):
    if os.name == "posix":
        signal_number = signal.SIGKILL if force else signal.SIGTERM
        try:
            os.killpg(process.pid, signal_number)
        except ProcessLookupError:
            pass
    elif force:
        process.kill()
    else:
        process.terminate()


def run_process(command, env, timeout, run_num):
    process = subprocess.Popen(
        command,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=(os.name == "posix"),
    )

    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        terminate_process(process)

        try:
            stdout, stderr = process.communicate(timeout=TIMEOUT_TERMINATION_GRACE)
        except subprocess.TimeoutExpired:
            terminate_process(process, force=True)
            stdout, stderr = process.communicate()

        output = (stdout or "") + "\n" + (stderr or "")
        raise RunTimeoutError(run_num, timeout, output)

    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


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


def run_once(exe, images, ix, iy, run_num, pin, per_host=False, timeout=DEFAULT_TIMEOUT):
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

    try:
        completed = run_process([exe], env, timeout, run_num)
    finally:
        # clean up temporary hostfile
        temp_hostfile_path = env.get("I_MPI_HYDRA_HOST_FILE")
        if temp_hostfile_path and os.path.exists(temp_hostfile_path):
            try:
                os.unlink(temp_hostfile_path)
            except Exception:
                pass

    output = completed.stdout + "\n" + completed.stderr

    if completed.returncode != 0:
        print(output)
        raise RuntimeError(f"run {run_num} failed with exit code {completed.returncode}")

    try:
        step_ms = parse_last(STEP_RE, output, "step time")
        mlups = parse_last(MLUPS_RE, output, "MLUPS")
        timing_spread = parse_timing_spread(output)
    except Exception as e:
        print("--- [ DEBUG: APP OUTPUT ON PARSE FAILURE ] ---")
        print(output)
        print("---------------------------------------------")
        raise

    return (
        step_ms,
        mlups,
        timing_spread,
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
    measurements = [
        (values[category][seconds_key], values[category][image_key])
        for values in timing_spreads
    ]
    measurements.sort(key=lambda measurement: measurement[0])

    return measurements[len(measurements) // 2]


def format_timing_measurement(seconds, image_id, width):
    image_text = "n/a" if image_id is None else str(image_id)
    return f"{seconds:.3f} ({image_text})".rjust(width)


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


def main():
    args = parse_args()

    if args.runs <= 0:
        sys.exit("error: --runs must be positive")
    if args.runs > MAX_RUNS:
        sys.exit(f"error: --runs must be <= {MAX_RUNS}")
    if args.timeout <= 0:
        sys.exit("error: --timeout must be positive")
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

    mlups_values = []
    timing_spreads = []

    print()
    print_header("benchmark script settings")
    print_param("executable", args.exe)
    print_param("runs", args.runs)
    print_param("timeout", f"{args.timeout} sec")
    print_param("images", args.images)
    print_param("image grid", f"{args.ix} x {args.iy}")
    print_param("sim size", f"{n_x} x {n_y}")
    print_pinning_settings(args.pin)
    print_param("mpi per host option", "enabled" if args.per_host else "disabled")
    print()

    runs_started_at = timestamp()

    for run_num in range(1, args.runs + 1):
        if run_num > 1:
            time.sleep(3.0) # extra time to clean up resources and cool down

        try:
            step_ms, mlups, timing_spread, output = run_once(
                args.exe, args.images, args.ix, args.iy, run_num, args.pin, args.per_host, args.timeout)
        except RunTimeoutError as error:
            print(f"{run_num:03d} | timed out after {args.timeout} sec")
            if error.output.strip():
                print(error.output)
            continue
        total_seconds = timing_spread["total"]["worst_seconds"]

        if not mlups_values:
            print_static_app_output(output)
            print_header(f"benchmark runs started at {runs_started_at}")

        print(
            f"{run_num:03d} | avg step time: {step_ms:.3f} ms | "
            f"total time: {total_seconds:.3f} sec | MLUPS: {int(mlups)}"
        )

        mlups_values.append(mlups)
        timing_spreads.append(timing_spread)

    if not mlups_values:
        sys.exit("error: no benchmark runs completed")

    print()
    mlups_stats = get_stats(mlups_values, higher_is_better=True)

    print_header("MLUPS metrics")
    print_param("median", f"{int(mlups_stats['median'])}")
    print_param("best", f"{int(mlups_stats['best'])}")
    print_param("worst", f"{int(mlups_stats['worst'])}")
    print_param("mean", f"{int(mlups_stats['mean'])}")
    print_param("stddev", f"{mlups_stats['stddev_percent']:.3f} %")

    print_timing_spread_medians(timing_spreads)


if __name__ == "__main__":
    main()
