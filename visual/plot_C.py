#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot velocity magnitude as streamlines ] ---

# run config
# options "shear_wave" (1), "couette_flow" (2), "poiseuille_flow" (3), "sliding_lid" (4)
SIM_MODE = 2
RUN_NAME = "run_001_CF"
DATA_NAME_X = "velocity_x"

_MODE_MAP = {1: "shear_wave", 2: "couette_flow", 3: "poiseuille_flow", 4: "sliding_lid"}
_SIM_MODE_RESOLVED = _MODE_MAP.get(SIM_MODE, SIM_MODE)
DATA_NAME_Y = "velocity_y"
PLOT_NAME = "velocity_mag_streamlines"

# step config
STEP_START = 0
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = None    # None -> uses export_interval from config.json
COLOR_LIMIT = None    # None -> uses max(|u|) across all selected steps
STREAM_STRIDE = 5
STREAM_DENSITY = 2.5
STREAM_LINEWIDTH = 1.5
STREAM_ARROWSIZE = 1.5

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = SCRIPT_DIR / _SIM_MODE_RESOLVED / RUN_NAME / "C"


def format_step_suffix(step: int, width: int = 9) -> str:
    return f"_{step:0{width}d}"


def load_config() -> dict:
    config_path = RUN_DIR / "config.json"
    with config_path.open("r", encoding="utf-8") as file:
        return json.load(file)


def get_steps(config: dict) -> list[int]:
    step_start = STEP_START
    step_end = config["N_STEPS"] if STEP_END is None else STEP_END
    step_stride = config["export_interval"] if STEP_STRIDE is None else STEP_STRIDE

    steps = list(range(step_start, step_end + 1, step_stride))
    if config["export_final_state"] and step_end not in steps:
        steps.append(step_end)

    return steps


def get_data_path(data_name: str, step: int) -> Path:
    return RUN_DIR / f"{data_name}{format_step_suffix(step)}.bin"


def get_file_dtype(config: dict) -> np.dtype:
    file_dtype = config.get("file_dtype", "real32")
    if file_dtype == "real32":
        return np.float32
    if file_dtype == "real64":
        return np.float64

    raise ValueError(f"unsupported file_dtype: {file_dtype}")


def load_field(path: Path, config: dict) -> np.ndarray:
    N_X = config["N_X"]
    N_Y = config["N_Y"]
    return np.fromfile(path, dtype=get_file_dtype(config)).reshape((N_Y, N_X))


def is_data_exported(config: dict) -> bool:
    return bool(config.get("export_u_x", False)) and bool(config.get("export_u_y", False))


def get_color_limit(steps: list[int], config: dict) -> float:
    if COLOR_LIMIT is not None:
        return float(COLOR_LIMIT)

    max_val = 0.0
    for step in steps:
        u_x_path = get_data_path(DATA_NAME_X, step)
        u_y_path = get_data_path(DATA_NAME_Y, step)
        if u_x_path.exists() and u_y_path.exists():
            u_x = load_field(u_x_path, config)
            u_y = load_field(u_y_path, config)
            velocity_mag = np.sqrt(u_x * u_x + u_y * u_y)
            max_val = max(max_val, float(np.max(velocity_mag)))

    return max(max_val, 1.0e-12)


def plot_step(step: int, config: dict, color_limit: float) -> None:
    u_x_path = get_data_path(DATA_NAME_X, step)
    u_y_path = get_data_path(DATA_NAME_Y, step)
    missing_paths = [path.name for path in [u_x_path, u_y_path] if not path.exists()]
    if missing_paths:
        print(f"skipped step {step:>9}: missing {', '.join(missing_paths)}")
        return

    u_x = load_field(u_x_path, config)
    u_y = load_field(u_y_path, config)
    velocity_mag = np.sqrt(u_x * u_x + u_y * u_y)

    x_grid = np.arange(config["N_X"], dtype=get_file_dtype(config))[::STREAM_STRIDE]
    y_grid = np.arange(config["N_Y"], dtype=get_file_dtype(config))[::STREAM_STRIDE]
    u_x_plot = u_x[::STREAM_STRIDE, ::STREAM_STRIDE]
    u_y_plot = u_y[::STREAM_STRIDE, ::STREAM_STRIDE]
    velocity_mag_plot = velocity_mag[::STREAM_STRIDE, ::STREAM_STRIDE]

    fig, ax = plt.subplots(figsize=(8, 5))
    stream = ax.streamplot(
        x_grid,
        y_grid,
        u_x_plot,
        u_y_plot,
        color=velocity_mag_plot,
        cmap="turbo",
        norm=matplotlib.colors.Normalize(vmin=0.0, vmax=color_limit),
        density=STREAM_DENSITY,
        linewidth=STREAM_LINEWIDTH,
        arrowsize=STREAM_ARROWSIZE,
    )

    fig.colorbar(stream.lines, ax=ax, label="|u|")
    ax.set_title(f"|u| streamlines at step {step}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(0, config["N_X"])
    ax.set_ylim(0, config["N_Y"])
    ax.set_aspect("equal")

    output_path = PLOT_DIR / f"{PLOT_NAME}{format_step_suffix(step)}.png"
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    config = load_config()
    steps = get_steps(config)
    color_limit = get_color_limit(steps, config)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"run directory:  {RUN_DIR}")
    print(f"plot directory: {PLOT_DIR}")
    print(f"data field:     {PLOT_NAME}")
    print(f"steps:          {steps}")
    print(f"color bounds:   [0, {color_limit:.6g}]")

    if not is_data_exported(config):
        print("warning: config does not mark velocity_x and velocity_y as exported")

    for step in steps:
        plot_step(step, config, color_limit)
