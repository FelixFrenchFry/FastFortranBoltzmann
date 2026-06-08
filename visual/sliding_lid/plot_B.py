#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot sliding lid velocity magnitude scalar field as heatmap ] ---

# run config
RUN_NAME = "run_001_SL"
DATA_NAME = "velocity_mag"

# step config
STEP_START = 0
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = None    # None -> uses export_interval from config.json
COLOR_LIMIT = None    # None -> uses max(field) across all selected steps

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = SCRIPT_DIR / "plots" / RUN_NAME / "B"


def format_step_suffix(step: int, width: int = 9) -> str:
    return f"_{step:0{width}d}"


def load_config() -> dict:
    config_path = RUN_DIR / "config.json"
    with config_path.open("r", encoding="utf-8") as file:
        return json.load(file)


def validate_config(config: dict) -> None:
    sim_mode = config.get("SIM_MODE")
    if sim_mode != "sliding_lid":
        raise ValueError(f"expected SIM_MODE 'sliding_lid', got {sim_mode!r}")


def get_steps(config: dict) -> list[int]:
    step_start = STEP_START
    step_end = config["N_STEPS"] if STEP_END is None else STEP_END
    step_stride = config["export_interval"] if STEP_STRIDE is None else STEP_STRIDE

    if step_stride <= 0:
        return []

    steps = list(range(step_start, step_end + 1, step_stride))
    if not config.get("export_initial_state", False):
        steps = [step for step in steps if step != 0]

    if config.get("export_final_state", False) and step_end not in steps:
        steps.append(step_end)

    return steps


def get_data_path(step: int) -> Path:
    return RUN_DIR / f"{DATA_NAME}{format_step_suffix(step)}.bin"


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
    return bool(config.get("export_u_mag", False))


def get_color_limit(steps: list[int], config: dict) -> float:
    if COLOR_LIMIT is not None:
        return float(COLOR_LIMIT)

    max_val = 0.0
    for step in steps:
        data_path = get_data_path(step)
        if data_path.exists():
            field = load_field(data_path, config)
            max_val = max(max_val, float(np.max(field)))

    return max(max_val, 1.0e-12)


def plot_step(step: int, config: dict, color_limit: float) -> None:
    data_path = get_data_path(step)
    if not data_path.exists():
        print(f"skipped step {step:>9}: missing {data_path.name}")
        return

    field = load_field(data_path, config)

    fig, ax = plt.subplots(figsize=(8, 5))
    image = ax.imshow(
        field,
        origin="lower",
        cmap="turbo",
        vmin=0.0,
        vmax=color_limit,
        extent=(0, config["N_X"], 0, config["N_Y"]),
        interpolation="nearest",
    )

    fig.colorbar(image, ax=ax, label="|u|")
    ax.set_title(f"|u| at step {step}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    output_path = PLOT_DIR / f"{DATA_NAME}{format_step_suffix(step)}.png"
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    config = load_config()
    validate_config(config)
    steps = get_steps(config)
    color_limit = get_color_limit(steps, config)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"run directory:  {RUN_DIR}")
    print(f"plot directory: {PLOT_DIR}")
    print(f"data field:     {DATA_NAME}")
    print(f"steps:          {steps}")
    print(f"color bounds:   [0, {color_limit:.6g}]")

    if not is_data_exported(config):
        print(f"warning: config does not mark {DATA_NAME!r} as exported")

    for step in steps:
        plot_step(step, config, color_limit)
