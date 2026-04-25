#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot X velocity scalar field ] ---

# run config
RUN_NAME = "run_001"
DATA_NAME = "velocity_x"

# step config
STEP_START = 0
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = None    # None -> uses export_interval from config.json
COLOR_LIMIT = None    # None -> uses max(abs(field)) across all selected steps

# path config
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = ROOT_DIR / "visual" / "plots" / RUN_NAME


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


def get_data_path(step: int) -> Path:
    return RUN_DIR / f"{DATA_NAME}{format_step_suffix(step)}.bin"


def load_field(path: Path, config: dict) -> np.ndarray:
    N_X = config["N_X"]
    N_Y = config["N_Y"]
    return np.fromfile(path, dtype=np.float32).reshape((N_Y, N_X))


def is_data_exported(config: dict) -> bool:
    export_key_by_data_name = {
        "density": "export_rho",
        "velocity_x": "export_u_x",
        "velocity_y": "export_u_y",
        "velocity_mag": "export_u_mag",
    }
    export_key = export_key_by_data_name[DATA_NAME]
    return bool(config.get(export_key, False))


def get_color_limit(steps: list[int], config: dict) -> float:
    if COLOR_LIMIT is not None:
        return float(COLOR_LIMIT)

    max_abs = 0.0
    for step in steps:
        data_path = get_data_path(step)
        if data_path.exists():
            field = load_field(data_path, config)
            max_abs = max(max_abs, float(np.max(np.abs(field))))

    return max(max_abs, 1.0e-12)


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
        cmap="seismic",
        vmin=-color_limit,
        vmax=color_limit,
        extent=(0, config["N_X"], 0, config["N_Y"]),
        interpolation="nearest",
    )

    fig.colorbar(image, ax=ax, label="u_x")
    ax.set_title(f"u_x at step {step}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    output_path = PLOT_DIR / f"{DATA_NAME}{format_step_suffix(step)}.png"
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
    print(f"data field:     {DATA_NAME}")
    print(f"steps:          {steps}")
    print(f"color bounds:   [{-color_limit:.6g}, {color_limit:.6g}]")

    if not is_data_exported(config):
        print(f"warning: config does not mark {DATA_NAME!r} as exported")

    for step in steps:
        plot_step(step, config, color_limit)
