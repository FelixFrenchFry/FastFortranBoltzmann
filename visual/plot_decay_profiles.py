#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot stacked 1d velocity profiles over time ] ---

# run config
RUN_NAME = "run_004_SW"
DATA_NAME = "velocity_x"

PLOT_NAME = "shear_wave_decay_profiles"
_METADATA_KEYS = ["rho_0", "omega", "u_max", "n_sin"]

PROFILE_INTERVAL = 100_000
PROFILE_LINE_WIDTH = 2.0
LEGEND = True

# step config
STEP_START = 0
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = None    # None -> uses export_interval from config.json

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = SCRIPT_DIR / "shear_wave" / RUN_NAME / "decay_profiles"


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
    export_endpoint_states = config.get("export_endpoint_states", config.get("export_final_state", False))
    if export_endpoint_states and step_end not in steps:
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
    if "export_macros" in config:
        return bool(config["export_macros"])

    return bool(config.get("export_u_x", False))


def format_metadata_value(value) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.6f}"

    return str(value)


def format_metadata_entries(config: dict, keys: list[str]) -> list[str]:
    return [f"{key}={format_metadata_value(config[key])}" for key in keys if key in config]


def get_metadata_text(config: dict) -> str:
    sim_entries = format_metadata_entries(config, _METADATA_KEYS)
    grid_entries = [
        f"N_X_TOTAL={format_metadata_value(config['N_X'])}",
        f"N_Y_TOTAL={format_metadata_value(config['N_Y'])}",
    ]

    return " | ".join(sim_entries + grid_entries)


def add_metadata_text(fig, config: dict):
    metadata_text = get_metadata_text(config)
    if not metadata_text:
        return None

    return fig.text(
        0.5,
        0.018,
        metadata_text,
        ha="center",
        va="bottom",
        fontsize=6.5,
        color="0.15",
    )


def load_sampled_profiles(steps: list[int], config: dict) -> tuple[list[int], list[np.ndarray]]:
    sampled_steps = []
    profiles = []

    for step in steps:
        if step % PROFILE_INTERVAL != 0 and step != steps[-1]:
            continue

        data_path = get_data_path(step)
        if not data_path.exists():
            print(f"skipped step {step:>9}: missing {data_path.name}")
            continue

        field = load_field(data_path, config)
        profile = np.mean(field, axis=1)
        sampled_steps.append(step)
        profiles.append(profile)

    if not sampled_steps:
        raise FileNotFoundError(f"no {DATA_NAME} fields found in {RUN_DIR}")

    return sampled_steps, profiles


def plot_decay_profiles(sampled_steps: list[int], profiles: list[np.ndarray], config: dict) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.55))
    fig.subplots_adjust(bottom=0.13)

    y_grid = np.arange(config["N_Y"], dtype=np.float64) + 0.5
    colors = plt.cm.viridis(np.linspace(0, 1, len(sampled_steps)))

    for step, profile, color in zip(sampled_steps, profiles, colors):
        ax.plot(
            profile,
            y_grid,
            linewidth=PROFILE_LINE_WIDTH,
            color=color,
            label=f"t = {step:,}",
        )

    ax.set_xlabel("u_x")
    ax.set_ylabel("y")
    ax.set_ylim(0, config["N_Y"])
    ax.grid(alpha=0.25)
    legend_obj = None
    if LEGEND:
        legend_obj = ax.legend(loc="upper right", fontsize=8)

    metadata_text = add_metadata_text(fig, config)
    extra_artists = [artist for artist in [metadata_text, legend_obj] if artist is not None]

    output_path = PLOT_DIR / f"{PLOT_NAME}.png"
    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
        bbox_extra_artists=extra_artists if extra_artists else None,
    )
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    config = load_config()
    steps = get_steps(config)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"run directory:  {RUN_DIR}")
    print(f"plot directory: {PLOT_DIR}")
    print(f"data field:     {DATA_NAME}")
    print(f"steps:          {steps}")

    if not is_data_exported(config):
        print(f"warning: config does not mark {DATA_NAME!r} as exported")

    sampled_steps, profiles = load_sampled_profiles(steps, config)
    plot_decay_profiles(sampled_steps, profiles, config)
