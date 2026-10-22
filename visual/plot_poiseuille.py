#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot measured poiseuille flow profile against analytical solution ] ---

# run config
RUN_NAME = "run_005_PF"
DATA_NAME = "velocity_x"

PLOT_NAME = "poiseuille_profile"
_METADATA_KEYS = ["rho_0", "omega", "rho_in", "rho_out"]

ANALYTICAL_COLOR = "black"
MEASUREMENT_STRIDE = 50
MEASUREMENT_INTERVAL = 200_000
MEASUREMENT_MARKER_SIZE = 6.0
ANALYTICAL_LINE_WIDTH = 2.5

# step config
STEP_START = MEASUREMENT_INTERVAL
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = MEASUREMENT_INTERVAL

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = SCRIPT_DIR / "poiseuille_flow" / RUN_NAME / "profile"


def format_step_suffix(step: int, width: int = 9) -> str:
    return f"_{step:0{width}d}"


def load_config() -> dict:
    config_path = RUN_DIR / "config.json"
    with config_path.open("r", encoding="utf-8") as file:
        return json.load(file)


def get_steps(config: dict) -> list[int]:
    step_start = STEP_START
    step_end = config["N_STEPS"] if STEP_END is None else STEP_END
    step_stride = STEP_STRIDE

    return list(range(step_start, step_end + 1, step_stride))


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


def get_kinematic_viscosity(config: dict) -> float:
    return (1.0 / config["omega"] - 0.5) / 3.0


def get_pressure_gradient(config: dict) -> float:
    speed_of_sound_squared = 1.0 / 3.0
    channel_length = config["N_X"]
    return speed_of_sound_squared * (config["rho_in"] - config["rho_out"]) / channel_length


def get_analytical_velocity_profile(config: dict) -> tuple[np.ndarray, np.ndarray]:
    y = np.arange(config["N_Y"], dtype=np.float64) + 0.5
    channel_height = float(config["N_Y"])
    viscosity = get_kinematic_viscosity(config)
    pressure_gradient = get_pressure_gradient(config)
    velocity_x = pressure_gradient / (2.0 * viscosity) * y * (channel_height - y)

    return y, velocity_x


def get_measurement_mask(size: int) -> np.ndarray:
    measurement_mask = np.arange(size) % MEASUREMENT_STRIDE == 0
    measurement_mask[0] = True
    measurement_mask[-1] = True
    return measurement_mask


def load_measured_profiles(steps: list[int], config: dict) -> dict[int, np.ndarray]:
    measured_profiles = {}

    for step in steps:
        data_path = get_data_path(step)
        if not data_path.exists():
            print(f"skipped step {step:>9}: missing {data_path.name}")
            continue

        field = load_field(data_path, config)
        measured_profiles[step] = np.mean(field, axis=1)

    if not measured_profiles:
        raise FileNotFoundError(f"no {DATA_NAME} fields found in {RUN_DIR}")

    return measured_profiles


def plot_profile(measured_profiles: dict[int, np.ndarray], config: dict) -> None:
    y, analytical_profile = get_analytical_velocity_profile(config)
    measurement_mask = get_measurement_mask(y.size)
    measurement_colors = plt.cm.PuRd(np.linspace(0.35, 0.75, len(measured_profiles)))

    fig, ax = plt.subplots(figsize=(8, 5.55))
    fig.subplots_adjust(bottom=0.13)

    ax.plot(
        y,
        analytical_profile,
        linewidth=ANALYTICAL_LINE_WIDTH,
        color=ANALYTICAL_COLOR,
        label=r"analytical: $u_x(y)=\frac{G}{2\nu}y(H-y)$",
    )
    for (step, measured_profile), color in zip(measured_profiles.items(), measurement_colors):
        if measured_profile.shape != analytical_profile.shape:
            raise ValueError("measured and analytical profiles have different shapes")

        ax.plot(
            y[measurement_mask],
            measured_profile[measurement_mask],
            "o",
            markersize=MEASUREMENT_MARKER_SIZE,
            color=color,
            label=f"measurement (step {step})",
        )
    ax.set_title("poiseuille flow velocity profile")
    ax.set_xlabel("y")
    ax.set_ylabel("velocity x")
    ax.grid(alpha=0.25)
    ax.legend()
    metadata_text = add_metadata_text(fig, config)

    output_path = PLOT_DIR / f"{PLOT_NAME}.png"
    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
        bbox_extra_artists=[metadata_text] if metadata_text is not None else None,
    )
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    config = load_config()
    steps = get_steps(config)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"run directory:      {RUN_DIR}")
    print(f"plot directory:     {PLOT_DIR}")
    print(f"data field:         {DATA_NAME}")
    print(f"steps:              {steps}")
    print(f"pressure gradient:  {get_pressure_gradient(config):.6g}")
    print(f"viscosity:          {get_kinematic_viscosity(config):.6g}")

    if not is_data_exported(config):
        print(f"warning: config does not mark {DATA_NAME!r} as exported")

    measured_profiles = load_measured_profiles(steps, config)
    print(f"measurement steps:  {list(measured_profiles)}")
    plot_profile(measured_profiles, config)
