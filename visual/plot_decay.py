#!/usr/bin/env python3
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



# --- [ plot measured shear wave decay against analytical solution ] ---

# run config
RUN_NAME = "run_004_SW"
DATA_NAME = "velocity_x"

PLOT_NAME = "shear_wave_decay"
_METADATA_KEYS = ["rho_0", "omega", "u_max", "n_sin"]

MEASUREMENT_COLOR = "#e78ac3"
ANALYTICAL_COLOR = "black"
MEASUREMENT_INTERVAL = 50_000
MEASUREMENT_MARKER_SIZE = 6.0
ANALYTICAL_LINE_WIDTH = 2.5

# step config
STEP_START = 0
STEP_END = None       # None -> uses N_STEPS from config.json
STEP_STRIDE = None    # None -> uses export_interval from config.json

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "output" / RUN_NAME
PLOT_DIR = SCRIPT_DIR / "shear_wave" / RUN_NAME / "decay"


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


def get_wave_number(config: dict) -> float:
    return 2.0 * np.pi * config["n_sin"] / config["N_Y"]


def get_kinematic_viscosity(config: dict) -> float:
    return (1.0 / config["omega"] - 0.5) / 3.0


def get_mode_amplitude(field: np.ndarray, config: dict) -> float:
    velocity_x_profile = np.mean(field, axis=1)
    y = np.arange(config["N_Y"], dtype=np.float64)
    mode = np.sin(get_wave_number(config) * y)
    normalization = np.dot(mode, mode)

    if normalization == 0.0:
        raise ValueError("shear wave mode has zero normalization")

    return float(np.dot(velocity_x_profile, mode) / normalization)


def get_analytical_amplitudes(steps: np.ndarray, config: dict) -> np.ndarray:
    wave_number = get_wave_number(config)
    viscosity = get_kinematic_viscosity(config)
    return config["u_max"] * np.exp(-viscosity * wave_number * wave_number * steps)


def load_measured_amplitudes(steps: list[int], config: dict) -> tuple[np.ndarray, np.ndarray]:
    measured_steps = []
    measured_amplitudes = []

    for step in steps:
        data_path = get_data_path(step)
        if not data_path.exists():
            print(f"skipped step {step:>9}: missing {data_path.name}")
            continue

        field = load_field(data_path, config)
        measured_steps.append(step)
        measured_amplitudes.append(get_mode_amplitude(field, config))

    if not measured_steps:
        raise FileNotFoundError(f"no {DATA_NAME} fields found in {RUN_DIR}")

    return np.asarray(measured_steps), np.asarray(measured_amplitudes)


def plot_decay(steps: np.ndarray, measured_amplitudes: np.ndarray, config: dict) -> None:
    analytical_amplitudes = get_analytical_amplitudes(steps, config)

    fig, ax = plt.subplots(figsize=(8, 5.55))
    fig.subplots_adjust(bottom=0.13)

    measurement_mask = steps % MEASUREMENT_INTERVAL == 0
    ax.plot(
        steps,
        analytical_amplitudes,
        linewidth=ANALYTICAL_LINE_WIDTH,
        color=ANALYTICAL_COLOR,
        label=r"analytical: $A(t)=A_0 e^{-\nu k^2 t}$",
    )
    ax.plot(
        steps[measurement_mask],
        measured_amplitudes[measurement_mask],
        "o",
        markersize=MEASUREMENT_MARKER_SIZE,
        color=MEASUREMENT_COLOR,
        label="measurement",
    )
    for step, measurement_amplitude, analytical_amplitude in zip(
        steps[measurement_mask],
        measured_amplitudes[measurement_mask],
        analytical_amplitudes[measurement_mask],
    ):
        percentage_difference = 100.0 * (measurement_amplitude - analytical_amplitude) / analytical_amplitude
        ax.annotate(
            f"{percentage_difference:+.5f}%",
            xy=(step, measurement_amplitude),
            xytext=(3, 3),
            textcoords="offset points",
            rotation=45,
            ha="left",
            va="bottom",
            fontsize=6,
            color="0.2",
        )
    ax.set_title("shear wave amplitude decay")
    ax.set_xlabel("timestep")
    ax.set_ylabel("amplitude")
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

    print(f"run directory:  {RUN_DIR}")
    print(f"plot directory: {PLOT_DIR}")
    print(f"data field:     {DATA_NAME}")
    print(f"steps:          {steps}")
    print(f"wave number:    {get_wave_number(config):.6g}")
    print(f"viscosity:      {get_kinematic_viscosity(config):.6g}")

    if not is_data_exported(config):
        print(f"warning: config does not mark {DATA_NAME!r} as exported")

    measured_steps, measured_amplitudes = load_measured_amplitudes(steps, config)
    plot_decay(measured_steps, measured_amplitudes, config)
