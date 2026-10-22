#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np



# --- [ plot one decomposition scaling benchmark result ] ---

# summary config
JOB_NAME = "A"  # options "A" or "B"
JOB_ID = None  # if "None", uses latest summary
HORIZONTAL_COLOR = (30/255, 140/255, 80/255)
SQUARED_COLOR = (50/255, 80/255, 170/255)
VERTICAL_COLOR = (160/255, 110/255, 30/255)
LAYOUTS = (
    ("1 x N decompositions", "1 x N (horizontal)", HORIZONTAL_COLOR),
    ("N x N decompositions", "M x M (squared)", SQUARED_COLOR),
    ("N x 1 decompositions", "N x 1 (vertical)", VERTICAL_COLOR),
)
INCLUDE_PARTIAL_MEASUREMENTS = True

# plot config
FIGURE_SIZE = (14, 6)
FIGURE_DPI = 600
GRID_ALPHA = 0.5
MEASUREMENT_MARKER_SIZE = 12.0
MEASUREMENT_LINE_WIDTH = 3.0
MEASUREMENT_LINE_STYLE = "-"

IDEAL_LINE_COLOR = "black"
IDEAL_LINE_STYLE = ":"
IDEAL_LINE_WIDTH = 3.0
IDEAL_BASELINE_IMAGES = 4
LEGEND_LOCATION = "lower right"

# font size config
X_AXIS_FONT_SIZE = 16
Y_AXIS_FONT_SIZE = 16
LEGEND_FONT_SIZE = 16

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[2]
SUMMARY_DIRS = {
    "A": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_A",
    "B": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_B",
}
SCRIPT_STEM = Path(__file__).stem
PLOT_DIR = SCRIPT_DIR / "plots"
PLOT_NAME = None  # None -> <SCRIPT_STEM>_<JOB_NAME>[_<JOB_ID>]


def resolve_summary_path(job_name: str, job_id: int | str | None = None) -> Path:
    if job_name not in SUMMARY_DIRS:
        raise ValueError(f"unsupported JOB_NAME: {job_name!r}")

    job_dir = SUMMARY_DIRS[job_name]
    if not job_dir.exists():
        raise FileNotFoundError(f"missing job directory: {job_dir}")

    if job_id is not None:
        target_path = job_dir / f"summary_{job_id}.md"
        if not target_path.exists():
            raise FileNotFoundError(f"missing summary for job id {job_id}: {target_path}")
        return target_path

    summary_files = list(job_dir.glob("summary_*.md"))
    if not summary_files:
        raise FileNotFoundError(f"no summary_*.md files found in {job_dir}")

    def extract_job_id(path: Path) -> int:
        parts = path.stem.split("_")
        if len(parts) >= 2 and parts[1].isdigit():
            return int(parts[1])
        return 0

    return max(summary_files, key=extract_job_id)


def split_table_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def load_summary(path: Path) -> dict[str, list[dict[str, object]]]:
    rows = {layout_name: [] for _, layout_name, _ in LAYOUTS}
    layout_by_heading = {heading: layout_name for heading, layout_name, _ in LAYOUTS}
    layout_by_heading["M x M decompositions"] = "M x M (squared)"
    layout_by_heading["N x N decompositions"] = "M x M (squared)"
    current_layout = None

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped_line = line.strip()

        if stripped_line.startswith("## "):
            current_layout = layout_by_heading.get(stripped_line[3:])
            continue

        if not stripped_line.startswith("|") or not stripped_line.endswith("|"):
            continue

        cells = split_table_row(stripped_line)
        if current_layout is None or len(cells) < 9 or cells[0] == "images":
            continue

        try:
            images = int(cells[0])
            mlups = int(cells[5])
        except ValueError:
            continue

        rows[current_layout].append(
            {
                "images": images,
                "mlups": mlups,
                "status": cells[8],
            }
        )

    return rows


def get_plot_rows(layout_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if INCLUDE_PARTIAL_MEASUREMENTS:
        return layout_rows

    return [row for row in layout_rows if row["status"] == "complete"]


def get_plot_name(summary_path: Path) -> str:
    if PLOT_NAME is not None:
        return PLOT_NAME

    if JOB_ID is not None:
        return f"{SCRIPT_STEM}_{JOB_NAME}_{JOB_ID}"

    return f"{SCRIPT_STEM}_{JOB_NAME}"


def plot_decomposition_scaling(
    ax,
    rows: dict[str, list[dict[str, object]]],
) -> None:
    for _, layout_name, color in LAYOUTS:
        layout_rows = get_plot_rows(rows[layout_name])
        if not layout_rows:
            continue

        layout_images = np.asarray([row["images"] for row in layout_rows])
        layout_mlups = np.asarray([row["mlups"] for row in layout_rows])
        ax.plot(
            layout_images,
            layout_mlups,
            marker="o",
            markersize=MEASUREMENT_MARKER_SIZE,
            linewidth=MEASUREMENT_LINE_WIDTH,
            linestyle=MEASUREMENT_LINE_STYLE,
            color=color,
            label=layout_name,
        )

    # capture limits from measured curves so ideal line does not expand them
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # calculate and plot ideal scaling line ('linear scaling') from 4-image baseline
    all_baseline_mlups = [
        row["mlups"]
        for layout_rows in rows.values()
        for row in get_plot_rows(layout_rows)
        if row["images"] == IDEAL_BASELINE_IMAGES
    ]
    if all_baseline_mlups:
        baseline_mlups = np.mean(all_baseline_mlups)
        all_images = [
            row["images"]
            for layout_rows in rows.values()
            for row in get_plot_rows(layout_rows)
        ]
        if all_images:
            x_min, x_max = min(all_images), max(all_images)
            x_linear = np.linspace(x_min, x_max, 200)
            y_linear = baseline_mlups * (x_linear / float(IDEAL_BASELINE_IMAGES))

            ax.plot(
                x_linear,
                y_linear,
                color=IDEAL_LINE_COLOR,
                linestyle=IDEAL_LINE_STYLE,
                linewidth=IDEAL_LINE_WIDTH,
                label="linear scaling",
                zorder=1,
            )

    # restore limits so linear line may run off-screen
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.tick_params(axis="x", labelsize=X_AXIS_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=Y_AXIS_FONT_SIZE)
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.grid(alpha=GRID_ALPHA)
    ax.legend(
        loc=LEGEND_LOCATION,
        ncol=1,
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
    )


def plot_selected_job() -> None:
    summary_path = resolve_summary_path(JOB_NAME, JOB_ID)

    rows = load_summary(summary_path)
    if not any(get_plot_rows(layout_rows) for layout_rows in rows.values()):
        raise ValueError(f"no measurements found in {summary_path}")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plot_decomposition_scaling(ax, rows)

    for layout_name in rows:
        partial_rows = [row for row in rows[layout_name] if row["status"] != "complete"]
        if partial_rows:
            partial_images = ", ".join(str(row["images"]) for row in partial_rows)
            print(f"warning: {JOB_NAME} {layout_name} includes incomplete measurements at {partial_images} images")

    output_path = PLOT_DIR / f"{get_plot_name(summary_path)}.png"
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    summary_path = resolve_summary_path(JOB_NAME, JOB_ID)
    print(f"plot directory: {PLOT_DIR}")
    print(f"summary:        {summary_path}")
    print(f"job:            bench_decomp_{JOB_NAME} (job id: {JOB_ID or summary_path.stem.replace('summary_', '')})")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_selected_job()

