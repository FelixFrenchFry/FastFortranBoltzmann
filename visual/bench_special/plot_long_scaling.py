#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import numpy as np



# --- [ plot long decomposition scaling benchmark result ] ---

# summary config
JOB_NAMES = ("A", "B", "C")  # options "A", "B", or "C"
JOB_IDS = {
    "A": None,
    "B": None,
    "C": None,
}  # if "None", uses latest summary for each job
JOB_COLORS = {
    "A": (40/255, 80/255, 180/255),
    "B": (40/255, 180/255, 80/255),
    "C": (250/255, 0/255, 50/255),
}
INCLUDE_PARTIAL_MEASUREMENTS = True

# plot config
FIGURE_SIZE = (14, 6)
FIGURE_DPI = 600
GRID_ALPHA = 0.5
MEASUREMENT_MARKER_SIZE = 10.0
MEASUREMENT_LINE_WIDTH = 3.0
MEASUREMENT_LINE_STYLE = "-"
USE_LOG_SCALE = False
TASK_LINE_STYLE = "-"
TASK_LINE_WIDTH = 2.0
TASK_LINE_ALPHA = 0.8
X_AXIS_TICK_STEP = 100
X_LIMIT_TICK_STEP = 25
X_LIMIT = None  # None -> uses the full measured x range
Y_LIMIT_MARGIN = 0.05

IDEAL_LINE_COLOR = "black"
IDEAL_LINE_STYLE = ":"
IDEAL_LINE_WIDTH = 3.0
IDEAL_BASELINE_IMAGES = 4
LEGEND_LOCATION = "upper left"

# font size config
X_AXIS_FONT_SIZE = 16
Y_AXIS_FONT_SIZE = 16
LEGEND_FONT_SIZE = 16

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[2]
SUMMARY_DIRS = {
    "A": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "long_A",
    "B": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "long_B",
    "C": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "long_C",
}
SCRIPT_STEM = Path(__file__).stem
PLOT_DIR = SCRIPT_DIR / "plots"
PLOT_NAME = None  # None -> <SCRIPT_STEM>_<JOB_NAMES>


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


def load_summary(path: Path) -> dict[str, object]:
    metadata = {}
    rows = []
    in_section = False

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped_line = line.strip()

        if stripped_line.startswith("## "):
            in_section = (stripped_line[3:].strip() == "1 x N decompositions")
            continue

        if not stripped_line.startswith("|") or not stripped_line.endswith("|"):
            continue

        cells = split_table_row(stripped_line)
        if not in_section:
            if len(cells) >= 2 and cells[0] in ("nodes", "total tasks"):
                metadata[cells[0]] = cells[1]
            continue

        if len(cells) < 9 or cells[0] == "images":
            continue

        try:
            images = int(cells[0])
        except ValueError:
            continue

        try:
            mlups = int(cells[5])
        except ValueError:
            mlups = None

        rows.append(
            {
                "images": images,
                "mlups": mlups,
                "status": cells[8],
            }
        )

    return {
        "metadata": metadata,
        "rows": rows,
    }


def get_plot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows_with_measurements = [row for row in rows if row["mlups"] is not None]
    if INCLUDE_PARTIAL_MEASUREMENTS:
        return rows_with_measurements

    return [row for row in rows_with_measurements if row["status"] == "complete"]


def get_job_label(metadata: dict[str, str]) -> str:
    nodes = metadata.get("nodes")
    total_tasks = get_total_tasks(metadata)
    if nodes and nodes != "-" and total_tasks is not None:
        return f"{nodes} nodes / {total_tasks} tasks"
    if nodes and nodes != "-":
        return f"{nodes} nodes"
    if total_tasks is not None:
        return f"{total_tasks} tasks"
    return "unknown resources"


def get_total_tasks(metadata: dict[str, str]) -> int | None:
    total_tasks = metadata.get("total tasks")
    if total_tasks is None:
        return None

    try:
        return int(total_tasks)
    except ValueError:
        return None


def get_plot_name() -> str:
    if PLOT_NAME is not None:
        return PLOT_NAME

    return f"{SCRIPT_STEM}_{'_'.join(JOB_NAMES)}"


def plot_long_scaling(
    ax,
    summaries: dict[str, dict[str, object]],
) -> None:
    if USE_LOG_SCALE:
        ax.set_xscale("log")
        ax.set_yscale("log")

    visible_y_values = []
    for job_name in JOB_NAMES:
        summary = summaries[job_name]
        rows = get_plot_rows(summary["rows"])
        if not rows:
            continue

        images = np.asarray([row["images"] for row in rows])
        mlups = np.asarray([row["mlups"] for row in rows])
        if X_LIMIT is not None:
            visible_y_values.extend(mlups[images <= X_LIMIT])
        ax.plot(
            images,
            mlups,
            marker="o",
            markersize=MEASUREMENT_MARKER_SIZE,
            linewidth=MEASUREMENT_LINE_WIDTH,
            linestyle=MEASUREMENT_LINE_STYLE,
            color=JOB_COLORS[job_name],
            label=get_job_label(summary["metadata"]),
        )

    # capture limits from measured curves so ideal line does not expand them
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # calculate and plot ideal scaling line ('linear scaling') from 4-image baseline
    all_baseline_mlups = [
        row["mlups"]
        for summary in summaries.values()
        for row in get_plot_rows(summary["rows"])
        if row["images"] == IDEAL_BASELINE_IMAGES
    ]
    if all_baseline_mlups:
        baseline_mlups = np.mean(all_baseline_mlups)
        all_images = [
            row["images"]
            for summary in summaries.values()
            for row in get_plot_rows(summary["rows"])
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

    for job_name in JOB_NAMES:
        total_tasks = get_total_tasks(summaries[job_name]["metadata"])
        if total_tasks is None:
            print(f"warning: missing total task count for long_{job_name}")
            continue

        color = JOB_COLORS[job_name]
        ax.axvline(
            total_tasks,
            color=color,
            linestyle=TASK_LINE_STYLE,
            linewidth=TASK_LINE_WIDTH,
            alpha=TASK_LINE_ALPHA,
            zorder=0,
        )

    if X_LIMIT is not None:
        if USE_LOG_SCALE:
            ax.set_xlim(right=X_LIMIT)
        else:
            ax.set_xlim(left=0, right=X_LIMIT)

        if visible_y_values:
            visible_y_min = float(np.min(visible_y_values))
            visible_y_max = float(np.max(visible_y_values))
            if USE_LOG_SCALE:
                ax.set_ylim(
                    bottom=visible_y_min / (1.0 + Y_LIMIT_MARGIN),
                    top=visible_y_max * (1.0 + Y_LIMIT_MARGIN),
                )
            else:
                y_padding = Y_LIMIT_MARGIN * (visible_y_max - visible_y_min)
                ax.set_ylim(
                    bottom=max(0.0, visible_y_min - y_padding),
                    top=visible_y_max + y_padding,
                )

    ax.tick_params(axis="x", labelsize=X_AXIS_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=Y_AXIS_FONT_SIZE)
    if not USE_LOG_SCALE:
        tick_step = X_LIMIT_TICK_STEP if X_LIMIT is not None else X_AXIS_TICK_STEP
        ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.grid(alpha=GRID_ALPHA)
    ax.legend(
        loc=LEGEND_LOCATION,
        ncol=1,
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
    )


def plot_selected_jobs(summary_paths: dict[str, Path]) -> None:
    summaries = {
        job_name: load_summary(summary_paths[job_name])
        for job_name in JOB_NAMES
    }
    if not any(get_plot_rows(summary["rows"]) for summary in summaries.values()):
        raise ValueError("no measurements found in selected summaries")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plot_long_scaling(ax, summaries)

    for job_name in JOB_NAMES:
        partial_rows = [
            row for row in summaries[job_name]["rows"]
            if row["status"] != "complete"
        ]
        if partial_rows:
            partial_images = ", ".join(str(row["images"]) for row in partial_rows)
            print(f"warning: long_{job_name} includes incomplete measurements at {partial_images} images")

    output_path = PLOT_DIR / f"{get_plot_name()}.png"
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    summary_paths = {
        job_name: resolve_summary_path(job_name, JOB_IDS[job_name])
        for job_name in JOB_NAMES
    }
    print(f"plot directory: {PLOT_DIR}")
    for job_name in JOB_NAMES:
        print(f"summary {job_name}:  {summary_paths[job_name]}")
    print(f"jobs:           {', '.join(f'long_{job_name}' for job_name in JOB_NAMES)}")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_selected_jobs(summary_paths)
