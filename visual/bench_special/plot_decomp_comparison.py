#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np



# --- [ plot decomposition comparison benchmark results as bar graph ] ---

# summary config
JOB_NAMES = ("D2", "D2")  # options ("C", "D") or ("C2", "D2")
JOB_IDS = (5972307, 5948625)  # if "None", uses latest summary for each job
JOB_LABELS = None  # None -> summary node/task metadata
INCLUDE_PARTIAL_MEASUREMENTS = True

# plot config
FIGURE_SIZE = (18, 6)
FIGURE_DPI = 600
GRID_ALPHA = 0.5
BAR_WIDTH = 0.38
BAR_EDGE_COLOR = "none"
HORIZONTAL_COLOR = np.array([30/255, 140/255, 80/255])
SQUARED_COLOR = np.array([50/255, 80/255, 170/255])
VERTICAL_COLOR = np.array([160/255, 110/255, 30/255])
BAR_HATCHES = ("", "//")
HATCH_EDGE_COLOR = "black"
PLACEHOLDER_COLOR = (110/255, 110/255, 110/255)
PLACEHOLDER_HATCH = "//"
PLACEHOLDER_HEIGHT_FRACTION = 0.03
PLACEHOLDER_TEXT = "n/a"

# font size config
X_AXIS_FONT_SIZE = 16
Y_AXIS_FONT_SIZE = 16
TICK_LABEL_FONT_SIZE = 12
LEGEND_FONT_SIZE = 12
LEGEND_LOCATION = "upper right"
LEGEND_COLUMNS = 1

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[2]
SUMMARY_DIRS = {
    "C": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_C",
    "C2": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_C2",
    "D": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_D",
    "D2": ROOT_DIR / "scripts" / "cluster" / "bench_special" / "bench_decomp_D2",
}
SCRIPT_STEM = Path(__file__).stem
PLOT_DIR = SCRIPT_DIR / "plots"
PLOT_NAME = None  # None -> derived from selected jobs and run ids


def validate_job_names() -> None:
    if len(JOB_NAMES) != 2:
        raise ValueError("JOB_NAMES must contain exactly two jobs")
    if len(JOB_IDS) != 2:
        raise ValueError("JOB_IDS must contain exactly two job ids")
    if (JOB_NAMES[0], str(JOB_IDS[0])) == (JOB_NAMES[1], str(JOB_IDS[1])):
        raise ValueError("JOB_NAMES and JOB_IDS must select two different runs")

    unsupported_jobs = [job_name for job_name in JOB_NAMES if job_name not in SUMMARY_DIRS]
    if unsupported_jobs:
        raise ValueError(f"unsupported JOB_NAMES: {unsupported_jobs}")


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
            heading = stripped_line[3:].strip()
            in_section = (heading == "domain decompositions")
            continue

        if not in_section:
            if stripped_line.startswith("|") and stripped_line.endswith("|"):
                cells = split_table_row(stripped_line)
                if len(cells) >= 2 and cells[0] in (
                    "nodes",
                    "tasks per node",
                ):
                    metadata[cells[0]] = cells[1]
            continue

        if not stripped_line.startswith("|") or not stripped_line.endswith("|"):
            continue

        cells = split_table_row(stripped_line)
        if len(cells) < 9 or cells[0] == "images":
            continue

        try:
            images = int(cells[0])
            decomposition = cells[1]
            mlups = int(cells[5])
        except ValueError:
            continue

        rows.append(
            {
                "images": images,
                "decomposition": decomposition,
                "mlups": mlups,
                "status": cells[8],
            }
        )

    return {
        "metadata": metadata,
        "rows": rows,
    }


def get_plot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if INCLUDE_PARTIAL_MEASUREMENTS:
        return rows

    return [row for row in rows if row["status"] == "complete"]


def get_job_run_suffix(job_index: int) -> str:
    if JOB_NAMES[0] != JOB_NAMES[1]:
        return ""

    job_id = JOB_IDS[job_index]
    if job_id is None:
        return " (latest)"
    return f" (job {job_id})"


def get_job_label(job_index: int) -> str:
    if JOB_LABELS is not None:
        return JOB_LABELS[job_index]

    return f"bench_decomp_{JOB_NAMES[job_index]}{get_job_run_suffix(job_index)}"


def get_job_legend_label(
    job_index: int,
    metadata: dict[str, str],
) -> str:
    if JOB_LABELS is not None:
        return JOB_LABELS[job_index]

    nodes = metadata.get("nodes")
    tasks_per_node = metadata.get("tasks per node")
    if nodes and nodes != "-" and tasks_per_node and tasks_per_node != "-":
        label = f"{nodes} nodes / {tasks_per_node} tasks per node"
    elif nodes and nodes != "-":
        label = f"{nodes} nodes"
    elif tasks_per_node and tasks_per_node != "-":
        label = f"{tasks_per_node} tasks per node"
    else:
        return get_job_label(job_index)

    return f"{label}{get_job_run_suffix(job_index)}"


def get_plot_name() -> str:
    if PLOT_NAME is not None:
        return PLOT_NAME

    if JOB_NAMES[0] == JOB_NAMES[1]:
        job_ids = tuple(
            "latest" if job_id is None else str(job_id)
            for job_id in JOB_IDS
        )
        return f"{SCRIPT_STEM}_{JOB_NAMES[0]}_{'_'.join(job_ids)}"

    return f"{SCRIPT_STEM}_{'_'.join(JOB_NAMES)}"


def get_decompositions(
    rows_by_job: tuple[list[dict[str, object]], ...],
) -> list[str]:
    decompositions = []
    for rows in rows_by_job:
        for row in rows:
            decomposition = row["decomposition"]
            if decomposition not in decompositions:
                decompositions.append(decomposition)
    return decompositions


def get_bar_colors(decompositions: list[str]) -> list[tuple[float, float, float]]:
    parsed = []
    for decomp in decompositions:
        p_str, q_str = decomp.split("x")
        parsed.append((int(p_str.strip()), int(q_str.strip())))

    mid_idx = min(range(len(parsed)), key=lambda i: abs(parsed[i][0] - parsed[i][1]))
    n_bars = len(decompositions)

    colors = []
    for i in range(n_bars):
        if i <= mid_idx:
            t = i / max(1, mid_idx)
            color = (1.0 - t) * HORIZONTAL_COLOR + t * SQUARED_COLOR
        else:
            t = (i - mid_idx) / max(1, n_bars - 1 - mid_idx)
            color = (1.0 - t) * SQUARED_COLOR + t * VERTICAL_COLOR
        colors.append(tuple(color))

    return colors


def plot_decomposition_comparison(
    ax,
    rows_by_job: tuple[list[dict[str, object]], ...],
    metadata_by_job: tuple[dict[str, str], ...],
) -> list[str]:
    decompositions = get_decompositions(rows_by_job)
    row_maps = [
        {row["decomposition"]: row for row in rows}
        for rows in rows_by_job
    ]
    x_pos = np.arange(len(decompositions))
    job_offsets = (
        -0.5 * BAR_WIDTH,
        0.5 * BAR_WIDTH,
    )
    measured_values = [
        row["mlups"]
        for row_map in row_maps
        for row in row_map.values()
    ]
    if not measured_values:
        raise ValueError("no measurements found in selected summaries")

    bar_colors = get_bar_colors(decompositions)
    placeholder_height = max(measured_values) * PLACEHOLDER_HEIGHT_FRACTION
    for job_index, row_map in enumerate(row_maps):
        measured_x = []
        measured_mlups = []
        measured_colors = []
        missing_x = []

        for index, decomposition in enumerate(decompositions):
            x_value = x_pos[index] + job_offsets[job_index]
            row = row_map.get(decomposition)
            if row is None:
                missing_x.append(x_value)
            else:
                measured_x.append(x_value)
                measured_mlups.append(row["mlups"])
                measured_colors.append(bar_colors[index])

        bar_hatch = BAR_HATCHES[job_index]
        bar_edge_color = HATCH_EDGE_COLOR if bar_hatch else BAR_EDGE_COLOR

        if measured_x:
            ax.bar(
                measured_x,
                measured_mlups,
                width=BAR_WIDTH,
                color=measured_colors,
                edgecolor=bar_edge_color,
                hatch=bar_hatch,
                label=get_job_legend_label(job_index, metadata_by_job[job_index]),
            )
        else:
            ax.bar(
                [],
                [],
                width=BAR_WIDTH,
                color=bar_colors[0],
                edgecolor=bar_edge_color,
                hatch=bar_hatch,
                label=get_job_legend_label(job_index, metadata_by_job[job_index]),
            )

        if missing_x:
            ax.bar(
                missing_x,
                [placeholder_height] * len(missing_x),
                width=BAR_WIDTH,
                color="none",
                edgecolor=PLACEHOLDER_COLOR,
                hatch=PLACEHOLDER_HATCH,
                label="_nolegend_",
            )
            for x_value in missing_x:
                ax.text(
                    x_value,
                    placeholder_height,
                    PLACEHOLDER_TEXT,
                    ha="center",
                    va="bottom",
                    fontsize=TICK_LABEL_FONT_SIZE,
                    color=PLACEHOLDER_COLOR,
                )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        decompositions,
        rotation=45,
        ha="right",
        fontsize=TICK_LABEL_FONT_SIZE,
    )
    ax.tick_params(axis="y", labelsize=Y_AXIS_FONT_SIZE)
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.grid(alpha=GRID_ALPHA, axis="y")
    ax.legend(
        loc=LEGEND_LOCATION,
        ncol=LEGEND_COLUMNS,
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
    )

    return decompositions


def plot_selected_jobs() -> None:
    summary_paths = tuple(
        resolve_summary_path(job_name, job_id)
        for job_name, job_id in zip(JOB_NAMES, JOB_IDS)
    )
    all_summaries = tuple(load_summary(path) for path in summary_paths)
    all_rows_by_job = tuple(summary["rows"] for summary in all_summaries)
    rows_by_job = tuple(get_plot_rows(rows) for rows in all_rows_by_job)
    if not any(rows_by_job):
        raise ValueError("no measurements found in selected summaries")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    metadata_by_job = tuple(summary["metadata"] for summary in all_summaries)
    decompositions = plot_decomposition_comparison(
        ax,
        rows_by_job,
        metadata_by_job,
    )

    for job_index, (all_rows, rows) in enumerate(zip(all_rows_by_job, rows_by_job)):
        row_decompositions = {row["decomposition"] for row in rows}
        missing_decompositions = [
            decomposition
            for decomposition in decompositions
            if decomposition not in row_decompositions
        ]
        if missing_decompositions:
            print(
                f"warning: {get_job_label(job_index)} is missing values at "
                + ", ".join(missing_decompositions)
            )

        partial_rows = [row for row in all_rows if row["status"] != "complete"]
        if partial_rows:
            partial_decompositions = ", ".join(
                row["decomposition"] for row in partial_rows
            )
            print(
                f"warning: {get_job_label(job_index)} includes incomplete "
                f"measurements at {partial_decompositions}"
            )

    output_path = PLOT_DIR / f"{get_plot_name()}.png"
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    validate_job_names()

    summary_paths = tuple(
        resolve_summary_path(job_name, job_id)
        for job_name, job_id in zip(JOB_NAMES, JOB_IDS)
    )
    print(f"plot directory: {PLOT_DIR}")
    for job_index, (job_name, summary_path) in enumerate(zip(JOB_NAMES, summary_paths)):
        job_id = JOB_IDS[job_index] or summary_path.stem.replace("summary_", "")
        print(f"summary {job_index + 1}:    {summary_path}")
        print(f"job {job_index + 1}:        bench_decomp_{job_name} (job id: {job_id})")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_selected_jobs()
