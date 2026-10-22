#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np



# --- [ plot decomposition shape benchmark result as bar graph ] ---

# summary config
JOB_NAME = "C"  # options "C", "C2", "D", or "D2"
JOB_ID = None  # if "None", uses latest summary
INCLUDE_PARTIAL_MEASUREMENTS = True

# plot config
FIGURE_SIZE = (18, 6)
FIGURE_DPI = 600
GRID_ALPHA = 0.5
BAR_WIDTH = 0.65
BAR_EDGE_COLOR = "none"

HORIZONTAL_COLOR = np.array([30/255, 140/255, 80/255])
SQUARED_COLOR = np.array([50/255, 80/255, 170/255])
VERTICAL_COLOR = np.array([160/255, 110/255, 30/255])

# font size config
X_AXIS_FONT_SIZE = 16
Y_AXIS_FONT_SIZE = 16
TICK_LABEL_FONT_SIZE = 12

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


def load_summary(path: Path) -> list[dict[str, object]]:
    rows = []
    in_section = False

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped_line = line.strip()

        if stripped_line.startswith("## "):
            heading = stripped_line[3:].strip()
            in_section = (heading == "domain decompositions")
            continue

        if not in_section:
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

    return rows


def get_plot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if INCLUDE_PARTIAL_MEASUREMENTS:
        return rows

    return [row for row in rows if row["status"] == "complete"]


def get_plot_name(summary_path: Path) -> str:
    if PLOT_NAME is not None:
        return PLOT_NAME

    if JOB_ID is not None:
        return f"{SCRIPT_STEM}_{JOB_NAME}_{JOB_ID}"

    return f"{SCRIPT_STEM}_{JOB_NAME}"


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


def plot_decomposition_shape(
    ax,
    rows: list[dict[str, object]],
) -> None:
    decompositions = [row["decomposition"] for row in rows]
    mlups = [row["mlups"] for row in rows]
    x_pos = np.arange(len(decompositions))
    bar_colors = get_bar_colors(decompositions)

    ax.bar(
        x_pos,
        mlups,
        width=BAR_WIDTH,
        color=bar_colors,
        edgecolor=BAR_EDGE_COLOR,
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


def plot_selected_job() -> None:
    summary_path = resolve_summary_path(JOB_NAME, JOB_ID)

    all_rows = load_summary(summary_path)
    rows = get_plot_rows(all_rows)
    if not rows:
        raise ValueError(f"no measurements found in {summary_path}")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plot_decomposition_shape(ax, rows)

    partial_rows = [row for row in all_rows if row["status"] != "complete"]
    if partial_rows:
        partial_decomps = ", ".join(row["decomposition"] for row in partial_rows)
        print(f"warning: {JOB_NAME} includes incomplete measurements at {partial_decomps}")

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
