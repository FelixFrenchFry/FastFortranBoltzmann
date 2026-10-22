#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import numpy as np



# --- [ plot bulk decomposition scaling benchmark comparison ] ---

# summary config
BENCHMARK_VARIANTS = (5, 4)  # first uses filled circles/solid, second uses open triangles/dotted
BENCHMARK_LETTERS = ("A", "B")  # 6- and 8-node configurations
BENCHMARK_COLORS = {
    "A": (40/255, 80/255, 180/255),
    "B": (40/255, 180/255, 80/255),
    "C": (250/255, 0/255, 50/255),
    "D": (165/255, 40/255, 115/255),
    "E": (255/255, 180/255, 0/255),
    "F": (0/255, 160/255, 180/255),
    "G": (89/255, 89/255, 89/255),
}
BENCHMARK_MARKERS = ("o", "^")
BENCHMARK_MARKER_FACE_COLORS = ("auto", "none")
BENCHMARK_LINE_STYLES = ("-", ":")
BATCH_LEGEND_TITLES = {
    5: "batch 5 (per-host MPI)",
    4: "batch 4 (Slurm-managed MPI)",
}
INCLUDE_PARTIAL_MEASUREMENTS = True

# plot config
FIGURE_SIZE = (14, 6)
FIGURE_DPI = 600
GRID_ALPHA = 0.5
MEASUREMENT_MARKER_SIZE = 10.0
MEASUREMENT_LINE_WIDTH = 3.0
USE_LOG_SCALE = False
TASK_LINE_COLOR = "black"
TASK_LINE_STYLE = "-"
TASK_LINE_WIDTH = 3.0
TASK_LINE_ALPHA = 0.8
X_AXIS_TICK_STEP = 50
Y_AXIS_TICK_STEP = 1000
X_LIMIT = None  # None -> uses the full measured x range

IDEAL_LINE_COLOR = "black"
IDEAL_LINE_STYLE = ":"
IDEAL_LINE_WIDTH = 3.0
IDEAL_BASELINE_IMAGES = 2
LEGEND_LOCATION = "upper left"
LEGEND_SECONDARY_LOCATION = "lower right"
REFERENCE_LEGEND_LOCATION = "upper right"
LEGEND_COLUMNS = 1

# font size config
X_AXIS_FONT_SIZE = 16
Y_AXIS_FONT_SIZE = 16
LEGEND_FONT_SIZE = 12
LEGEND_TITLE_FONT_SIZE = 12

# path config
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = Path(__file__).resolve().parents[2]
SUMMARY_ROOT_DIR = ROOT_DIR / "scripts" / "cluster" / "bench_bulk"
SCRIPT_STEM = Path(__file__).stem
PLOT_DIR = SCRIPT_DIR / "plots"
PLOT_NAME = None  # None -> <SCRIPT_STEM>_<BENCHMARK_VARIANTS>


def validate_benchmark_variants() -> None:
    if len(BENCHMARK_VARIANTS) != 2:
        raise ValueError("BENCHMARK_VARIANTS must contain exactly two variants")
    if len(set(BENCHMARK_VARIANTS)) != 2:
        raise ValueError("BENCHMARK_VARIANTS must contain two different variants")
    if any(variant not in (1, 2, 3, 4, 5) for variant in BENCHMARK_VARIANTS):
        raise ValueError("BENCHMARK_VARIANTS must contain values between 1 and 5")


def get_benchmark_names(variant: int) -> tuple[str, ...]:
    suffix = "" if variant == 1 else str(variant)
    benchmark_names = tuple(
        f"{benchmark_letter}{suffix}"
        for benchmark_letter in BENCHMARK_LETTERS
    )
    return tuple(
        benchmark_name
        for benchmark_name in benchmark_names
        if list(
            (SUMMARY_ROOT_DIR / f"bench_{benchmark_name}").glob(
                "summary_*.md"
            )
        )
    )


def resolve_summary_path(
    benchmark_name: str,
    job_id: int | str | None = None,
) -> Path:
    benchmark_dir = SUMMARY_ROOT_DIR / f"bench_{benchmark_name}"
    if not benchmark_dir.exists():
        raise FileNotFoundError(f"missing benchmark directory: {benchmark_dir}")

    if job_id is not None:
        target_path = benchmark_dir / f"summary_{job_id}.md"
        if not target_path.exists():
            raise FileNotFoundError(
                f"missing summary for job id {job_id}: {target_path}"
            )
        return target_path

    summary_files = list(benchmark_dir.glob("summary_*.md"))
    if not summary_files:
        raise FileNotFoundError(
            f"no summary_*.md files found in {benchmark_dir}"
        )

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
            if len(cells) >= 2 and cells[0] in (
                "nodes",
                "tasks per node",
                "total tasks",
                "partition",
                "MPI per-host option",
            ):
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


def get_comparable_benchmark_letters(
    summaries: dict[int, dict[str, dict[str, object]]],
) -> tuple[str, ...]:
    comparable_letters = []
    for benchmark_letter in BENCHMARK_LETTERS:
        comparable = True
        for variant in BENCHMARK_VARIANTS:
            suffix = "" if variant == 1 else str(variant)
            benchmark_name = f"{benchmark_letter}{suffix}"
            summary = summaries[variant].get(benchmark_name)
            if summary is None or not get_plot_rows(summary["rows"]):
                comparable = False
                break
        if comparable:
            comparable_letters.append(benchmark_letter)

    return tuple(comparable_letters)


def get_total_tasks(metadata: dict[str, str]) -> int | None:
    total_tasks = metadata.get("total tasks")
    if total_tasks is None:
        return None

    try:
        return int(total_tasks)
    except ValueError:
        return None


def get_benchmark_label(metadata: dict[str, str]) -> str:
    nodes = metadata.get("nodes")
    total_tasks = get_total_tasks(metadata)
    if nodes and nodes != "-" and total_tasks is not None:
        return f"{nodes} nodes / {total_tasks} tasks"
    if nodes and nodes != "-":
        return f"{nodes} nodes"
    if total_tasks is not None:
        return f"{total_tasks} tasks"
    return "unknown resources"


def get_benchmark_color(benchmark_name: str):
    return BENCHMARK_COLORS[benchmark_name[0]]


def get_plot_name() -> str:
    if PLOT_NAME is not None:
        return PLOT_NAME

    variants = "_".join(str(variant) for variant in BENCHMARK_VARIANTS)
    return f"{SCRIPT_STEM}_{variants}"


def plot_bulk_comparison(
    ax,
    summaries: dict[int, dict[str, dict[str, object]]],
) -> None:
    if USE_LOG_SCALE:
        ax.set_xscale("log")
        ax.set_yscale("log")

    all_summaries = [
        summary
        for batch_summaries in summaries.values()
        for summary in batch_summaries.values()
    ]
    comparable_benchmark_letters = get_comparable_benchmark_letters(summaries)
    legend_handles = [[], []]
    legend_labels = [[], []]
    for batch_index, variant in enumerate(BENCHMARK_VARIANTS):
        for benchmark_name, summary in summaries[variant].items():
            if benchmark_name[0] not in comparable_benchmark_letters:
                continue

            rows = get_plot_rows(summary["rows"])
            if not rows:
                continue

            images = np.asarray([row["images"] for row in rows])
            mlups = np.asarray([row["mlups"] for row in rows])
            line, = ax.plot(
                images,
                mlups,
                marker=BENCHMARK_MARKERS[batch_index],
                markersize=MEASUREMENT_MARKER_SIZE,
                linewidth=MEASUREMENT_LINE_WIDTH,
                linestyle=BENCHMARK_LINE_STYLES[batch_index],
                color=get_benchmark_color(benchmark_name),
                markerfacecolor=BENCHMARK_MARKER_FACE_COLORS[batch_index],
                label="_nolegend_",
            )
            legend_handles[batch_index].append(line)
            legend_labels[batch_index].append(
                get_benchmark_label(summary["metadata"])
            )

    # keep axis limits based on every available configuration, including hidden curves
    all_plot_rows = [
        row
        for summary in all_summaries
        for row in get_plot_rows(summary["rows"])
    ]
    if all_plot_rows:
        ax.update_datalim(
            np.column_stack(
                (
                    [row["images"] for row in all_plot_rows],
                    [row["mlups"] for row in all_plot_rows],
                )
            )
        )
        ax.autoscale_view()

    # capture limits from measured curves so ideal line does not expand them
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # calculate and plot ideal scaling line ('linear scaling') from baseline image count
    all_baseline_mlups = [
        row["mlups"]
        for summary in all_summaries
        for row in get_plot_rows(summary["rows"])
        if row["images"] == IDEAL_BASELINE_IMAGES
    ]
    ideal_line = None
    if all_baseline_mlups:
        baseline_mlups = np.mean(all_baseline_mlups)
        all_images = [
            row["images"]
            for summary in all_summaries
            for row in get_plot_rows(summary["rows"])
        ]
        if all_images:
            x_min, x_max = min(all_images), max(all_images)
            x_linear = np.linspace(x_min, x_max, 200)
            y_linear = baseline_mlups * (x_linear / float(IDEAL_BASELINE_IMAGES))

            ideal_line, = ax.plot(
                x_linear,
                y_linear,
                color=IDEAL_LINE_COLOR,
                linestyle=IDEAL_LINE_STYLE,
                linewidth=IDEAL_LINE_WIDTH,
                label="_nolegend_",
                zorder=1,
            )

    # restore limits so linear line may run off-screen
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    task_counts = sorted(
        {
            total_tasks
            for summary in all_summaries
            for total_tasks in [get_total_tasks(summary["metadata"])]
            if total_tasks is not None
        }
    )
    for total_tasks in task_counts:
        ax.axvline(
            total_tasks,
            color=TASK_LINE_COLOR,
            linestyle=TASK_LINE_STYLE,
            linewidth=TASK_LINE_WIDTH,
            alpha=TASK_LINE_ALPHA,
            zorder=0,
        )

    if X_LIMIT is not None:
        ax.set_xlim(right=X_LIMIT)

    ax.tick_params(axis="x", labelsize=X_AXIS_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=Y_AXIS_FONT_SIZE)
    if not USE_LOG_SCALE:
        ax.xaxis.set_major_locator(MultipleLocator(X_AXIS_TICK_STEP))
        ax.yaxis.set_major_locator(MultipleLocator(Y_AXIS_TICK_STEP))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.grid(alpha=GRID_ALPHA)
    for batch_index, location in enumerate(
        (LEGEND_LOCATION, LEGEND_SECONDARY_LOCATION)
    ):
        if not legend_handles[batch_index]:
            continue
        batch_legend = ax.legend(
            handles=legend_handles[batch_index],
            labels=legend_labels[batch_index],
            loc=location,
            ncol=LEGEND_COLUMNS,
            frameon=True,
            fontsize=LEGEND_FONT_SIZE,
            title=BATCH_LEGEND_TITLES.get(
                BENCHMARK_VARIANTS[batch_index],
                f"batch {BENCHMARK_VARIANTS[batch_index]}",
            ),
            title_fontsize=LEGEND_TITLE_FONT_SIZE,
        )
        ax.add_artist(batch_legend)

    if ideal_line is not None:
        ax.legend(
            handles=[ideal_line],
            labels=["linear scaling"],
            loc=REFERENCE_LEGEND_LOCATION,
            ncol=1,
            frameon=True,
            fontsize=LEGEND_FONT_SIZE,
        )


def plot_selected_benchmarks(
    summary_paths: dict[int, dict[str, Path]],
) -> None:
    summaries = {
        variant: {
            benchmark_name: load_summary(summary_path)
            for benchmark_name, summary_path in variant_paths.items()
        }
        for variant, variant_paths in summary_paths.items()
    }
    all_summaries = [
        summary
        for variant_summaries in summaries.values()
        for summary in variant_summaries.values()
    ]
    if not any(get_plot_rows(summary["rows"]) for summary in all_summaries):
        raise ValueError("no measurements found in selected benchmark summaries")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plot_bulk_comparison(ax, summaries)

    for variant, variant_summaries in summaries.items():
        for benchmark_name, summary in variant_summaries.items():
            partial_rows = [
                row for row in summary["rows"]
                if row["status"] != "complete"
            ]
            if partial_rows:
                partial_images = ", ".join(str(row["images"]) for row in partial_rows)
                print(
                    f"warning: bench_{benchmark_name} (batch {variant}) includes "
                    f"incomplete measurements at {partial_images} images"
                )

    output_path = PLOT_DIR / f"{get_plot_name()}.png"
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    print(f"saved plot: {output_path}")


if __name__ == "__main__":
    validate_benchmark_variants()

    enabled_benchmarks = {
        variant: get_benchmark_names(variant)
        for variant in BENCHMARK_VARIANTS
    }
    for variant, benchmark_names in enabled_benchmarks.items():
        if not benchmark_names:
            raise ValueError(
                f"no summaries found for benchmark variant {variant}"
            )

    summary_paths = {
        variant: {
            benchmark_name: resolve_summary_path(benchmark_name)
            for benchmark_name in benchmark_names
        }
        for variant, benchmark_names in enabled_benchmarks.items()
    }
    print(f"plot directory: {PLOT_DIR}")
    print(
        "benchmark variants: "
        + ", ".join(str(variant) for variant in BENCHMARK_VARIANTS)
    )
    for variant, benchmark_names in enabled_benchmarks.items():
        print(
            f"batch {variant}:        "
            + ", ".join(f"bench_{name}" for name in benchmark_names)
        )
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_selected_benchmarks(summary_paths)
