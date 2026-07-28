################
#### CONFIG ####
################

# ---- DATA PATH ----
from pathlib import Path

DATA_PATH_PREFIX = "../../database/exports/"

DATA_PATH_SUFFIX = "barcode01_agilent/metrics_view_cluster/at_decoding.csv"

DATA_PATH = Path(DATA_PATH_PREFIX) / Path(DATA_PATH_SUFFIX)

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}


# ---- DATA SPECS ----

EXPECTED_COLUMN_NAMES_IN_THE_CSV = [
    "label",
    "dec_run_id",
    "item_id",
    "cluster_id",
    "count_used_for_consensus",
]

LABELS_NAME_MAP = {
    "barcode01_agilent_alignment_decoding": "Read pool from supplier 1 (Agilent)"
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/cluster_size_distribution_at_decoding.png"
OUTPUT_PATH_SVG = "../plots/cluster_size_distribution_at_decoding.svg"

# ---- NAMES ----

FIGURE_TITLE = "Cluster size distributions for items JPEG DNA and JPEG DNA delta G over multiple runs"

FIGURE_DESCRIPTION = f"Average cluster size distribution at decoding over all runs of label {'barcode01_agilent_alignment_decoding'}."

PLOT_NAME_MAP = {
    "cluster_size_distribution": "Cluster size distribution (average over the runs with standard deviation)",
}

X_AXIS_NAME_MAP = {
    "cluster_size": "Coverage (size of the cluster)",
}

Y_AXIS_NAME_MAP = {
    "cluster_size_distribution": "Number of clusters having the given coverage",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green (fallback)
}

PLOT_STYLE_MAP = {
    "cluster_size_distribution": {
        "alpha": 0.85,
        "bar_width": 0.42,
        "edgecolor": "black",
        "linewidth": 0.5,
        "error_capsize": 3,
        "error_linewidth": 1.1,
    },
    "average_coverage_marker": {
        "linestyle": "--",
        "linewidth": 2.0,
        "alpha": 0.9,
    }
}


##############
#### CODE ####
##############

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def _validate_columns(df: pd.DataFrame, expected_columns: list[str]) -> None:
    missing_columns = [column for column in expected_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing expected columns in CSV: {missing_columns}")


def _compute_distribution_stats(
    df: pd.DataFrame,
) -> tuple[list[int], dict[int, np.ndarray], dict[int, np.ndarray], dict[int, float]]:
    df = df.copy()
    df["count_used_for_consensus"] = (
        pd.to_numeric(df["count_used_for_consensus"], errors="coerce").fillna(0).astype(int)
    )

    run_item_cluster_counts = (
        df.groupby(["dec_run_id", "item_id", "count_used_for_consensus"])
        .size()
        .rename("n_refs")
        .reset_index()
    )

    run_item_histograms = run_item_cluster_counts.pivot_table(
        index=["dec_run_id", "item_id"],
        columns="count_used_for_consensus",
        values="n_refs",
        fill_value=0,
        aggfunc="sum",
    )

    expected_totals = df.groupby(["dec_run_id", "item_id"]).size().sort_index()
    actual_totals = run_item_histograms.sum(axis=1).sort_index()
    if not actual_totals.equals(expected_totals):
        raise ValueError("Inconsistent histogram totals: sum over cluster sizes differs from reference counts.")

    cluster_sizes = sorted(int(size) for size in run_item_histograms.columns.tolist())

    mean_by_item: dict[int, np.ndarray] = {}
    std_by_item: dict[int, np.ndarray] = {}
    avg_coverage_by_item: dict[int, float] = {}

    for item_id in sorted(df["item_id"].unique().tolist()):
        if item_id not in ITEM_ID_NAME_MAP:
            continue
        item_histograms = run_item_histograms.xs(item_id, level="item_id")
        if isinstance(item_histograms, pd.Series):
            item_histograms = item_histograms.to_frame().T
        item_histograms = item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)
        mean_by_item[item_id] = item_histograms.mean(axis=0).to_numpy(dtype=float)
        std_by_item[item_id] = item_histograms.std(axis=0, ddof=0).to_numpy(dtype=float)

    # Per-item average coverage from cluster view: run mean over references, then mean over runs.
    run_level_mean_coverage = (
        df.groupby(["dec_run_id", "item_id"], as_index=False)["count_used_for_consensus"].mean()
    )
    for item_id in sorted(df["item_id"].unique().tolist()):
        if item_id not in ITEM_ID_NAME_MAP:
            continue
        item_run_means = run_level_mean_coverage[
            run_level_mean_coverage["item_id"] == item_id
        ]["count_used_for_consensus"]
        avg_coverage_by_item[item_id] = float(item_run_means.mean())

    return cluster_sizes, mean_by_item, std_by_item, avg_coverage_by_item


def main() -> None:
    data_path = _resolve_path_from_script(DATA_PATH)
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_svg = _resolve_path_from_script(Path(OUTPUT_PATH_SVG))

    df = pd.read_csv(data_path)
    _validate_columns(df, EXPECTED_COLUMN_NAMES_IN_THE_CSV)

    selected_labels = set(LABELS_NAME_MAP.keys())
    df = df[df["label"].isin(selected_labels)]
    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())]

    if df.empty:
        raise ValueError("No data available after filtering by label and item_id.")

    cluster_sizes, mean_by_item, std_by_item, avg_coverage_by_item = _compute_distribution_stats(df)

    style = PLOT_STYLE_MAP["cluster_size_distribution"]
    fig, ax = plt.subplots(figsize=(13, 7))

    x = np.array(cluster_sizes, dtype=float)
    item_ids = sorted(mean_by_item.keys())
    n_items = len(item_ids)
    bar_width = float(style["bar_width"])
    marker_style = PLOT_STYLE_MAP["average_coverage_marker"]

    for idx, item_id in enumerate(item_ids):
        offsets = x + (idx - (n_items - 1) / 2.0) * bar_width
        ax.bar(
            offsets,
            mean_by_item[item_id],
            width=bar_width,
            yerr=std_by_item[item_id],
            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
            alpha=style["alpha"],
            edgecolor=style["edgecolor"],
            linewidth=style["linewidth"],
            capsize=style["error_capsize"],
            error_kw={"elinewidth": style["error_linewidth"]},
            label=ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}"),
        )
        if item_id in avg_coverage_by_item:
            avg_cov = avg_coverage_by_item[item_id]
            ax.axvline(
                avg_cov,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                linestyle=marker_style["linestyle"],
                linewidth=marker_style["linewidth"],
                alpha=marker_style["alpha"],
                label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} mean cov={avg_cov:.2f}",
            )

    ax.set_title(PLOT_NAME_MAP["cluster_size_distribution"])
    ax.set_xlabel(X_AXIS_NAME_MAP["cluster_size"])
    ax.set_ylabel(Y_AXIS_NAME_MAP["cluster_size_distribution"])
    ax.set_xticks(x)
    ax.set_xticklabels([str(size) for size in cluster_sizes])
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    ax.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.94, FIGURE_DESCRIPTION, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.05, 0.97, 0.88))

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=200)
    fig.savefig(output_path_svg)
    plt.close(fig)


if __name__ == "__main__":
    main()
