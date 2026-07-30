################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

DB_PATH = "../../database/dspl.db"

RUN_LABEL = "barcode01_agilent_alignment_decoding"

ITEM_IDS_TO_PLOT = [0, 1]

MIN_PERFECT_DECODED_RATIO = 0.5

VISUAL_INFINITY_FACTOR = 1.08

SQL_QUERY = """
WITH eligible_run_item AS (
    SELECT
            m.exp_id,
            m.dec_run_id,
            m.item_id
    FROM metrics_view_run_cluster_at_decoding m
    JOIN decoding_run_label_record l
        ON l.exp_id = m.exp_id
    AND l.dec_run_id = m.dec_run_id
    WHERE l.label = ?
        AND m.item_id IN ({item_placeholders})
    GROUP BY m.exp_id, m.dec_run_id, m.item_id
    HAVING SUM(
        CASE
            WHEN NULLIF(TRIM(CAST(m.dec_pass_id AS TEXT)), '') IS NOT NULL THEN 1
            ELSE 0
        END
    ) > 0
)
SELECT
        m.dec_run_id,
        m.item_id,
    m.region_id,
    m.position_id,
    m.count_used_for_consensus,
    m.count_at_first_decoding,
    m.hamming_dist_to_ref
FROM metrics_view_run_cluster_at_decoding m
JOIN eligible_run_item e
    ON e.exp_id = m.exp_id
   AND e.dec_run_id = m.dec_run_id
   AND e.item_id = m.item_id
ORDER BY m.dec_run_id ASC, m.item_id ASC
"""

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/cluster_size_distribution_at_decoding.png"
OUTPUT_PATH_SVG = "../plots/cluster_size_distribution_at_decoding.svg"

# ---- NAMES ----

FIGURE_TITLE = "Cluster size distributions for items JPEG DNA and JPEG DNA delta G over multiple runs"

FIGURE_DESCRIPTION = (
    f"Average normalized cluster size distribution at decoding over all runs of label {RUN_LABEL}. "
    "For each item, $\\mu_{\\text{cov}}$ and $\\sigma_{\\text{cov}}$ are computed from the aggregated "
    "inter-run relative-frequency distribution over coverage bins."
)

PLOT_NAME_MAP = {
    "cluster_size_distribution": "Cluster size distribution (inter-run mean relative frequency with standard error)",
    "cluster_size_at_first_decoding": "Cluster size at first decoding (zeros shown at infinity)",
}

X_AXIS_NAME_MAP = {
    "cluster_size": "Coverage (size of the cluster)",
    "cluster_size_at_first_decoding": "Cluster size at first decoding (integer bins)",
}

Y_AXIS_NAME_MAP = {
    "cluster_size_distribution": "Relative frequency at image decoding time",
    "cluster_size_at_first_decoding": "Number of references",
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

import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def _compute_distribution_stats(
    df: pd.DataFrame,
) -> tuple[list[int], dict[int, np.ndarray], dict[int, np.ndarray], dict[int, tuple[float, float]]]:
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

    # Normalize per (run, item) so each row is a proper distribution summing to 1.
    run_item_histograms = run_item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)
    run_item_relative_histograms = run_item_histograms.div(
        run_item_histograms.sum(axis=1),
        axis=0,
    )

    mean_by_item: dict[int, np.ndarray] = {}
    se_by_item: dict[int, np.ndarray] = {}
    coverage_moments_by_item: dict[int, tuple[float, float]] = {}

    for item_id in sorted(df["item_id"].unique().tolist()):
        if item_id not in ITEM_ID_NAME_MAP:
            continue
        item_histograms = run_item_relative_histograms.xs(item_id, level="item_id")
        if isinstance(item_histograms, pd.Series):
            item_histograms = item_histograms.to_frame().T
        item_histograms = item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)
        mean_distribution = item_histograms.mean(axis=0).to_numpy(dtype=float)
        n_runs = item_histograms.shape[0]
        if n_runs > 1:
            std_distribution = item_histograms.std(axis=0, ddof=1).to_numpy(dtype=float)
            se_distribution = std_distribution / np.sqrt(float(n_runs))
        else:
            se_distribution = np.zeros_like(mean_distribution)

        mean_by_item[item_id] = mean_distribution
        se_by_item[item_id] = se_distribution

        support = np.array(cluster_sizes, dtype=float)
        total_prob = float(mean_distribution.sum())
        if total_prob <= 0:
            coverage_moments_by_item[item_id] = (0.0, 0.0)
        else:
            normalized_distribution = mean_distribution / total_prob
            mu_cov = float(np.sum(support * normalized_distribution))
            var_cov = float(np.sum(((support - mu_cov) ** 2) * normalized_distribution))
            sigma_cov = float(np.sqrt(max(var_cov, 0.0)))
            coverage_moments_by_item[item_id] = (mu_cov, sigma_cov)

    return cluster_sizes, mean_by_item, se_by_item, coverage_moments_by_item


def _compute_first_decoding_reference_stats(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df["count_at_first_decoding"] = pd.to_numeric(df["count_at_first_decoding"], errors="coerce")
    df["count_at_first_decoding"] = df["count_at_first_decoding"].fillna(0)

    reference_level_df = (
        df.groupby(["item_id", "region_id", "position_id"], as_index=False)
        .agg(
            mean_positive_count_at_first_decoding=(
                "count_at_first_decoding",
                lambda s: pd.to_numeric(s, errors="coerce").loc[
                    pd.to_numeric(s, errors="coerce") > 0
                ].mean(),
            ),
            n_runs_total=("count_at_first_decoding", "size"),
            n_runs_zero=(
                "count_at_first_decoding",
                lambda s: int((pd.to_numeric(s, errors="coerce") <= 0).sum()),
            ),
        )
    )
    reference_level_df["zero_ratio"] = reference_level_df["n_runs_zero"] / reference_level_df["n_runs_total"]

    included_df = reference_level_df[reference_level_df["zero_ratio"] < MIN_PERFECT_DECODED_RATIO].copy()
    excluded_df = reference_level_df[reference_level_df["zero_ratio"] >= MIN_PERFECT_DECODED_RATIO].copy()

    return reference_level_df, included_df, excluded_df


def main() -> None:
    db_path = _resolve_path_from_script(Path(DB_PATH))
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_svg = _resolve_path_from_script(Path(OUTPUT_PATH_SVG))

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    query = SQL_QUERY.format(item_placeholders=item_placeholders)
    query_params = [RUN_LABEL, *ITEM_IDS_TO_PLOT]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())]

    if df.empty:
        raise ValueError("No data returned by SQL query. Check RUN_LABEL and ITEM_IDS_TO_PLOT.")

    cluster_sizes, mean_by_item, se_by_item, coverage_moments_by_item = _compute_distribution_stats(df)
    reference_level_df, included_reference_df, excluded_reference_df = _compute_first_decoding_reference_stats(df)

    style = PLOT_STYLE_MAP["cluster_size_distribution"]
    fig, (ax_item_distribution, ax_first_decoding) = plt.subplots(
        2,
        1,
        figsize=(13, 13),
        sharex=False,
        gridspec_kw={"height_ratios": [1.0, 1.0]},
    )

    x = np.array(cluster_sizes, dtype=float)
    item_ids = sorted(mean_by_item.keys())
    n_items = len(item_ids)
    bar_width = float(style["bar_width"])
    marker_style = PLOT_STYLE_MAP["average_coverage_marker"]

    for idx, item_id in enumerate(item_ids):
        offsets = x + (idx - (n_items - 1) / 2.0) * bar_width
        ax_item_distribution.bar(
            offsets,
            mean_by_item[item_id],
            width=bar_width,
            yerr=se_by_item[item_id],
            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
            alpha=style["alpha"],
            edgecolor=style["edgecolor"],
            linewidth=style["linewidth"],
            capsize=style["error_capsize"],
            error_kw={"elinewidth": style["error_linewidth"]},
            label=ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}"),
        )
        if item_id in coverage_moments_by_item:
            mu_cov, sigma_cov = coverage_moments_by_item[item_id]
            ax_item_distribution.axvline(
                mu_cov,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                linestyle=marker_style["linestyle"],
                linewidth=marker_style["linewidth"],
                alpha=marker_style["alpha"],
                label=(
                    f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} "
                    f"($\\mu_{{\\text{{cov}}}}={mu_cov:.2f}, \\sigma_{{\\text{{cov}}}}={sigma_cov:.2f}$)"
                ),
            )

    ax_item_distribution.set_title(PLOT_NAME_MAP["cluster_size_distribution"])
    ax_item_distribution.set_xlabel(X_AXIS_NAME_MAP["cluster_size"])
    ax_item_distribution.set_ylabel(Y_AXIS_NAME_MAP["cluster_size_distribution"])
    ax_item_distribution.set_xticks(x)
    ax_item_distribution.set_xticklabels([str(size) for size in cluster_sizes])
    ax_item_distribution.grid(axis="y", alpha=0.25, linestyle="--")
    ax_item_distribution.legend()

    finite_first_decoding = included_reference_df["mean_positive_count_at_first_decoding"].dropna()
    max_finite_first_decoding = (
        int(np.ceil(float(finite_first_decoding.max()))) if not finite_first_decoding.empty else 1
    )
    visual_infinity_x = float(max_finite_first_decoding + 2)

    first_decoding_rows: list[pd.DataFrame] = []

    for item_id in sorted(ITEM_ID_NAME_MAP.keys()):
        item_included = included_reference_df[included_reference_df["item_id"] == item_id].copy()
        item_excluded = excluded_reference_df[excluded_reference_df["item_id"] == item_id].copy()

        if not item_included.empty:
            item_included["bin_upper"] = (
                np.ceil(item_included["mean_positive_count_at_first_decoding"]).astype(int)
            )
            item_included = (
                item_included.groupby("bin_upper", as_index=False)
                .agg(n_points=("mean_positive_count_at_first_decoding", "size"))
            )
            item_included["item_id"] = item_id
            item_included["kind"] = "included"
            item_included["bin_center"] = item_included["bin_upper"].astype(float)
            first_decoding_rows.append(item_included)

        if not item_excluded.empty:
            first_decoding_rows.append(
                pd.DataFrame(
                    [
                        {
                            "bin_upper": visual_infinity_x,
                            "n_points": len(item_excluded),
                            "item_id": item_id,
                            "kind": "excluded",
                            "bin_center": visual_infinity_x,
                        }
                    ]
                )
            )

    if first_decoding_rows:
        first_decoding_df = pd.concat(first_decoding_rows, ignore_index=True)
        item_ids_in_bins = sorted(first_decoding_df["item_id"].unique().tolist())
        n_items_in_bins = len(item_ids_in_bins)
        group_width = 0.86
        bar_width = group_width / max(n_items_in_bins, 1)

        for idx, item_id in enumerate(item_ids_in_bins):
            item_binned = first_decoding_df[first_decoding_df["item_id"] == item_id].sort_values("bin_center")
            x_positions = (
                item_binned["bin_center"].to_numpy(dtype=float)
                + (idx - (n_items_in_bins - 1) / 2.0) * bar_width
            )
            included_mask = item_binned["kind"].eq("included")
            excluded_mask = item_binned["kind"].eq("excluded")
            if included_mask.any():
                ax_first_decoding.bar(
                    x_positions[included_mask.to_numpy()],
                    item_binned.loc[included_mask, "n_points"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.75,
                    edgecolor="black",
                    linewidth=0.4,
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} included",
                )
            if excluded_mask.any():
                ax_first_decoding.bar(
                    x_positions[excluded_mask.to_numpy()],
                    item_binned.loc[excluded_mask, "n_points"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.95,
                    edgecolor="black",
                    linewidth=0.4,
                    hatch="//",
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} excluded",
                )

    ax_first_decoding.set_title(PLOT_NAME_MAP["cluster_size_at_first_decoding"])
    ax_first_decoding.set_xlabel(X_AXIS_NAME_MAP["cluster_size_at_first_decoding"])
    ax_first_decoding.set_ylabel(Y_AXIS_NAME_MAP["cluster_size_at_first_decoding"])
    ax_first_decoding.set_xlim(0.5, visual_infinity_x + 1.0)
    ax_first_decoding.grid(True, alpha=0.25, linestyle="--")
    first_decoding_ticks = sorted(set(range(1, max_finite_first_decoding + 1)) | {int(visual_infinity_x)})
    ax_first_decoding.set_xticks(first_decoding_ticks)
    ax_first_decoding.set_xticklabels([str(x_tick) for x_tick in first_decoding_ticks])
    ax_first_decoding.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.94, FIGURE_DESCRIPTION, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.05, 0.97, 0.88))

    plt.show()

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=200)
    fig.savefig(output_path_svg)
    plt.close(fig)


if __name__ == "__main__":
    main()
