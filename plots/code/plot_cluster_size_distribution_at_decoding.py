################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

# "barcode01_agilent" | "barcode02_dynegene" | "barcode03_genscript" | "barcode04_six_images"
BARCODE = "barcode04_six_images"

DB_PATH = f"/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/{BARCODE}.db"

RUN_LABEL = f"{BARCODE}_alignment_decoding"

if BARCODE == "barcode04_six_images":
    ITEM_IDS_TO_PLOT = [
        0,
        1,
        3,
        4,
        5
    ]

    ITEM_ID_NAME_MAP = {
        0: "Chest",
        1: "Woman",
        2: "Burger",
        3: "Bird",
        4: "Night",
        5: "Day"
    }

else:
    ITEM_IDS_TO_PLOT = [
        0,
        1,
    ]

    ITEM_ID_NAME_MAP = {
        0: "JPEGDNA-reference",
        1: "JPEGDNA-delta-G",
    }

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = f"../plots/{BARCODE}/cluster_size_distribution_at_decoding.png"
OUTPUT_PATH_PDF = f"../plots/{BARCODE}/cluster_size_distribution_at_decoding.pdf"


VISUAL_INFINITY_FACTOR = 1.08
FIRST_DECODING_BIN_WIDTH = 0.10
MAX_XTICKS_SECOND_PLOT = 18

# First subplot coverage display window (inclusive).
COVERAGE_DISPLAY_MIN = 0
COVERAGE_DISPLAY_MAX = 50
MAX_XTICKS_FIRST_PLOT = 16

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

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*) AS n_runs_total
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""



# ---- NAMES ----

FIGURE_TITLE = f"Cluster size distributions for items JPEG DNA and JPEG DNA delta G over multiple runs ({BARCODE})"

FIGURE_DESCRIPTION_BASE = (
    f"Average normalized cluster size distribution at decoding over all runs of label {RUN_LABEL}. "
    "For each item, $\\mu_{\\text{cov}}$ and $\\sigma_{\\text{cov}}$ are computed from the aggregated "
    "normalized distribution over coverage bins."
)

PLOT_NAME_MAP = {
    "cluster_size_distribution": "Reference Coverage distribution at $\\mathbf{image\\ decoding}$",
    "cluster_size_at_first_decoding": "Reference Coverage distribution at $\\mathbf{reference\\ decoding}$ (dropped-out references at $\\infty$)",
}

X_AXIS_NAME_MAP = {
    "cluster_size": "Coverage (i.e. size of the reference's cluster)",
    "cluster_size_at_first_decoding": "Coverage of a reference at its first perfect decoding",
}

Y_AXIS_NAME_MAP = {
    "cluster_size_distribution": "Normalized number of references",
    "cluster_size_at_first_decoding": "Normalized number of references in bin",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green
    3: "#ff7f0e",  # orange
    4: "#9467bd",  # purple
    5: "#8c564b",  # brown
}

PLOT_STYLE_MAP = {
    "cluster_size_distribution": {
        "alpha": 0.85,
        "bar_width": 0.5,
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
) -> tuple[list[int], dict[int, np.ndarray], dict[int, tuple[float, float, float]]]:
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

    run_item_histograms = run_item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)

    normalized_by_item: dict[int, np.ndarray] = {}
    coverage_moments_by_item: dict[int, tuple[float, float, float]] = {}

    for item_id in sorted(df["item_id"].unique().tolist()):
        if item_id not in ITEM_ID_NAME_MAP:
            continue
        item_histograms = run_item_histograms.xs(item_id, level="item_id")
        if isinstance(item_histograms, pd.Series):
            item_histograms = item_histograms.to_frame().T
        item_histograms = item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)

        total_refs_item = float(item_histograms.to_numpy(dtype=float).sum())
        if total_refs_item <= 0:
            normalized_distribution = np.zeros(len(cluster_sizes), dtype=float)
        else:
            aggregated_counts = item_histograms.sum(axis=0).to_numpy(dtype=float)
            normalized_distribution = aggregated_counts / total_refs_item

        normalized_by_item[item_id] = normalized_distribution

        support = np.array(cluster_sizes, dtype=float)
        total_prob = float(normalized_distribution.sum())
        if total_prob <= 0:
            coverage_moments_by_item[item_id] = (0.0, 0.0, 0.0)
        else:
            prob_distribution = normalized_distribution / total_prob
            mu_cov = float(np.sum(support * prob_distribution))
            var_cov = float(np.sum(((support - mu_cov) ** 2) * prob_distribution))
            sigma_cov = float(np.sqrt(max(var_cov, 0.0)))
            cdf = np.cumsum(prob_distribution)
            median_idx = int(np.searchsorted(cdf, 0.5, side="left"))
            median_idx = min(max(median_idx, 0), len(support) - 1)
            m_cov = float(support[median_idx])
            coverage_moments_by_item[item_id] = (mu_cov, sigma_cov, m_cov)

    return cluster_sizes, normalized_by_item, coverage_moments_by_item


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
    reference_level_df["zero_ratio_percent"] = 100.0 * reference_level_df["zero_ratio"]

    included_df = reference_level_df[
        reference_level_df["n_runs_zero"] < reference_level_df["n_runs_total"]
    ].copy()
    excluded_df = reference_level_df[
        reference_level_df["n_runs_zero"] == reference_level_df["n_runs_total"]
    ].copy()

    return reference_level_df, included_df, excluded_df


def main() -> None:
    db_path = _resolve_path_from_script(Path(DB_PATH))
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_PDF))

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    query = SQL_QUERY.format(item_placeholders=item_placeholders)
    query_params = [RUN_LABEL, *ITEM_IDS_TO_PLOT]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)
        n_runs_total = int(
            conn.execute(SQL_QUERY_LABEL_RUNS_TOTAL, (RUN_LABEL,)).fetchone()[0]
        )

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())]

    if df.empty:
        raise ValueError("No data returned by SQL query. Check RUN_LABEL and ITEM_IDS_TO_PLOT.")

    n_runs_displayed = int(df["dec_run_id"].nunique())
    figure_description = (
        f"{FIGURE_DESCRIPTION_BASE} "
        f"Label '{RUN_LABEL}': decoding runs studied={n_runs_total}, displayed={n_runs_displayed}."
    )

    cluster_sizes, normalized_by_item, coverage_moments_by_item = _compute_distribution_stats(df)
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
    if COVERAGE_DISPLAY_MIN > COVERAGE_DISPLAY_MAX:
        raise ValueError("COVERAGE_DISPLAY_MIN must be <= COVERAGE_DISPLAY_MAX.")
    x_in_range_mask = (x >= float(COVERAGE_DISPLAY_MIN)) & (x <= float(COVERAGE_DISPLAY_MAX))
    if not np.any(x_in_range_mask):
        raise ValueError(
            "No coverage bins in first subplot range. "
            "Adjust COVERAGE_DISPLAY_MIN/COVERAGE_DISPLAY_MAX."
        )
    x_first_plot = x[x_in_range_mask]

    item_ids = sorted(normalized_by_item.keys())
    n_items = len(item_ids)
    group_width_first_plot = 1.0
    bar_width = group_width_first_plot / max(n_items, 1)
    marker_style = PLOT_STYLE_MAP["average_coverage_marker"]

    for idx, item_id in enumerate(item_ids):
        mu_cov, sigma_cov, m_cov = coverage_moments_by_item.get(item_id, (0.0, 0.0, 0.0))
        offsets = x_first_plot + (idx - (n_items - 1) / 2.0) * bar_width
        ax_item_distribution.bar(
            offsets,
            normalized_by_item[item_id][x_in_range_mask],
            width=bar_width,
            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
            alpha=style["alpha"],
            edgecolor="none",
            linewidth=0.0,
            label=(
                f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} "
                f"($\\mu_{{\\text{{cov}}}}={mu_cov:.2f}, "
                f"\\sigma_{{\\text{{cov}}}}={sigma_cov:.2f}, "
                f"m_{{\\text{{cov}}}}={m_cov:.2f}$)"
            ),
        )
        if item_id in coverage_moments_by_item:
            ax_item_distribution.axvline(
                m_cov,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                linestyle=marker_style["linestyle"],
                linewidth=marker_style["linewidth"],
                alpha=marker_style["alpha"],
            )

    ax_item_distribution.set_title(PLOT_NAME_MAP["cluster_size_distribution"])
    ax_item_distribution.set_xlabel(X_AXIS_NAME_MAP["cluster_size"])
    ax_item_distribution.set_ylabel(Y_AXIS_NAME_MAP["cluster_size_distribution"])
    first_plot_ticks = x_first_plot.astype(int).tolist()
    if first_plot_ticks:
        tick_step = max(1, int(np.ceil(len(first_plot_ticks) / float(MAX_XTICKS_FIRST_PLOT))))
        reduced_ticks = first_plot_ticks[::tick_step]
        if reduced_ticks[-1] != first_plot_ticks[-1]:
            reduced_ticks.append(first_plot_ticks[-1])
        ax_item_distribution.set_xticks(reduced_ticks)
        ax_item_distribution.set_xticklabels([str(size) for size in reduced_ticks])
    ax_item_distribution.set_xlim(float(COVERAGE_DISPLAY_MIN) - 0.5, float(COVERAGE_DISPLAY_MAX) + 0.5)
    ax_item_distribution.grid(axis="y", alpha=0.25, linestyle="--")
    ax_item_distribution.legend()

    if FIRST_DECODING_BIN_WIDTH <= 0:
        raise ValueError("FIRST_DECODING_BIN_WIDTH must be > 0.")

    finite_first_decoding = included_reference_df["mean_positive_count_at_first_decoding"].dropna()
    max_finite_first_decoding = (
        int(np.ceil(float(finite_first_decoding.max()))) if not finite_first_decoding.empty else 1
    )
    max_finite_bin_edge = (
        np.ceil(max_finite_first_decoding / FIRST_DECODING_BIN_WIDTH) * FIRST_DECODING_BIN_WIDTH
    )
    visual_infinity_x = float(max_finite_bin_edge + 2.0 * FIRST_DECODING_BIN_WIDTH)
    total_refs_by_item = reference_level_df.groupby("item_id").size().to_dict()

    first_decoding_rows: list[pd.DataFrame] = []

    for item_id in sorted(ITEM_ID_NAME_MAP.keys()):
        item_included = included_reference_df[included_reference_df["item_id"] == item_id].copy()
        item_excluded = excluded_reference_df[excluded_reference_df["item_id"] == item_id].copy()

        if not item_included.empty:
            item_included = item_included[np.isfinite(item_included["mean_positive_count_at_first_decoding"])].copy()
        if not item_included.empty:
            item_included["bin_center"] = (
                np.floor(
                    item_included["mean_positive_count_at_first_decoding"] / FIRST_DECODING_BIN_WIDTH
                )
                * FIRST_DECODING_BIN_WIDTH
                + (FIRST_DECODING_BIN_WIDTH / 2.0)
            )
            item_included = (
                item_included.groupby("bin_center", as_index=False)
                .agg(n_refs=("mean_positive_count_at_first_decoding", "size"))
            )
            total_refs_item = float(total_refs_by_item.get(item_id, 0))
            if total_refs_item <= 0:
                continue
            item_included["n_refs_normalized"] = item_included["n_refs"] / total_refs_item
            item_included["item_id"] = item_id
            item_included["kind"] = "included"
            first_decoding_rows.append(item_included)

        if not item_excluded.empty:
            total_refs_item = float(total_refs_by_item.get(item_id, 0))
            if total_refs_item <= 0:
                continue
            first_decoding_rows.append(
                pd.DataFrame(
                    [
                        {
                            "n_refs": len(item_excluded),
                            "n_refs_normalized": len(item_excluded) / total_refs_item,
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
        group_width = FIRST_DECODING_BIN_WIDTH
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
                    item_binned.loc[included_mask, "n_refs_normalized"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.75,
                    edgecolor="none",
                    linewidth=0.0,
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} decoded at least once",
                )
            if excluded_mask.any():
                ax_first_decoding.bar(
                    x_positions[excluded_mask.to_numpy()],
                    item_binned.loc[excluded_mask, "n_refs_normalized"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.95,
                    edgecolor="none",
                    linewidth=0.0,
                    hatch="//",
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} dropped-out (never decoded)",
                )

    ax_first_decoding.set_title(PLOT_NAME_MAP["cluster_size_at_first_decoding"])
    ax_first_decoding.set_xlabel(X_AXIS_NAME_MAP["cluster_size_at_first_decoding"])
    ax_first_decoding.set_ylabel(Y_AXIS_NAME_MAP["cluster_size_at_first_decoding"])
    if first_decoding_rows:
        finite_centers = first_decoding_df.loc[
            first_decoding_df["bin_center"] < visual_infinity_x,
            "bin_center",
        ]
        if not finite_centers.empty:
            left_xlim = float(finite_centers.min() - FIRST_DECODING_BIN_WIDTH / 2.0)
        else:
            left_xlim = float(visual_infinity_x - FIRST_DECODING_BIN_WIDTH)
    else:
        left_xlim = float(visual_infinity_x - FIRST_DECODING_BIN_WIDTH)
    ax_first_decoding.set_xlim(left_xlim, visual_infinity_x + FIRST_DECODING_BIN_WIDTH)
    ax_first_decoding.grid(True, alpha=0.25, linestyle="--")

    if first_decoding_rows:
        finite_tick_candidates = sorted(
            {
                float(x_tick)
                for x_tick in first_decoding_df.loc[
                    first_decoding_df["bin_center"] < visual_infinity_x,
                    "bin_center",
                ].tolist()
            }
        )
    else:
        finite_tick_candidates = []

    if finite_tick_candidates:
        tick_step = max(1, int(np.ceil(len(finite_tick_candidates) / float(MAX_XTICKS_SECOND_PLOT))))
        finite_ticks = finite_tick_candidates[::tick_step]
        if finite_ticks[-1] != finite_tick_candidates[-1]:
            finite_ticks.append(finite_tick_candidates[-1])
    else:
        finite_ticks = []

    first_decoding_ticks = [*finite_ticks, visual_infinity_x]
    ax_first_decoding.set_xticks(first_decoding_ticks)
    ax_first_decoding.set_xticklabels(
        [
            "∞"
            if np.isclose(x_tick, visual_infinity_x)
            else (f"{x_tick:.1f}" if not np.isclose(x_tick, round(x_tick)) else str(int(round(x_tick))))
            for x_tick in first_decoding_ticks
        ]
    )
    ax_first_decoding.axvline(
        visual_infinity_x,
        color="black",
        linewidth=1.0,
        linestyle=":",
        label="∞ = dropped-out references",
    )
    ax_first_decoding.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.94, figure_description, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.05, 0.97, 0.88))

    plt.show()

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=200)
    fig.savefig(output_path_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()
