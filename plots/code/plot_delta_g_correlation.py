################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

DB_PATH = "../../database/dspl.db"

RUN_LABEL = "barcode01_agilent_alignment_decoding"

ITEM_IDS_TO_PLOT = [0, 1]

#USE_SINGLE_RUN_ID_INSTEAD = "decoding_1"  # None | str
USE_SINGLE_RUN_ID_INSTEAD = None

MIN_PERFECT_DECODED_RATIO = 0.5

VISUAL_INFINITY_FACTOR = 1.08

DELTA_G_PRECISION = 1.5

MIN_POINTS_PER_BIN_ACCEPTED = 10
MIN_POINTS_PER_BIN_RATIO = 2

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
        AND (? IS NULL OR m.dec_run_id = ?)
    GROUP BY m.exp_id, m.dec_run_id, m.item_id
    HAVING SUM(
        CASE
            WHEN NULLIF(TRIM(CAST(m.dec_pass_id AS TEXT)), '') IS NOT NULL THEN 1
            ELSE 0
        END
    ) > 0
)
SELECT
        m.item_id,
        m.region_id,
        m.position_id,
        m.delta_g,
    m.count_at_first_decoding,
    m.hamming_dist_to_ref
FROM metrics_view_run_cluster_at_decoding m
JOIN eligible_run_item e
    ON e.exp_id = m.exp_id
   AND e.dec_run_id = m.dec_run_id
   AND e.item_id = m.item_id
ORDER BY m.item_id ASC, m.region_id ASC, m.position_id ASC
"""

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/delta_g_correlation.png"
OUTPUT_PATH_SVG = "../plots/delta_g_correlation.svg"

# ---- NAMES ----

FIGURE_TITLE = "Cluster size at decoding against delta G."

FIGURE_DESCRIPTION = f"Cluster size at decoding against delta G, for label {RUN_LABEL}."

PLOT_NAME_MAP = {
    "delta_g_correlation": "Coverage at first correct decoding against delta G (zeros shown at infinity)",
    "delta_g_binned_regression": "Binned mean positive coverage over delta G intervals (included references)",
    "delta_g_excluded_ratio": "Excluded references ratio by delta G bin",
}

X_AXIS_NAME_MAP = {
    "delta_g": "Delta G",
}

Y_AXIS_NAME_MAP = {
    "count_at_first_decoding": "Coverage at first correct decoding (mean of positive runs)",
    "count_at_first_decoding_binned": "Mean positive coverage in bin",
    "excluded_ratio_percent": "Excluded references (%)",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green (fallback)
}

PLOT_STYLE_MAP = {
    "delta_g_correlation": {
        "marker": "o",
        "s": 18,
        "alpha": 0.45,
        "edgecolors": "none",
    },
    "delta_g_correlation_excluded": {
        "marker": "x",
        "s": 24,
        "alpha": 0.80,
        "linewidths": 1.2,
    },
    "delta_g_binned_regression": {
        "group_width_ratio": 0.86,
        "alpha": 0.85,
        "edgecolor": "black",
        "linewidth": 0.4,
    }
}


##############
#### CODE ####
##############

import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def main() -> None:
    db_path = _resolve_path_from_script(Path(DB_PATH))
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_svg = _resolve_path_from_script(Path(OUTPUT_PATH_SVG))

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    query = SQL_QUERY.format(item_placeholders=item_placeholders)
    query_params = [RUN_LABEL, *ITEM_IDS_TO_PLOT, USE_SINGLE_RUN_ID_INSTEAD, USE_SINGLE_RUN_ID_INSTEAD]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()

    if df.empty:
        raise ValueError("No data returned by SQL query. Check RUN_LABEL and ITEM_IDS_TO_PLOT.")

    df["delta_g"] = pd.to_numeric(df["delta_g"], errors="coerce")
    df["count_at_first_decoding"] = (
        pd.to_numeric(df["count_at_first_decoding"], errors="coerce").fillna(0)
    )
    df = df.dropna(subset=["delta_g"])

    # One point per reference and item: mean of strictly positive runs only.
    reference_level_df = (
        df.groupby(["item_id", "region_id", "position_id"], as_index=False)
        .agg(
            delta_g=("delta_g", "first"),
            n_runs_total=("count_at_first_decoding", "size"),
            n_runs_zero=("count_at_first_decoding", lambda s: int((pd.to_numeric(s, errors="coerce") <= 0).sum())),
            mean_positive_coverage=(
                "count_at_first_decoding",
                lambda s: pd.to_numeric(s, errors="coerce").loc[pd.to_numeric(s, errors="coerce") > 0].mean(),
            ),
        )
    )
    reference_level_df["zero_ratio"] = (
        reference_level_df["n_runs_zero"] / reference_level_df["n_runs_total"]
    )

    accepted_reference_level_df = reference_level_df[
        reference_level_df["zero_ratio"] < MIN_PERFECT_DECODED_RATIO
    ].copy()
    rejected_reference_level_df = reference_level_df[
        reference_level_df["zero_ratio"] >= MIN_PERFECT_DECODED_RATIO
    ].copy()

    rejected_reference_level_df["item_name"] = rejected_reference_level_df["item_id"].map(
        ITEM_ID_NAME_MAP
    )
    rejected_reference_level_df["rejection_reason"] = (
        "zero_ratio_above_threshold"
    )

    print(
        "Rejected references summary:"
        f" total={len(rejected_reference_level_df)},"
        f" accepted={len(accepted_reference_level_df)},"
        f" threshold={MIN_PERFECT_DECODED_RATIO:.2f}"
    )
    if not rejected_reference_level_df.empty:
        for item_id in sorted(ITEM_ID_NAME_MAP.keys()):
            n_item_rejected = int((rejected_reference_level_df["item_id"] == item_id).sum())
            if n_item_rejected > 0:
                print(
                    f"  item_id={item_id}"
                    f" ({ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')})"
                    f": rejected={n_item_rejected}"
                )

    fig, (ax_scatter, ax_hist, ax_ratio) = plt.subplots(
        3,
        1,
        figsize=(11.5, 13.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.3, 1.3]},
    )
    scatter_style = PLOT_STYLE_MAP["delta_g_correlation"]
    excluded_scatter_style = PLOT_STYLE_MAP["delta_g_correlation_excluded"]
    hist_style = PLOT_STYLE_MAP["delta_g_binned_regression"]
    binned_rows: list[pd.DataFrame] = []
    ratio_rows: list[pd.DataFrame] = []

    finite_means = accepted_reference_level_df["mean_positive_coverage"].dropna()
    visual_infinity_y = (
        float(finite_means.max()) * VISUAL_INFINITY_FACTOR if not finite_means.empty else 1.0
    )

    for item_id, item_name in ITEM_ID_NAME_MAP.items():
        accepted_item_df = accepted_reference_level_df[
            accepted_reference_level_df["item_id"] == item_id
        ]
        excluded_item_df = rejected_reference_level_df[
            rejected_reference_level_df["item_id"] == item_id
        ]

        if not accepted_item_df.empty:
            ax_scatter.scatter(
                accepted_item_df["delta_g"],
                accepted_item_df["mean_positive_coverage"],
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                label=f"{item_name} accepted (n={len(accepted_item_df)})",
                marker=scatter_style["marker"],
                s=scatter_style["s"],
                alpha=scatter_style["alpha"],
                edgecolors=scatter_style["edgecolors"],
            )

        if not excluded_item_df.empty:
            ax_scatter.scatter(
                excluded_item_df["delta_g"],
                np.full(len(excluded_item_df), visual_infinity_y),
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                label=f"{item_name} excluded (n={len(excluded_item_df)})",
                marker=excluded_scatter_style["marker"],
                s=excluded_scatter_style["s"],
                alpha=excluded_scatter_style["alpha"],
                linewidths=excluded_scatter_style["linewidths"],
            )

        # Upper-floor binning: values in (u-precision, u] belong to bin upper bound u.
        item_binned = accepted_item_df.copy()
        if not item_binned.empty:
            item_binned["bin_upper"] = (
                np.ceil(item_binned["delta_g"] / DELTA_G_PRECISION) * DELTA_G_PRECISION
            )
            item_binned = (
                item_binned.groupby("bin_upper", as_index=False)
                .agg(
                    n_points=("delta_g", "size"),
                    mean_coverage=("mean_positive_coverage", "mean"),
                )
            )
            item_binned = item_binned[
                item_binned["n_points"] >= MIN_POINTS_PER_BIN_ACCEPTED
            ].copy()
            if not item_binned.empty:
                item_binned["item_id"] = item_id
                item_binned["bin_center"] = item_binned["bin_upper"] - (DELTA_G_PRECISION / 2.0)
                binned_rows.append(item_binned)
    if binned_rows:
        binned_df = pd.concat(binned_rows, ignore_index=True)
        item_ids_in_bins = sorted(binned_df["item_id"].unique().tolist())
        n_items_in_bins = len(item_ids_in_bins)
        group_width = DELTA_G_PRECISION * float(hist_style["group_width_ratio"])
        bar_width = group_width / max(n_items_in_bins, 1)

        for idx, item_id in enumerate(item_ids_in_bins):
            item_binned = binned_df[binned_df["item_id"] == item_id].sort_values("bin_center")
            x_positions = (
                item_binned["bin_center"].to_numpy(dtype=float)
                + (idx - (n_items_in_bins - 1) / 2.0) * bar_width
            )
            ax_hist.bar(
                x_positions,
                item_binned["mean_coverage"],
                width=bar_width,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                alpha=hist_style["alpha"],
                edgecolor=hist_style["edgecolor"],
                linewidth=hist_style["linewidth"],
                label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} bins>= {MIN_POINTS_PER_BIN_ACCEPTED}",
                align="center",
            )

    ax_scatter.set_title(PLOT_NAME_MAP["delta_g_correlation"])
    ax_scatter.set_ylabel(Y_AXIS_NAME_MAP["count_at_first_decoding"])
    ax_scatter.grid(True, alpha=0.25, linestyle="--")
    ax_scatter.axhline(visual_infinity_y, color="black", linewidth=1.0, linestyle=":")
    ax_scatter.set_ylim(0, visual_infinity_y * 1.08)
    ax_scatter.legend()

    ax_hist.set_title(PLOT_NAME_MAP["delta_g_binned_regression"])
    ax_hist.set_ylabel(Y_AXIS_NAME_MAP["count_at_first_decoding_binned"])
    ax_hist.grid(True, alpha=0.25, linestyle="--")
    if binned_rows:
        ax_hist.legend()

    for item_id, _item_name in ITEM_ID_NAME_MAP.items():
        accepted_item_df = accepted_reference_level_df[
            accepted_reference_level_df["item_id"] == item_id
        ].copy()
        rejected_item_df = rejected_reference_level_df[
            rejected_reference_level_df["item_id"] == item_id
        ].copy()

        if accepted_item_df.empty and rejected_item_df.empty:
            continue

        if not accepted_item_df.empty:
            accepted_item_df["bin_upper"] = (
                np.ceil(accepted_item_df["delta_g"] / DELTA_G_PRECISION) * DELTA_G_PRECISION
            )
            accepted_counts = (
                accepted_item_df.groupby("bin_upper", as_index=False)
                .size()
                .rename(columns={"size": "n_accepted"})
            )
        else:
            accepted_counts = pd.DataFrame(columns=["bin_upper", "n_accepted"])

        if not rejected_item_df.empty:
            rejected_item_df["bin_upper"] = (
                np.ceil(rejected_item_df["delta_g"] / DELTA_G_PRECISION) * DELTA_G_PRECISION
            )
            rejected_counts = (
                rejected_item_df.groupby("bin_upper", as_index=False)
                .size()
                .rename(columns={"size": "n_excluded"})
            )
        else:
            rejected_counts = pd.DataFrame(columns=["bin_upper", "n_excluded"])

        merged_counts = accepted_counts.merge(rejected_counts, on="bin_upper", how="outer").fillna(0)
        merged_counts["n_accepted"] = merged_counts["n_accepted"].astype(int)
        merged_counts["n_excluded"] = merged_counts["n_excluded"].astype(int)
        merged_counts["n_total"] = merged_counts["n_accepted"] + merged_counts["n_excluded"]
        merged_counts = merged_counts[merged_counts["n_total"] >= MIN_POINTS_PER_BIN_RATIO].copy()
        if merged_counts.empty:
            continue

        merged_counts["excluded_ratio_percent"] = (
            100.0 * merged_counts["n_excluded"] / merged_counts["n_total"]
        )
        merged_counts["bin_center"] = merged_counts["bin_upper"] - (DELTA_G_PRECISION / 2.0)
        merged_counts["item_id"] = item_id
        ratio_rows.append(merged_counts)

    if ratio_rows:
        ratio_df = pd.concat(ratio_rows, ignore_index=True)
        item_ids_in_ratio = sorted(ratio_df["item_id"].unique().tolist())
        n_items_in_ratio = len(item_ids_in_ratio)
        group_width_ratio = DELTA_G_PRECISION * float(hist_style["group_width_ratio"])
        bar_width_ratio = group_width_ratio / max(n_items_in_ratio, 1)

        for idx, item_id in enumerate(item_ids_in_ratio):
            item_ratio = ratio_df[ratio_df["item_id"] == item_id].sort_values("bin_center")
            x_positions = (
                item_ratio["bin_center"].to_numpy(dtype=float)
                + (idx - (n_items_in_ratio - 1) / 2.0) * bar_width_ratio
            )
            ax_ratio.bar(
                x_positions,
                item_ratio["excluded_ratio_percent"],
                width=bar_width_ratio,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                alpha=hist_style["alpha"],
                edgecolor=hist_style["edgecolor"],
                linewidth=hist_style["linewidth"],
                label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} bins>= {MIN_POINTS_PER_BIN_RATIO}",
                align="center",
            )

    ax_ratio.set_title(PLOT_NAME_MAP["delta_g_excluded_ratio"])
    ax_ratio.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
    ax_ratio.set_ylabel(Y_AXIS_NAME_MAP["excluded_ratio_percent"])
    ax_ratio.set_ylim(0, 100)
    ax_ratio.grid(True, alpha=0.25, linestyle="--")
    if ratio_rows:
        ax_ratio.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.948, FIGURE_DESCRIPTION, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.04, 0.97, 0.9))

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=250)
    fig.savefig(output_path_svg)
    plt.close(fig)


if __name__ == "__main__":
    main()
