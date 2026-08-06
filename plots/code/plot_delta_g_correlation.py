################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

# "barcode01_agilent" | "barcode02_dynegene" | "barcode03_genscript" | "barcode04_six_images"
BARCODE = "barcode03_genscript"

DB_PATH = f"/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/{BARCODE}.db"

RUN_LABEL = f"{BARCODE}_alignment_decoding"

#RUN_LABEL = "test_label"

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
        2: "motif-paircode",
    }

#USE_SINGLE_RUN_ID_INSTEAD = "decoding_1"  # None | str
USE_SINGLE_RUN_ID_INSTEAD = None

# Exclusion rule:
# a reference is excluded when its zero-decoding run ratio is >= this threshold.
# Ratio scale is [0, 1]. Example: 0.5 means "exclude references with at least 50% zero runs".
# And 1.0 means "exclude references with 100% zero runs only", i.e. dropped out references.
MAX_ZERO_RUN_RATIO_FOR_INCLUSION = 0.5

VISUAL_INFINITY_FACTOR = 1.08

DELTA_G_PRECISION = 0.33333

MIN_POINTS_PER_BIN_ACCEPTED = 10
MIN_POINTS_PER_BIN_RATIO = 2

DISPLAY_LIN_REG = True

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
    m.dec_run_id,
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

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*) AS n_runs_total
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""

# ---- OUTPUT PATH ----

OUTPUT_PATH_FIG_COVERAGE_PNG = f"../plots/{BARCODE}/delta_g_scatter_and_effect_on_coverage_at_oligo_decoding.png"
OUTPUT_PATH_FIG_COVERAGE_PDF = f"../plots/{BARCODE}/delta_g_scatter_and_effect_on_coverage_at_oligo_decoding.pdf"

OUTPUT_PATH_FIG_DROPOUT_PNG = f"../plots/{BARCODE}/delta_g_effect_on_underdecoded_references.png"
OUTPUT_PATH_FIG_DROPOUT_PDF = f"../plots/{BARCODE}/delta_g_effect_on_underdecoded_references.pdf"

# ---- NAMES ----

FIGURE_COVERAGE_TITLE = f"Reference Coverage at reference decoding against delta G ({BARCODE})."

FIGURE_COVERAGE_DESCRIPTION_BASE = (
    "Included and drop out references form a full partition (100% of selected references).\n"
)

FIGURE_DROPOUT_TITLE = f"Delta G effect on Underdecoded Reference Probability ({BARCODE})."

FIGURE_DROPOUT_DESCRIPTION_BASE = (
    "$\\mathbf{DEFINITION}$: a reference is said to be $\\mathbf{underdecoded}$ iff it is perfectly decoded during less than "
    f"{MAX_ZERO_RUN_RATIO_FOR_INCLUSION*100}% of the decoding runs."
    "Otherwise, it is said to be $\\mathbf{included}$.\n"

)

PLOT_NAME_MAP = {
    "delta_g_correlation": "Reference Coverage at first perfect decoding against delta G ($\\mathbf{drop\\ out}$ = never decoded references shown at infinity on the top)",
    "delta_g_binned_regression": "Binned mean positive coverage over delta G intervals ($\\mathbf{included}$ references)",
    "delta_g_excluded_ratio": "$\\mathbf{Underdecoded\\ references}$ ratio by Delta G bin",
}

X_AXIS_NAME_MAP = {
    "delta_g": "Delta G",
}

Y_AXIS_NAME_MAP = {
    "count_at_first_decoding": "Reference Coverage at first correct decoding (mean over positive runs)",
    "count_at_first_decoding_binned": "Mean coverage in bin",
    "excluded_ratio_probability": "Probability for a reference to be $\\mathbf{underdecoded}$",
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
    },
    "delta_g_dropout_probability": {
        "bar_width_ratio": 1.0,
        "alpha": 0.9,
        "edgecolor": "none",
        "linewidth": 0.0,
    },
    "delta_g_dropout_lin_reg": {
        "linestyle": "-",
        "linewidth": 2.0,
        "alpha": 0.95,
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
    output_path_fig_coverage_png = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_COVERAGE_PNG))
    output_path_fig_coverage_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_COVERAGE_PDF))
    output_path_fig_dropout_png = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_DROPOUT_PNG))
    output_path_fig_dropout_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_DROPOUT_PDF))

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    query = SQL_QUERY.format(item_placeholders=item_placeholders)
    query_params = [RUN_LABEL, *ITEM_IDS_TO_PLOT, USE_SINGLE_RUN_ID_INSTEAD, USE_SINGLE_RUN_ID_INSTEAD]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)
        n_runs_total = int(
            conn.execute(SQL_QUERY_LABEL_RUNS_TOTAL, (RUN_LABEL,)).fetchone()[0]
        )

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()

    if df.empty:
        raise ValueError("No data returned by SQL query. Check RUN_LABEL and ITEM_IDS_TO_PLOT.")

    n_runs_displayed = int(df["dec_run_id"].nunique())
    volumetry_suffix = (
        f" Label '{RUN_LABEL}': decoding runs studied={n_runs_total}, displayed={n_runs_displayed}."
    )
    figure_coverage_description = FIGURE_COVERAGE_DESCRIPTION_BASE + volumetry_suffix
    figure_dropout_description = FIGURE_DROPOUT_DESCRIPTION_BASE + volumetry_suffix

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
        reference_level_df["zero_ratio"] < MAX_ZERO_RUN_RATIO_FOR_INCLUSION
    ].copy()
    rejected_reference_level_df = reference_level_df[
        reference_level_df["zero_ratio"] >= MAX_ZERO_RUN_RATIO_FOR_INCLUSION
    ].copy()

    # Figure 1 semantics: drop out means "never decoded" only.
    never_decoded_reference_level_df = reference_level_df[
        reference_level_df["n_runs_zero"] == reference_level_df["n_runs_total"]
    ].copy()
    decoded_at_least_once_reference_level_df = reference_level_df[
        reference_level_df["n_runs_zero"] < reference_level_df["n_runs_total"]
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
        f" threshold={MAX_ZERO_RUN_RATIO_FOR_INCLUSION:.2f}"
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

    fig_coverage, (ax_scatter, ax_hist) = plt.subplots(
        2,
        1,
        figsize=(11.5, 10.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.3]},
    )
    scatter_style = PLOT_STYLE_MAP["delta_g_correlation"]
    excluded_scatter_style = PLOT_STYLE_MAP["delta_g_correlation_excluded"]
    hist_style = PLOT_STYLE_MAP["delta_g_binned_regression"]
    binned_rows: list[pd.DataFrame] = []

    finite_means = decoded_at_least_once_reference_level_df["mean_positive_coverage"].dropna()
    visual_infinity_y = (
        float(finite_means.max()) * VISUAL_INFINITY_FACTOR if not finite_means.empty else 1.0
    )

    for item_id, item_name in ITEM_ID_NAME_MAP.items():
        accepted_item_df = decoded_at_least_once_reference_level_df[
            decoded_at_least_once_reference_level_df["item_id"] == item_id
        ]
        excluded_item_df = never_decoded_reference_level_df[
            never_decoded_reference_level_df["item_id"] == item_id
        ]

        if not accepted_item_df.empty:
            ax_scatter.scatter(
                accepted_item_df["delta_g"],
                accepted_item_df["mean_positive_coverage"],
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                label=f"{item_name} "+"$\\mathbf{included}$"+ f"references (count={len(accepted_item_df)})",
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
                label=f"{item_name} "+"$\\mathbf{drop\\ out}$"+ f" references (count={len(excluded_item_df)})",
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
    ax_scatter.axhline(
        visual_infinity_y,
        color="black",
        linewidth=1.0,
        linestyle=":",
        label="Dropped-out references",
    )
    finite_scatter_values = decoded_at_least_once_reference_level_df["mean_positive_coverage"].dropna().to_numpy(dtype=float)
    if finite_scatter_values.size > 0:
        y_min_data = float(np.min(finite_scatter_values))
        y_max_data = float(np.max(finite_scatter_values))
        y_span = max(y_max_data - y_min_data, 1e-9)
        y_lower = y_min_data - (0.08 * y_span)
        y_upper = visual_infinity_y + (0.08 * y_span)
        if y_upper <= visual_infinity_y:
            y_upper = visual_infinity_y + 0.5
        ax_scatter.set_ylim(y_lower, y_upper)
    else:
        ax_scatter.set_ylim(visual_infinity_y - 1.0, visual_infinity_y + 0.5)

    y_ticks = [
        float(tick)
        for tick in ax_scatter.get_yticks()
        if float(tick) < float(visual_infinity_y)
    ]
    y_ticks.append(float(visual_infinity_y))
    y_ticks = sorted(set(y_ticks))
    y_labels = ["∞" if abs(float(tick) - float(visual_infinity_y)) <= 1e-9 else f"{tick:g}" for tick in y_ticks]
    ax_scatter.set_yticks(y_ticks)
    ax_scatter.set_yticklabels(y_labels)
    ax_scatter.legend()

    ax_hist.set_title(PLOT_NAME_MAP["delta_g_binned_regression"])
    ax_hist.set_ylabel(Y_AXIS_NAME_MAP["count_at_first_decoding_binned"])
    ax_hist.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
    ax_hist.grid(True, alpha=0.25, linestyle="--")
    if binned_rows:
        ax_hist.legend()

    fig_coverage.suptitle(FIGURE_COVERAGE_TITLE, fontsize=16, y=0.985)
    fig_coverage.text(
        0.5,
        0.948,
        figure_coverage_description,
        ha="center",
        va="top",
        wrap=True,
        fontsize=11,
    )
    fig_coverage.tight_layout(rect=(0.03, 0.04, 0.97, 0.9))

    plt.show()

    output_path_fig_coverage_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_fig_coverage_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig_coverage.savefig(output_path_fig_coverage_png, dpi=250)
    fig_coverage.savefig(output_path_fig_coverage_pdf)
    plt.close(fig_coverage)

    ratio_rows_by_item: dict[int, pd.DataFrame] = {}

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

        merged_counts["excluded_ratio_probability"] = (
            merged_counts["n_excluded"] / merged_counts["n_total"]
        )
        merged_counts["bin_center"] = merged_counts["bin_upper"] - (DELTA_G_PRECISION / 2.0)
        merged_counts["item_id"] = item_id
        ratio_rows_by_item[item_id] = merged_counts

    item_ids_for_dropout = [item_id for item_id in ITEM_ID_NAME_MAP if item_id in ITEM_IDS_TO_PLOT]
    if not item_ids_for_dropout:
        raise ValueError("No item id selected for dropout plot. Check ITEM_IDS_TO_PLOT.")

    fig_dropout, axes_dropout = plt.subplots(
        len(item_ids_for_dropout),
        1,
        figsize=(11.5, 4.8 * len(item_ids_for_dropout)),
        sharex=False,
    )
    if len(item_ids_for_dropout) == 1:
        axes_dropout = [axes_dropout]

    dropout_style = PLOT_STYLE_MAP["delta_g_dropout_probability"]
    lin_reg_style = PLOT_STYLE_MAP["delta_g_dropout_lin_reg"]

    for ax_item_dropout, item_id in zip(axes_dropout, item_ids_for_dropout):
        item_name = ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}")
        item_ratio = ratio_rows_by_item.get(item_id)
        n_underdecoded_plotted = int(item_ratio["n_excluded"].sum()) if item_ratio is not None and not item_ratio.empty else 0
        if item_ratio is None or item_ratio.empty:
            ax_item_dropout.text(
                0.5,
                0.5,
                "No bins pass threshold",
                ha="center",
                va="center",
                transform=ax_item_dropout.transAxes,
            )
        else:
            bar_width_ratio = DELTA_G_PRECISION * float(dropout_style["bar_width_ratio"])
            x_data = item_ratio["bin_center"].to_numpy(dtype=float)
            y_data = item_ratio["excluded_ratio_probability"].to_numpy(dtype=float)
            ax_item_dropout.bar(
                x_data,
                y_data,
                width=bar_width_ratio,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                alpha=dropout_style["alpha"],
                edgecolor=dropout_style["edgecolor"],
                linewidth=dropout_style["linewidth"],
                label=f"{item_name} bins>= {MIN_POINTS_PER_BIN_RATIO}",
                align="center",
            )

            if DISPLAY_LIN_REG and len(x_data) >= 2:
                slope, intercept = np.polyfit(x_data, y_data, deg=1)
                y_pred = (slope * x_data) + intercept
                ss_res = float(np.sum((y_data - y_pred) ** 2))
                ss_tot = float(np.sum((y_data - float(np.mean(y_data))) ** 2))
                r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

                x_line = np.linspace(float(np.min(x_data)), float(np.max(x_data)), 200)
                y_line = (slope * x_line) + intercept
                ax_item_dropout.plot(
                    x_line,
                    y_line,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    linestyle=lin_reg_style["linestyle"],
                    linewidth=lin_reg_style["linewidth"],
                    alpha=lin_reg_style["alpha"],
                    label=f"lin reg: coef={slope:.4f}, R²={r_squared:.4f}",
                )

            ax_item_dropout.legend()

        ax_item_dropout.set_title(
            f"{PLOT_NAME_MAP['delta_g_excluded_ratio']} - {item_name} "
            f"(N underdecoded refs plotted={n_underdecoded_plotted})"
        )
        ax_item_dropout.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
        ax_item_dropout.set_ylabel(Y_AXIS_NAME_MAP["excluded_ratio_probability"])
        ax_item_dropout.grid(True, alpha=0.25, linestyle="--")

    fig_dropout.suptitle(FIGURE_DROPOUT_TITLE, fontsize=16, y=0.985)
    fig_dropout.text(
        0.5,
        0.955,
        figure_dropout_description,
        ha="center",
        va="top",
        wrap=True,
        fontsize=11,
    )
    fig_dropout.tight_layout(rect=(0.03, 0.04, 0.97, 0.92))

    plt.show()

    output_path_fig_dropout_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_fig_dropout_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig_dropout.savefig(output_path_fig_dropout_png, dpi=250)
    fig_dropout.savefig(output_path_fig_dropout_pdf)
    plt.close(fig_dropout)


if __name__ == "__main__":
    main()
