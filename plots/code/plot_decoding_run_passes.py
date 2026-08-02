################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

DB_PATH = "../../database/barcode01_agilent.db"

DEC_RUN_ID_TO_PLOT = "decoding_164"

ITEM_IDS_TO_PLOT = [0, 1, 2]

# pass_index | coverage | n_tot_reads | estimated_sequencing_duration | run_duration
X_AXIS_KEY = "estimated_sequencing_duration"

FILTER_OUT_NULL_PSNR_ROWS = True

SQL_QUERY = """
SELECT
        m.exp_id,
        m.dec_run_id,
        m.item_id,
        m.dec_pass_id,
        m.coverage,
        m.hamming_distance_normalized,
        m.perfectly_decoded_payload_ratio,
        m.fpdpc,
        m.psnr,
        p.N_trimmed_reads_cum_all_items,
        p.N_decoding_reads_cum_all_items,
        p.run_duration,
        p.estimated_sequencing_duration
FROM metrics_view_run_item_pass m
JOIN metrics_view_run_pass p
    ON p.exp_id = m.exp_id
   AND p.dec_run_id = m.dec_run_id
   AND p.dec_pass_id = m.dec_pass_id
WHERE m.dec_run_id = ?
    AND m.item_id IN ({item_placeholders})
ORDER BY m.item_id ASC, m.dec_pass_id ASC
"""

ITEM_ID_NAME_MAP = {
    0: "JPEGDNA-reference",
    1: "JPEGDNA-delta-G",
    2: "Motif-paircode"
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/decoding_run_passes.png"
OUTPUT_PATH_SVG = "../plots/decoding_run_passes.svg"

# ---- NAMES ----

FIGURE_TITLE = "Decoding run passes metrics"

FIGURE_DESCRIPTION_BASE = (
    f"Pass-level metrics for decoding run {DEC_RUN_ID_TO_PLOT}. "
    "The X axis is configurable and uses the same pass ordering for every subplot."
)

PLOT_NAME_MAP = {
    "psnr": "PSNR",
    "hamming_distance_normalized": "Normalized Hamming distance",
    "perfectly_decoded_payload_ratio": "Perfectly decoded payload ratio",
    "first_time_perfectly_decoded_payload_count": "First-time perfectly decoded payload count",
}

X_AXIS_NAME_MAP = {
    "pass_index": "Pass index",
    "coverage": "Coverage",
    "n_tot_reads": "Total trimmed reads (cumulative across items)",
    "estimated_sequencing_duration": "Estimated sequencing duration",
    "run_duration": "Run duration",
}

Y_AXIS_NAME_MAP = {
    "psnr": "PSNR",
    "hamming_distance_normalized": "Normalized Hamming distance",
    "perfectly_decoded_payload_ratio": "Perfectly decoded payload ratio (%)",
    "first_time_perfectly_decoded_payload_count": "First perfectly decoded payload count (fpdpc)",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",
    1: "#d62728",
    2: "#2ca02c",
}

PLOT_STYLE_MAP = {
    "psnr": {
        "linestyle": "-",
        "linewidth": 2.2,
        "marker": "o",
        "markersize": 5.5,
    },
    "hamming_distance_normalized": {
        "linestyle": "-",
        "linewidth": 2.2,
        "marker": "o",
        "markersize": 5.5,
    },
    "perfectly_decoded_payload_ratio": {
        "linestyle": "-",
        "linewidth": 2.2,
        "marker": "o",
        "markersize": 5.5,
    },
    "first_time_perfectly_decoded_payload_count": {
        "alpha": 0.78,
        "edgecolor": "black",
        "linewidth": 1,
    },
}


##############
#### CODE ####
##############

import sqlite3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


X_AXIS_COLUMN_MAP = {
    "pass_index": ("dec_pass_id", "item_specific"),
    "coverage": ("coverage", "item_specific"),
    "n_tot_reads": ("N_trimmed_reads_cum_all_items", "shared"),
    "estimated_sequencing_duration": ("estimated_sequencing_duration", "shared"),
    "run_duration": ("run_duration", "shared"),
}

METRIC_COLUMN_MAP = {
    "psnr": "psnr",
    "hamming_distance_normalized": "hamming_distance_normalized",
    "perfectly_decoded_payload_ratio": "perfectly_decoded_payload_ratio",
    "first_time_perfectly_decoded_payload_count": "fpdpc",
}


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def _compute_min_positive_spacing(values: pd.Series) -> float:
    unique_values = np.sort(pd.unique(pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)))
    if unique_values.size < 2:
        return 1.0
    diffs = np.diff(unique_values)
    positive_diffs = diffs[diffs > 0]
    return float(positive_diffs.min()) if positive_diffs.size > 0 else 1.0


def _get_x_column_and_mode() -> tuple[str, str]:
    if X_AXIS_KEY not in X_AXIS_COLUMN_MAP:
        raise ValueError(f"Unsupported X_AXIS_KEY: {X_AXIS_KEY}")
    return X_AXIS_COLUMN_MAP[X_AXIS_KEY]


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    db_path = (base_dir / DB_PATH).resolve()
    output_path_png = (_resolve_path_from_script(Path(OUTPUT_PATH_PNG))).resolve()
    output_path_svg = (_resolve_path_from_script(Path(OUTPUT_PATH_SVG))).resolve()

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    query = SQL_QUERY.format(item_placeholders=item_placeholders)
    query_params = [DEC_RUN_ID_TO_PLOT, *ITEM_IDS_TO_PLOT]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)

    if df.empty:
        raise ValueError("No data returned by SQL query. Check DEC_RUN_ID_TO_PLOT and ITEM_IDS_TO_PLOT.")

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()
    if df.empty:
        raise ValueError("No rows left after filtering ITEM_IDS_TO_PLOT to known item ids.")

    x_column, x_mode = _get_x_column_and_mode()

    for column in [
        "dec_pass_id",
        "coverage",
        "hamming_distance_normalized",
        "perfectly_decoded_payload_ratio",
        "fpdpc",
        "psnr",
        "N_trimmed_reads_cum_all_items",
        "N_decoding_reads_cum_all_items",
        "run_duration",
        "estimated_sequencing_duration",
    ]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    if FILTER_OUT_NULL_PSNR_ROWS:
        df = df[df["psnr"].notna()].copy()

    df["perfectly_decoded_payload_ratio"] = df["perfectly_decoded_payload_ratio"] * 100.0
    df["fpdpc"] = df["fpdpc"].fillna(0)

    df = df.sort_values(["item_id", "dec_pass_id"]).reset_index(drop=True)

    n_items_displayed = int(df["item_id"].nunique())
    figure_description = (
        f"{FIGURE_DESCRIPTION_BASE} "
        f"Displayed items={n_items_displayed}; x-axis={X_AXIS_NAME_MAP[X_AXIS_KEY]}."
    )

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))
    axes_flat = axes.flatten()
    plot_keys = list(PLOT_NAME_MAP.keys())

    for idx, plot_key in enumerate(plot_keys):
        ax = axes_flat[idx]
        style = PLOT_STYLE_MAP[plot_key]
        y_column = METRIC_COLUMN_MAP[plot_key]

        if plot_key == "first_time_perfectly_decoded_payload_count":
            if x_mode == "shared":
                shared_spacing = _compute_min_positive_spacing(df[x_column])
                group_width = 0.8 * shared_spacing
                bar_width = group_width / max(len(ITEM_IDS_TO_PLOT), 1)
            else:
                bar_width = 0.8 * _compute_min_positive_spacing(df[x_column])

            for item_index, item_id in enumerate(ITEM_IDS_TO_PLOT):
                item_df = df[df["item_id"] == item_id].reset_index(drop=True)
                if item_df.empty:
                    continue

                x_values = item_df[x_column].to_numpy(dtype=float)
                y_values = item_df[y_column].to_numpy(dtype=float)

                if x_mode == "shared":
                    offsets = (item_index - (len(ITEM_IDS_TO_PLOT) - 1) / 2.0) * bar_width
                    x_positions = x_values + offsets
                else:
                    x_positions = x_values

                ax.bar(
                    x_positions,
                    y_values,
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=style["alpha"],
                    edgecolor=style["edgecolor"],
                    linewidth=style["linewidth"],
                    label=ITEM_ID_NAME_MAP.get(item_id, f"item {item_id}"),
                )
        else:
            for item_id in ITEM_IDS_TO_PLOT:
                item_df = df[df["item_id"] == item_id].reset_index(drop=True)
                if item_df.empty:
                    continue

                x_values = item_df[x_column]
                y_values = item_df[y_column]

                ax.plot(
                    x_values,
                    y_values,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    label=ITEM_ID_NAME_MAP.get(item_id, f"item {item_id}"),
                    **style,
                )

        ax.set_title(PLOT_NAME_MAP[plot_key])
        ax.set_xlabel(X_AXIS_NAME_MAP[X_AXIS_KEY])
        ax.set_ylabel(Y_AXIS_NAME_MAP[plot_key])
        ax.grid(True, alpha=0.25)
        ax.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.98)
    fig.text(
        0.5,
        0.945,
        figure_description,
        ha="center",
        va="top",
        fontsize=10,
        wrap=True,
    )
    fig.subplots_adjust(top=0.84, hspace=0.35, wspace=0.25)

    plt.show()

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_svg, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()