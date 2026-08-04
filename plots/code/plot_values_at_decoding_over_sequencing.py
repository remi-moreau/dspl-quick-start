
################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

DB_PATH = "../../database/barcode03_genscript.db"

RUN_LABEL = "barcode03_genscript_alignment_decoding"
#RUN_LABEL = "test_label"

ITEM_IDS_TO_PLOT = [0, 1]

DISPLAY_METRIC_MEANS = True

SQL_QUERY = """
SELECT
        m.item_id,
        m.dec_run_id,
        m.estimated_sequencing_duration_at_decoding,
        m.coverage_at_decoding,
        m.hamming_distance_normalized_at_decoding,
        m.perfectly_decoded_payload_ratio_at_decoding
FROM metrics_view_run_item_at_decoding m
JOIN decoding_run_label_record l
    ON l.exp_id = m.exp_id
 AND l.dec_run_id = m.dec_run_id
WHERE l.label = ?
    AND m.item_id IN ({item_placeholders})
    AND m.pass_at_decoding IS NOT NULL
ORDER BY m.item_id ASC, m.dec_run_id ASC
"""

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*) AS n_runs_total
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""

ITEM_ID_NAME_MAP = {
    0: "JPEGDNA",
    1: "JPEGDNA-delta-G",
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/values_at_decoding_over_sequencing.png"
OUTPUT_PATH_SVG = "../plots/values_at_decoding_over_sequencing.svg"
OUTPUT_PATH_PDF = "../plots/values_at_decoding_over_sequencing.pdf"

# ---- NAMES ----

FIGURE_TITLE = "Metrics at image decoding over sequencing (successive runs)"

FIGURE_DESCRIPTION_BASE = "Different image decoding-level metrics over sequencing. "\
    + "The X axis is the numbers of the successive runs. A new run starts as soon as the previous one succeed."\

PLOT_NAME_MAP = {
    "coverage": "Average coverage at image decoding",
    "hamming_distance_normalized_at_decoding": "Normalized Hamming distance at image decoding", 
    "perfectly_decoded_payload_ratio_at_decoding": "Perfectly decoded payload ratio at image decoding",
    "estimated_sequencing_duration_at_decoding": "Sequencing duration of the decoding run at image decoding",
}

X_AXIS_NAME_MAP = {
    "run_number": "Run",
}

Y_AXIS_NAME_MAP = {
    "coverage": "Normalized average coverage",
    "hamming_distance_normalized_at_decoding": "Normalized Hamming distance", 
    "perfectly_decoded_payload_ratio_at_decoding": "Perfectly decoded payload ratio (%)",
    "estimated_sequencing_duration_at_decoding": "Sequencing duration (minutes)",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green (fallback)
}

PLOT_STYLE_MAP = {
    "coverage": {
        "linestyle": "-",
        "linewidth": 1.8,
        "marker": "o",
        "markersize": 4,
    },
    "hamming_distance_normalized_at_decoding": {
        "linestyle": "-",
        "linewidth": 1.8,
        "marker": "o",
        "markersize": 4,
    },
    "perfectly_decoded_payload_ratio_at_decoding": {
        "linestyle": "-",
        "linewidth": 1.8,
        "marker": "o",
        "markersize": 4,
    },
    "estimated_sequencing_duration_at_decoding": {
        "linestyle": "-",
        "linewidth": 1.8,
        "marker": "o",
        "markersize": 4,
    },
    "metric_mean_marker": {
        "linestyle": "--",
        "linewidth": 1.8,
        "alpha": 0.9,
    },
}


##############
#### CODE ####
##############

import sqlite3
import matplotlib.pyplot as plt
import pandas as pd


METRIC_COLUMN_MAP = {
    "coverage": "coverage_at_decoding",
    "hamming_distance_normalized_at_decoding": "hamming_distance_normalized_at_decoding",
    "perfectly_decoded_payload_ratio_at_decoding": "perfectly_decoded_payload_ratio_at_decoding",
    "estimated_sequencing_duration_at_decoding": "estimated_sequencing_duration_at_decoding",
}


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    db_path = (base_dir / DB_PATH).resolve()
    output_path_png = (base_dir / OUTPUT_PATH_PNG).resolve()
    output_path_svg = (base_dir / OUTPUT_PATH_SVG).resolve()
    output_path_pdf = (base_dir / OUTPUT_PATH_PDF).resolve()

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

    if df.empty:
        raise ValueError("No data returned by SQL query. Check RUN_LABEL and ITEM_IDS_TO_PLOT.")

    n_runs_displayed = int(df["dec_run_id"].nunique())
    figure_description = (
        f"{FIGURE_DESCRIPTION_BASE} "
        f"Label '{RUN_LABEL}': decoding runs studied={n_runs_total}, displayed={n_runs_displayed}."
    )

    for column in METRIC_COLUMN_MAP.values():
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # Convert ratio to percent and sequencing duration from seconds to minutes.
    df["perfectly_decoded_payload_ratio_at_decoding"] = (
        df["perfectly_decoded_payload_ratio_at_decoding"] * 100.0
    )
    df["estimated_sequencing_duration_at_decoding"] = (
        df["estimated_sequencing_duration_at_decoding"] / 60.0
    )

    # Means are computed after unit conversions so marker values match displayed curves.
    average_by_metric_by_item: dict[str, dict[int, float]] = {}
    for metric_key, metric_column in METRIC_COLUMN_MAP.items():
        metric_means: dict[int, float] = {}
        for item_id in ITEM_ID_NAME_MAP:
            item_values = df[df["item_id"] == item_id][metric_column].dropna()
            if item_values.empty:
                continue
            metric_means[item_id] = float(item_values.mean())
        average_by_metric_by_item[metric_key] = metric_means

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))
    axes_flat = axes.flatten()
    metric_keys = list(PLOT_NAME_MAP.keys())

    for idx, metric_key in enumerate(metric_keys):
        ax = axes_flat[idx]
        metric_column = METRIC_COLUMN_MAP[metric_key]
        style = PLOT_STYLE_MAP.get(metric_key, {})

        for item_id in ITEM_ID_NAME_MAP:
            item_df = df[df["item_id"] == item_id].reset_index(drop=True)
            if item_df.empty:
                continue

            x_values = range(1, len(item_df) + 1)
            y_values = item_df[metric_column]

            ax.plot(
                x_values,
                y_values,
                color=ITEM_ID_COLOR_MAP.get(item_id),
                label=ITEM_ID_NAME_MAP.get(item_id, f"item {item_id}"),
                **style,
            )

        if DISPLAY_METRIC_MEANS:
            marker_style = PLOT_STYLE_MAP["metric_mean_marker"]
            for item_id, mean_value in average_by_metric_by_item.get(metric_key, {}).items():
                ax.axhline(
                    y=mean_value,
                    color=ITEM_ID_COLOR_MAP.get(item_id),
                    linestyle=marker_style["linestyle"],
                    linewidth=marker_style["linewidth"],
                    alpha=marker_style["alpha"],
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item {item_id}')} mean={mean_value:.2f}",
                )

        ax.set_title(PLOT_NAME_MAP[metric_key])
        ax.set_xlabel(X_AXIS_NAME_MAP["run_number"])
        ax.set_ylabel(Y_AXIS_NAME_MAP[metric_key])
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
    output_path_pdf.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_svg, bbox_inches="tight")
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()





