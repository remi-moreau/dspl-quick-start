
################
#### CONFIG ####
################

# ---- DATA PATH ----
from pathlib import Path

DATA_PATH_PREFIX = "../../database/exports/"

DATA_PATH_SUFFIX = "barcode01_agilent/metrics_view_run_item/at_decoding.csv"

DATA_PATH = Path(DATA_PATH_PREFIX) / Path(DATA_PATH_SUFFIX)

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}


# ---- DATA STRUCTURE ----

EXPECTED_COLUMN_NAMES_IN_THE_CSV = [
    "item_id",
    "dec_run_id",
    "estimated_sequencing_duration_at_decoding",
    "coverage_at_decoding",
    "hamming_distance_normalized_at_decoding",
    "perfectly_decoded_payload_ratio_at_decoding",
]


# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/values_at_decoding_over_sequencing.png"
OUTPUT_PATH_SVG = "../plots/values_at_decoding_over_sequencing.svg"

# ---- NAMES ----

FIGURE_TITLE = "Metrics at decoding over sequencing"

FIGURE_DESCRIPTION = "Different metrics at decoding over sequencing. "\
    + "The X axis is the number of the considered run, which increases"\
    + " with sequencing time."

PLOT_NAME_MAP = {
    "coverage": "Average coverage at decoding over sequencing",
    "hamming_distance_normalized_at_decoding": "Normalized Hamming distance at decoding", 
    "perfectly_decoded_payload_ratio_at_decoding": "Perfectly decoded payload ratio at decoding",
    "estimated_sequencing_duration_at_decoding": "Estimated sequencing duration of the decoding run",
}

X_AXIS_NAME_MAP = {
    "run_number": "Run",
}

Y_AXIS_NAME_MAP = {
    "coverage": "Norma. avg. coverage",
    "hamming_distance_normalized_at_decoding": "Norm. Hamming distance", 
    "perfectly_decoded_payload_ratio_at_decoding": "Perfectly decoded payload ratio (%)",
    "estimated_sequencing_duration_at_decoding": "Est. sequencing duration (minutes)",
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
        "linewidth": 2.2,
        "marker": "o",
        "markersize": 6,
    },
    "hamming_distance_normalized_at_decoding": {
        "linestyle": "--",
        "linewidth": 2.0,
        "marker": "s",
        "markersize": 5.5,
    },
    "perfectly_decoded_payload_ratio_at_decoding": {
        "linestyle": "-.",
        "linewidth": 2.0,
        "marker": "^",
        "markersize": 6,
    },
    "estimated_sequencing_duration_at_decoding": {
        "linestyle": ":",
        "linewidth": 2.2,
        "marker": "D",
        "markersize": 5.5,
    },
    "coverage_mean_marker": {
        "linestyle": "--",
        "linewidth": 1.8,
        "alpha": 0.9,
    },
}


##############
#### CODE ####
##############

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
    data_path = (base_dir / DATA_PATH).resolve()
    output_path_png = (base_dir / OUTPUT_PATH_PNG).resolve()
    output_path_svg = (base_dir / OUTPUT_PATH_SVG).resolve()

    df = pd.read_csv(data_path)
    missing_columns = [
        column for column in EXPECTED_COLUMN_NAMES_IN_THE_CSV if column not in df.columns
    ]
    if missing_columns:
        raise ValueError(
            "Missing expected columns in CSV: " + ", ".join(sorted(missing_columns))
        )

    for column in METRIC_COLUMN_MAP.values():
        df[column] = pd.to_numeric(df[column], errors="coerce")

    average_coverage_by_item = {
        item_id: float(
            df[df["item_id"] == item_id]["coverage_at_decoding"].mean()
        )
        for item_id in ITEM_ID_NAME_MAP
        if not df[df["item_id"] == item_id]["coverage_at_decoding"].dropna().empty
    }

    # Convert ratio to percent and sequencing duration from seconds to minutes.
    df["perfectly_decoded_payload_ratio_at_decoding"] = (
        df["perfectly_decoded_payload_ratio_at_decoding"] * 100.0
    )
    df["estimated_sequencing_duration_at_decoding"] = (
        df["estimated_sequencing_duration_at_decoding"] / 60.0
    )

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

        if metric_key == "coverage":
            marker_style = PLOT_STYLE_MAP["coverage_mean_marker"]
            for item_id, avg_cov in average_coverage_by_item.items():
                ax.axhline(
                    y=avg_cov,
                    color=ITEM_ID_COLOR_MAP.get(item_id),
                    linestyle=marker_style["linestyle"],
                    linewidth=marker_style["linewidth"],
                    alpha=marker_style["alpha"],
                    label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item {item_id}')} mean={avg_cov:.2f}",
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
        FIGURE_DESCRIPTION,
        ha="center",
        va="top",
        fontsize=10,
        wrap=True,
    )
    fig.subplots_adjust(top=0.84, hspace=0.35, wspace=0.25)

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_svg, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()





