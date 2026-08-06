################
#### CONFIG ####
################

# ---- LABELS / DATABASES ----
from pathlib import Path

# Map: label -> database path containing this label.
# Up to 4 labels can be plotted simultaneously.
LABEL_DB_PATH_MAP = {
    "barcode01_agilent_alignment_decoding": "/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/barcode01_agilent.db",
    "barcode02_dynegene_alignment_decoding": "/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/barcode02_dynegene.db",
    "barcode03_genscript_alignment_decoding": "/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/barcode03_genscript.db",
    "barcode04_six_images_alignment_decoding": "/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/barcode04_six_images.db",
}


LABEL_ITEM_IDS_TO_PLOT_MAP = {
    "barcode01_agilent_alignment_decoding": [0, ],
    "barcode02_dynegene_alignment_decoding": [0, ],
    "barcode03_genscript_alignment_decoding": [0, ],
    "barcode04_six_images_alignment_decoding": [3, ],
}

# pass_index | coverage | n_tot_reads | estimated_sequencing_duration | run_duration
X_AXIS_KEY = "coverage"

# Metrics to average and display.
METRICS_TO_PLOT = [
    "psnr",
    "hamming_distance_normalized",
    "perfectly_decoded_payload_ratio",
    "first_time_perfectly_decoded_payload_count",
]

# Interpolation grid settings for the direct averaging approach.
# Grid domain is [COMMON_X_ORIGIN, mean(end_x across runs)] for each item/metric.
COMMON_X_ORIGIN = 0.0
GRID_N_POINTS = 220

# Data quality / display toggles.
FILTER_OUT_NULL_PSNR_ROWS = True
PLOT_INDIVIDUAL_CURVES = False
PLOT_AVERAGE_CURVE = True
PLOT_SE_BAND = True
PLOT_STD_BAND = False

MIN_POINTS_PER_CURVE = 2

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
JOIN decoding_run_label_record l
    ON l.exp_id = m.exp_id
   AND l.dec_run_id = m.dec_run_id
WHERE l.label = ?
    AND m.item_id IN ({item_placeholders})
ORDER BY m.dec_run_id ASC, m.item_id ASC, m.dec_pass_id ASC
"""

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*) AS n_runs_total
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""

# ITEM_ID_NAME_MAP = {
#     0: "JPEGDNA-reference",
#     1: "JPEGDNA-delta-G",
#     2: "Motif-paircode",
# }


ITEM_ID_NAME_MAP = {
    0: "JPEGDNA-ref-bird",
    1: "Woman",
    2: "Burger",
    3: "JPEGDNA-ref-bird",
    4: "Night",
    5: "Day"
}


# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/AVG_decoding_run_passes.png"
OUTPUT_PATH_PDF = "../plots/AVG_decoding_run_passes.pdf"

# ---- NAMES ----

FIGURE_TITLE = "Average pass-level metrics over decoding runs"

FIGURE_DESCRIPTION_BASE = (
    "Average curves over runs, using interpolation on a shared X grid. "
    "Grid domain per item/metric is [0, mean of run-specific curve end X]."
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
    # Base color + 4 variants per item to distinguish up to 5 runs plotted together.
    0: ["#1f77b4", "#55bbff", "#3b93cc", "#61aad9", "#8ac2e6"],
    1: ["#d62728", "#c93a3b", "#bb4d4e", "#ae5f61", "#a17274"],
    2: ["#2ca02c", "#41aa41", "#56b456", "#6bbe6b", "#80c880"],
    3: ["#ff7f0e", "#ff9132", "#ffa457", "#ffb87d", "#ffcca3"],
    4: ["#9467bd", "#a27cc8", "#b091d3", "#bea7de", "#ccbee9"],
    5: ["#17becf", "#39c8d6", "#5bd2dd", "#7ddce4", "#9fe7ec"],
}


PLOT_STYLE_MAP = {
    "individual_curve": {
        "linestyle": "-",
        "linewidth": 1.1,
        "alpha": 0.35,
    },
    "average_curve": {
        "linestyle": "-",
        "linewidth": 2.5,
        "alpha": 0.98,
    },
    "se_band": {
        "alpha": 0.18,
    },
    "std_band": {
        "alpha": 0.10,
    },
}


##############
#### CODE ####
##############

import sqlite3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure


X_AXIS_COLUMN_MAP = {
    "pass_index": "dec_pass_id",
    "coverage": "coverage",
    "n_tot_reads": "N_trimmed_reads_cum_all_items",
    "estimated_sequencing_duration": "estimated_sequencing_duration",
    "run_duration": "run_duration",
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


def _get_x_column() -> str:
    if X_AXIS_KEY not in X_AXIS_COLUMN_MAP:
        raise ValueError(f"Unsupported X_AXIS_KEY: {X_AXIS_KEY}")
    return X_AXIS_COLUMN_MAP[X_AXIS_KEY]


def _get_item_base_color(item_id: int) -> str:
    variants = ITEM_ID_COLOR_MAP.get(item_id, ["#7f7f7f"])
    return variants[0] if variants else "#7f7f7f"


def _get_label_variant_color(item_id: int, label_index: int) -> str:
    variants = ITEM_ID_COLOR_MAP.get(item_id, ["#7f7f7f"])
    if not variants:
        return "#7f7f7f"
    return variants[label_index % len(variants)]


def _prepare_xy_curve(run_item_df: pd.DataFrame, x_col: str, y_col: str) -> tuple[np.ndarray, np.ndarray] | None:
    curve = run_item_df[[x_col, y_col]].copy()
    curve = curve.dropna(subset=[x_col, y_col])
    if curve.empty:
        return None

    # Merge duplicate x-values within one run by averaging y.
    curve = curve.groupby(x_col, as_index=False)[y_col].mean()
    order = np.argsort(curve[x_col].to_numpy(dtype=float))
    curve = curve.iloc[order].reset_index(drop=True)

    if len(curve) < MIN_POINTS_PER_CURVE:
        return None

    x_values = curve[x_col].to_numpy(dtype=float)
    y_values = curve[y_col].to_numpy(dtype=float)
    return x_values, y_values


def _interpolate_on_grid(
    x_values: np.ndarray,
    y_values: np.ndarray,
    x_grid: np.ndarray,
) -> np.ndarray:
    y_grid = np.full_like(x_grid, np.nan, dtype=float)
    in_support = (x_grid >= x_values.min()) & (x_grid <= x_values.max())
    if np.any(in_support):
        y_grid[in_support] = np.interp(x_grid[in_support], x_values, y_values)
    return y_grid


def _make_figure_axes(n_metrics: int) -> tuple[Figure, np.ndarray]:
    if n_metrics <= 2:
        fig, axes = plt.subplots(1, n_metrics, figsize=(7.5 * n_metrics, 4.8))
        axes_flat = np.array(axes).reshape(-1)
        return fig, axes_flat

    ncols = 2
    nrows = int(np.ceil(n_metrics / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5 * nrows))
    axes_flat = np.array(axes).reshape(-1)
    return fig, axes_flat


def main() -> None:
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_PDF))

    if not LABEL_DB_PATH_MAP:
        raise ValueError("LABEL_DB_PATH_MAP must contain at least one label -> database mapping.")
    if not LABEL_ITEM_IDS_TO_PLOT_MAP:
        raise ValueError("LABEL_ITEM_IDS_TO_PLOT_MAP must contain at least one label -> item ids mapping.")
    if len(LABEL_DB_PATH_MAP) > 4:
        raise ValueError("At most 4 labels are supported simultaneously.")
    if not METRICS_TO_PLOT:
        raise ValueError("METRICS_TO_PLOT must contain at least one metric key.")
    if GRID_N_POINTS < 2:
        raise ValueError("GRID_N_POINTS must be >= 2.")

    unknown_metrics = [m for m in METRICS_TO_PLOT if m not in METRIC_COLUMN_MAP]
    if unknown_metrics:
        raise ValueError(f"Unsupported metrics in METRICS_TO_PLOT: {unknown_metrics}")

    missing_labels_in_item_map = [
        label_name for label_name in LABEL_DB_PATH_MAP if label_name not in LABEL_ITEM_IDS_TO_PLOT_MAP
    ]
    if missing_labels_in_item_map:
        raise ValueError(
            "Each label in LABEL_DB_PATH_MAP must be present in LABEL_ITEM_IDS_TO_PLOT_MAP. "
            f"Missing: {missing_labels_in_item_map}"
        )

    label_item_ids_map: dict[str, list[int]] = {}
    for label_name in LABEL_DB_PATH_MAP:
        raw_item_ids = LABEL_ITEM_IDS_TO_PLOT_MAP.get(label_name, [])
        unique_item_ids = sorted({int(item_id) for item_id in raw_item_ids})
        if not unique_item_ids:
            raise ValueError(
                f"LABEL_ITEM_IDS_TO_PLOT_MAP['{label_name}'] must contain at least one item id."
            )
        label_item_ids_map[label_name] = unique_item_ids

    selected_item_ids_union = sorted(
        {item_id for item_ids in label_item_ids_map.values() for item_id in item_ids}
    )
    unknown_item_ids = sorted(set(selected_item_ids_union) - set(ITEM_ID_NAME_MAP.keys()))
    if unknown_item_ids:
        raise ValueError(
            "All selected item ids must be declared in ITEM_ID_NAME_MAP. "
            f"Unknown item ids: {unknown_item_ids}"
        )

    x_column = _get_x_column()

    frames: list[pd.DataFrame] = []
    runs_total_by_label: dict[str, int] = {}

    for label_name, relative_db_path in LABEL_DB_PATH_MAP.items():
        db_path = _resolve_path_from_script(Path(relative_db_path))
        label_item_ids = label_item_ids_map[label_name]
        item_placeholders = ",".join(["?"] * len(label_item_ids))
        query = SQL_QUERY.format(item_placeholders=item_placeholders)
        query_params = [label_name, *label_item_ids]
        with sqlite3.connect(db_path) as conn:
            label_df = pd.read_sql_query(query, conn, params=query_params)
            n_runs_total = int(conn.execute(SQL_QUERY_LABEL_RUNS_TOTAL, (label_name,)).fetchone()[0])
        runs_total_by_label[label_name] = n_runs_total

        if label_df.empty:
            print(f"[WARN] No rows for label '{label_name}' in DB '{db_path}'.")
            continue

        label_df = label_df.copy()
        label_df["label_name"] = label_name
        frames.append(label_df)

    if not frames:
        raise ValueError("No data returned for any label in LABEL_DB_PATH_MAP.")

    df = pd.concat(frames, ignore_index=True)

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()
    if df.empty:
        raise ValueError("No rows left after filtering ITEM_IDS_TO_PLOT to known item ids.")

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

    df = df.sort_values(["label_name", "dec_run_id", "item_id", "dec_pass_id"]).reset_index(drop=True)

    configured_labels_ordered = list(LABEL_DB_PATH_MAP.keys())
    displayed_labels_ordered = [
        label_name
        for label_name in configured_labels_ordered
        if label_name in set(df["label_name"].dropna().unique().tolist())
    ]
    label_order = {label_name: idx for idx, label_name in enumerate(displayed_labels_ordered)}

    displayed_item_ids_by_label: dict[str, list[int]] = {}
    for label_name in displayed_labels_ordered:
        configured_item_ids = label_item_ids_map[label_name]
        displayed_item_ids_by_label[label_name] = [
            item_id
            for item_id in configured_item_ids
            if ((df["label_name"] == label_name) & (df["item_id"] == item_id)).any()
        ]

    displayed_runs_by_label: dict[str, int] = {}
    for label_name in displayed_labels_ordered:
        label_runs = df[df["label_name"] == label_name]["dec_run_id"].nunique()
        displayed_runs_by_label[label_name] = int(label_runs)

    runs_summary = ", ".join(
        [
            f"{label_name}: studied={runs_total_by_label.get(label_name, 0)}, displayed={displayed_runs_by_label.get(label_name, 0)}"
            for label_name in displayed_labels_ordered
        ]
    )
    figure_description = (
        f"{FIGURE_DESCRIPTION_BASE} "
        f"Labels displayed={len(displayed_labels_ordered)}. {runs_summary}. "
        f"X grid origin={COMMON_X_ORIGIN}, grid points={GRID_N_POINTS}."
    )

    fig, axes_flat = _make_figure_axes(len(METRICS_TO_PLOT))

    for metric_idx, metric_key in enumerate(METRICS_TO_PLOT):
        ax = axes_flat[metric_idx]
        y_column = METRIC_COLUMN_MAP[metric_key]

        for label_name in displayed_labels_ordered:
            for item_id in displayed_item_ids_by_label.get(label_name, []):
                label_item_df = df[
                    (df["label_name"] == label_name) & (df["item_id"] == item_id)
                ].copy()
                if label_item_df.empty:
                    continue

                curves: list[tuple[str, np.ndarray, np.ndarray]] = []
                label_run_ids = sorted(label_item_df["dec_run_id"].dropna().unique().tolist())
                for run_id in label_run_ids:
                    run_item_df = label_item_df[label_item_df["dec_run_id"] == run_id]
                    prepared = _prepare_xy_curve(run_item_df, x_column, y_column)
                    if prepared is None:
                        continue
                    x_values, y_values = prepared
                    curves.append((run_id, x_values, y_values))

                if not curves:
                    continue

                end_x_values = np.array([curve[1].max() for curve in curves], dtype=float)
                x_grid_end = float(np.mean(end_x_values))
                if x_grid_end <= COMMON_X_ORIGIN:
                    continue

                x_grid = np.linspace(COMMON_X_ORIGIN, x_grid_end, GRID_N_POINTS)

                interpolated_rows = []
                label_index = label_order[label_name]
                label_color = _get_label_variant_color(item_id, label_index)
                for run_id, x_values, y_values in curves:
                    y_grid = _interpolate_on_grid(x_values, y_values, x_grid)
                    interpolated_rows.append(y_grid)

                    if PLOT_INDIVIDUAL_CURVES:
                        line_style = PLOT_STYLE_MAP["individual_curve"]
                        ax.plot(
                            x_grid,
                            y_grid,
                            color=label_color,
                            linestyle=line_style["linestyle"],
                            linewidth=line_style["linewidth"],
                            alpha=line_style["alpha"],
                        )

                y_matrix = np.vstack(interpolated_rows)
                n_contrib = np.sum(~np.isnan(y_matrix), axis=0)
                y_mean = np.nanmean(y_matrix, axis=0)

                y_se = np.full_like(y_mean, np.nan, dtype=float)
                y_std = np.full_like(y_mean, np.nan, dtype=float)
                valid_for_se = n_contrib > 1
                if np.any(valid_for_se):
                    std_ddof = np.nanstd(y_matrix[:, valid_for_se], axis=0, ddof=1)
                    y_std[valid_for_se] = std_ddof
                    y_se[valid_for_se] = std_ddof / np.sqrt(n_contrib[valid_for_se])

                if PLOT_STD_BAND:
                    std_band_style = PLOT_STYLE_MAP["std_band"]
                    ax.fill_between(
                        x_grid,
                        y_mean - y_std,
                        y_mean + y_std,
                        color=label_color,
                        alpha=std_band_style["alpha"],
                    )

                if PLOT_SE_BAND:
                    band_style = PLOT_STYLE_MAP["se_band"]
                    ax.fill_between(
                        x_grid,
                        y_mean - y_se,
                        y_mean + y_se,
                        color=label_color,
                        alpha=band_style["alpha"],
                    )

                if PLOT_AVERAGE_CURVE:
                    avg_style = PLOT_STYLE_MAP["average_curve"]
                    ax.plot(
                        x_grid,
                        y_mean,
                        color=label_color,
                        linestyle=avg_style["linestyle"],
                        linewidth=avg_style["linewidth"],
                        alpha=avg_style["alpha"],
                        label=(
                            f"{ITEM_ID_NAME_MAP.get(item_id, f'item {item_id}')} "
                            f"| {label_name} AVG (n_runs={len(curves)})"
                        ),
                    )

        ax.set_title(PLOT_NAME_MAP[metric_key])
        ax.set_xlabel(X_AXIS_NAME_MAP[X_AXIS_KEY])
        ax.set_ylabel(Y_AXIS_NAME_MAP[metric_key])
        ax.grid(True, alpha=0.25)
        ax.legend()

    for idx in range(len(METRICS_TO_PLOT), len(axes_flat)):
        axes_flat[idx].axis("off")

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
    output_path_pdf.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
