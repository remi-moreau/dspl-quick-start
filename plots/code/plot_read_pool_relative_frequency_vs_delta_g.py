################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

# "barcode01_agilent" | "barcode02_dynegene" | "barcode03_genscript" | "barcode04_six_images"
BARCODE = "barcode04_six_images"

if BARCODE == "barcode04_six_images":
    DB_PATH = f"/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/{BARCODE}.db"
    EXP_ID = "six_images"
    READ_POOL_ID = BARCODE
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
    DB_PATH = f"/media/remi-moreau/Seagate Expansion Drive/4_EXPERIENCES_PRO_STAGES/2026_Stage_3A_CNRS_I3S_MEDIACODING/5_DATA/2026-08_dspl_databases/{BARCODE}.db"
    EXP_ID = "synthesis_bench"
    READ_POOL_ID = BARCODE
    ITEM_IDS_TO_PLOT = [
        0,
        1,
        2,
    ]

    ITEM_ID_NAME_MAP = {
        0: "JPEGDNA-reference",
        1: "JPEGDNA-delta-G",
        2: "motif-paircode",
    }

# ---- PLOT CONFIG ----
ITEM_DISTINCTION = False

# If set to an integer N, only reads with position_in_read_pool < N are used.
# If None, no read-count limit is applied.
MAX_READ_COUNT = None

# ITEM_DISTINCTION behavior:
# - True: per-item mode. Relative frequency denominator is total reads for the item.
# - False: global mode. Relative frequency denominator is total reads in the selected
#   read-pool window (MAX_READ_COUNT), regardless of item.
# In both modes, ITEM_IDS_TO_PLOT is an inclusion filter for references.

DELTA_G_PRECISION = 1.0
MIN_POINTS_PER_BIN = 10

DISPLAY_LIN_REG = True

# ---- PERFORMANCE ----

SQLITE_CACHE_SIZE_KIB = 262144
SQLITE_TEMP_STORE = "MEMORY"
FORCE_READ_POOL_INDEX = True

READ_POOL_INDEX_NAME = "idx_read_exp_pool_fastq_pos"

SQL_QUERY_REF_COUNTS_AND_STATS = """
WITH selected_aligned_reads AS (
    SELECT
        r.exp_id,
        r.read_pool_id,
        r.enc_run_id,
        r.item_id,
        r.region_id,
        r.position_id
    FROM read r {read_index_hint}
    WHERE r.exp_id = ?
      AND r.read_pool_id = ?
      {max_read_where}
      AND r.item_id IN ({item_placeholders})
      AND r.enc_run_id IS NOT NULL
      AND r.region_id IS NOT NULL
      AND r.position_id IS NOT NULL
),
ref_counts AS (
    SELECT
        a.exp_id,
        a.enc_run_id,
        a.item_id,
        a.region_id,
        a.position_id,
        COUNT(*) AS n_reads_on_ref
    FROM selected_aligned_reads a
    GROUP BY a.exp_id, a.enc_run_id, a.item_id, a.region_id, a.position_id
),
ref_counts_with_optional_delta_g AS (
    SELECT
        c.item_id,
        c.region_id,
        c.position_id,
        ep.delta_g,
        c.n_reads_on_ref
    FROM ref_counts c
    LEFT JOIN encoded_payload ep
      ON ep.exp_id = c.exp_id
     AND ep.enc_run_id = c.enc_run_id
     AND ep.item_id = c.item_id
     AND ep.region_id = c.region_id
     AND ep.position_id = c.position_id
),
selected_pool_stats AS (
    SELECT
        COUNT(*) AS n_reads_selected,
        COUNT(DISTINCT r.fastq_index) AS n_fastq_selected
    FROM read r {read_index_hint}
    WHERE r.exp_id = ?
      AND r.read_pool_id = ?
      {max_read_where}
),
full_pool_stats AS (
    SELECT
        COUNT(*) AS n_reads_total,
        COUNT(DISTINCT r.fastq_index) AS n_fastq_total
    FROM read r
    WHERE r.exp_id = ?
      AND r.read_pool_id = ?
)
SELECT
    g.item_id,
    g.region_id,
    g.position_id,
    g.delta_g,
    g.n_reads_on_ref,
    AVG(g.n_reads_on_ref) OVER () AS mean_reads_per_ref_global,
    AVG(g.n_reads_on_ref) OVER (PARTITION BY g.item_id) AS mean_reads_per_ref_item,
    sps.n_reads_selected,
    sps.n_fastq_selected,
    fps.n_reads_total,
    fps.n_fastq_total
FROM ref_counts_with_optional_delta_g g
CROSS JOIN selected_pool_stats sps
CROSS JOIN full_pool_stats fps
ORDER BY g.item_id ASC, g.region_id ASC, g.position_id ASC
"""

# ---- OUTPUT PATH ----

OUTPUT_PATH_FIG_DISTRIBUTION_PNG = f"../plots/{BARCODE}/normalized_read_pool_coverage_distribution.png"
OUTPUT_PATH_FIG_DISTRIBUTION_PDF = f"../plots/{BARCODE}/normalized_read_pool_coverage_distribution.pdf"

OUTPUT_PATH_FIG_SCATTER_PNG = f"../plots/{BARCODE}/normalized_read_pool_coverage_vs_delta_g.png"
OUTPUT_PATH_FIG_SCATTER_PDF = f"../plots/{BARCODE}/normalized_read_pool_coverage_vs_delta_g.pdf"

OUTPUT_PATH_FIG_BINNED_PNG = f"../plots/{BARCODE}/mean_normalized_read_pool_coverage_in_bins_vs_delta_g.png"
OUTPUT_PATH_FIG_BINNED_PDF = f"../plots/{BARCODE}/mean_normalized_read_pool_coverage_in_bins_vs_delta_g.pdf"

# ---- NAMES ----

FIGURE_DISTRIBUTION_TITLE = f"Reference count distribution over Normalized Read Pool Coverage ({BARCODE})"
FIGURE_SCATTER_TITLE = f"Normalized Read Pool Coverage against delta G ({BARCODE})"
FIGURE_BINNED_TITLE = f"Mean Normalized Read Pool Coverage in delta G bins ({BARCODE})"

FIGURE_DESCRIPTION_BASE = (
    f"Read pool {READ_POOL_ID} (exp_id={EXP_ID}): per-reference Normalized Read Pool Coverage "
    "= reads_on_reference / mean(reads_on_reference)."
)

PLOT_NAME_MAP = {
    "distribution": "Reference counts over Normalized Read Pool Coverage",
    "scatter": "Normalized Read Pool Coverage against delta G for each reference",
    "binned_mean": "Binned mean Normalized Read Pool Coverage over delta G intervals",
    "binned_global": "Global binned mean Normalized Read Pool Coverage over delta G intervals",
}

X_AXIS_NAME_MAP = {
    "normalized_read_pool_coverage": "Normalized Read Pool Coverage",
    "delta_g": "Delta G",
}

Y_AXIS_NAME_MAP = {
    "n_references": "Number of references",
    "normalized_read_pool_coverage": "Normalized Read Pool Coverage",
    "normalized_read_pool_coverage_binned": "Mean Normalized Read Pool Coverage in bin",
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
    "distribution": {
        "alpha": 0.7,
        "edgecolor": "black",
        "linewidth": 0.4,
    },
    "scatter": {
        "marker": "o",
        "s": 18,
        "alpha": 0.5,
        "edgecolors": "none",
    },
    "lin_reg": {
        "linestyle": "-",
        "linewidth": 3.0,
        "alpha": 0.95,
    },
    "binned": {
        "bar_width_ratio": 0.92,
        "alpha": 0.9,
        "edgecolor": "black",
        "linewidth": 0.4,
    },
}


##############
#### CODE ####
##############

import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def _open_read_connection(db_path: Path) -> sqlite3.Connection:
    # Read-only URI + memory temp store improve large analytical scans.
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.execute(f"PRAGMA cache_size = {-int(SQLITE_CACHE_SIZE_KIB)}")
    conn.execute(f"PRAGMA temp_store = {SQLITE_TEMP_STORE}")
    return conn


def _compute_spearman_stats(x_data: np.ndarray, y_data: np.ndarray) -> tuple[float | None, float | None, str | None]:
    if len(x_data) < 2 or len(y_data) < 2:
        return None, None, "insufficient points"
    if np.unique(x_data).size < 2 or np.unique(y_data).size < 2:
        return None, None, "constant input"

    result = spearmanr(x_data, y_data)
    rho_raw = getattr(result, "statistic", getattr(result, "correlation", None))
    p_raw = getattr(result, "pvalue", None)
    if rho_raw is None or p_raw is None:
        return None, None, "undefined"

    rho_value = float(rho_raw)
    p_value = float(p_raw)
    if np.isnan(rho_value) or np.isnan(p_value):
        return None, None, "undefined"

    return rho_value, p_value, None


def main() -> None:
    db_path = _resolve_path_from_script(Path(DB_PATH))
    output_path_fig_distribution_png = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_DISTRIBUTION_PNG))
    output_path_fig_distribution_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_DISTRIBUTION_PDF))
    output_path_fig_scatter_png = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_SCATTER_PNG))
    output_path_fig_scatter_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_SCATTER_PDF))
    output_path_fig_binned_png = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_BINNED_PNG))
    output_path_fig_binned_pdf = _resolve_path_from_script(Path(OUTPUT_PATH_FIG_BINNED_PDF))

    if not ITEM_IDS_TO_PLOT:
        raise ValueError("ITEM_IDS_TO_PLOT must contain at least one item id.")

    max_read_count_value = MAX_READ_COUNT
    if max_read_count_value is not None:
        max_read_count_value = int(max_read_count_value)
        if max_read_count_value <= 0:
            raise ValueError("MAX_READ_COUNT must be a positive integer or None.")

    item_placeholders = ",".join(["?"] * len(ITEM_IDS_TO_PLOT))
    read_index_hint = f"INDEXED BY {READ_POOL_INDEX_NAME}" if FORCE_READ_POOL_INDEX else ""
    max_read_where = "AND r.position_in_read_pool < ?" if max_read_count_value is not None else ""
    max_read_params = [max_read_count_value] if max_read_count_value is not None else []

    query = SQL_QUERY_REF_COUNTS_AND_STATS.format(
        item_placeholders=item_placeholders,
        read_index_hint=read_index_hint,
        max_read_where=max_read_where,
    )
    query_params = [
        EXP_ID,
        READ_POOL_ID,
        *max_read_params,
        *ITEM_IDS_TO_PLOT,
        EXP_ID,
        READ_POOL_ID,
        *max_read_params,
        EXP_ID,
        READ_POOL_ID,
    ]

    with _open_read_connection(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)

    if df.empty:
        raise ValueError(
            "No aligned reference data found. Check EXP_ID, READ_POOL_ID and ITEM_IDS_TO_PLOT."
        )

    n_reads_total = int(df["n_reads_total"].iloc[0]) if df["n_reads_total"].notna().any() else 0
    n_fastq = int(df["n_fastq_total"].iloc[0]) if df["n_fastq_total"].notna().any() else 0
    n_reads_selected = int(df["n_reads_selected"].iloc[0]) if df["n_reads_selected"].notna().any() else 0
    n_fastq_selected = int(df["n_fastq_selected"].iloc[0]) if df["n_fastq_selected"].notna().any() else 0
    mode_text = "per-item" if ITEM_DISTINCTION else "global"
    max_read_text = "None" if max_read_count_value is None else str(max_read_count_value)
    denominator_text = (
        "mean reads_on_reference computed per item"
        if ITEM_DISTINCTION
        else "mean reads_on_reference computed globally"
    )
    figure_description = (
        f"{FIGURE_DESCRIPTION_BASE} "
        f"\nMode={mode_text}; denominator={denominator_text}."
        f"\nRead filter: ITEM_IDS_TO_PLOT applied; MAX_READ_COUNT={max_read_text}."
        f"\nRead pool volumetry: reads_total={n_reads_total}, reads_selected={n_reads_selected}, "
        f"fastq_total={n_fastq}, fastq_selected={n_fastq_selected}."
    )

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()

    df["delta_g"] = pd.to_numeric(df["delta_g"], errors="coerce")
    df["n_reads_on_ref"] = pd.to_numeric(df["n_reads_on_ref"], errors="coerce")
    df["mean_reads_per_ref_global"] = pd.to_numeric(df["mean_reads_per_ref_global"], errors="coerce")
    df["mean_reads_per_ref_item"] = pd.to_numeric(df["mean_reads_per_ref_item"], errors="coerce")

    if ITEM_DISTINCTION:
        df["normalized_read_pool_coverage"] = (
            df["n_reads_on_ref"] / df["mean_reads_per_ref_item"]
        )
    else:
        df["normalized_read_pool_coverage"] = (
            df["n_reads_on_ref"] / df["mean_reads_per_ref_global"]
        )

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["normalized_read_pool_coverage"])

    item_ids_to_render = [item_id for item_id in ITEM_IDS_TO_PLOT if item_id in ITEM_ID_NAME_MAP]
    if not item_ids_to_render:
        raise ValueError("No configured item_id is present in ITEM_ID_NAME_MAP.")

    # -------------------------
    # Figure 1: Distribution
    # -------------------------
    fig_distribution, ax_distribution = plt.subplots(1, 1, figsize=(11.5, 5.2))
    distribution_style = PLOT_STYLE_MAP["distribution"]

    all_values = df["normalized_read_pool_coverage"].to_numpy(dtype=float)
    if all_values.size < 2:
        bins_for_hist: int | list[float] = 10
    else:
        bins_for_hist = np.histogram_bin_edges(all_values, bins="fd").tolist()

    if ITEM_DISTINCTION:
        for item_id in item_ids_to_render:
            item_df = df[df["item_id"] == item_id]
            if item_df.empty:
                continue
            ax_distribution.hist(
                item_df["normalized_read_pool_coverage"].to_numpy(dtype=float),
                bins=bins_for_hist,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                alpha=distribution_style["alpha"],
                edgecolor=distribution_style["edgecolor"],
                linewidth=distribution_style["linewidth"],
                label=f"{ITEM_ID_NAME_MAP.get(item_id, f'item_id={item_id}')} refs (n={len(item_df)})",
            )
        ax_distribution.legend()
    else:
        ax_distribution.hist(
            all_values,
            bins=bins_for_hist,
            color="#4c4c4c",
            alpha=distribution_style["alpha"],
            edgecolor=distribution_style["edgecolor"],
            linewidth=distribution_style["linewidth"],
            label=f"All selected refs (n={len(df)})",
        )
        ax_distribution.legend()

    ax_distribution.set_title(PLOT_NAME_MAP["distribution"])
    ax_distribution.set_xlabel(X_AXIS_NAME_MAP["normalized_read_pool_coverage"])
    ax_distribution.set_ylabel(Y_AXIS_NAME_MAP["n_references"])
    ax_distribution.grid(True, alpha=0.25, linestyle="--")

    fig_distribution.suptitle(FIGURE_DISTRIBUTION_TITLE, fontsize=15, y=0.98)
    fig_distribution.text(0.5, 0.94, figure_description, ha="center", va="top", wrap=True, fontsize=10)
    fig_distribution.tight_layout(rect=(0.03, 0.04, 0.97, 0.9))

    # Display this new figure first.
    plt.show()

    output_path_fig_distribution_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_fig_distribution_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig_distribution.savefig(output_path_fig_distribution_png, dpi=250)
    fig_distribution.savefig(output_path_fig_distribution_pdf)
    plt.close(fig_distribution)

    df_with_delta = df.dropna(subset=["delta_g"]).copy()
    if df_with_delta.empty:
        print(
            "No delta_g data available for selected references. "
            "Generated only Figure 1 (distribution)."
        )
        return

    df = df_with_delta

    # -------------------------
    # Figure 2: Scatter vs delta_g
    # -------------------------
    fig_scatter, scatter_ax = plt.subplots(1, 1, figsize=(11.5, 5.5))
    scatter_style = PLOT_STYLE_MAP["scatter"]
    lin_reg_style = PLOT_STYLE_MAP["lin_reg"]
    spearman_lines: list[str] = []

    for item_id in item_ids_to_render:
        item_name = ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}")
        item_df = df[df["item_id"] == item_id].copy()
        if item_df.empty:
            continue

        x_data = item_df["delta_g"].to_numpy(dtype=float)
        y_data = item_df["normalized_read_pool_coverage"].to_numpy(dtype=float)

        scatter_ax.scatter(
            x_data,
            y_data,
            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
            label=f"{item_name} refs (n={len(item_df)})",
            marker=scatter_style["marker"],
            s=scatter_style["s"],
            alpha=scatter_style["alpha"],
            edgecolors=scatter_style["edgecolors"],
        )

        if ITEM_DISTINCTION:
            rho_value, p_value, spearman_issue = _compute_spearman_stats(x_data, y_data)
            if spearman_issue is None:
                line_stdout = (
                    f"Spearman {item_name}: n={len(item_df)}, "
                    f"rho={rho_value:.6f}, p-value={p_value:.3e}"
                )
                line_figure = f"{item_name}: rho={rho_value:.3f}, p={p_value:.3e}, n={len(item_df)}"
            else:
                line_stdout = f"Spearman {item_name}: n={len(item_df)}, undefined ({spearman_issue})"
                line_figure = f"{item_name}: undefined ({spearman_issue}), n={len(item_df)}"
            print(line_stdout)
            spearman_lines.append(line_figure)

        if ITEM_DISTINCTION and DISPLAY_LIN_REG and len(item_df) >= 2 and np.std(x_data) > 0:
            slope, intercept = np.polyfit(x_data, y_data, deg=1)
            y_pred = (slope * x_data) + intercept
            ss_res = float(np.sum((y_data - y_pred) ** 2))
            ss_tot = float(np.sum((y_data - float(np.mean(y_data))) ** 2))
            r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

            x_line = np.linspace(float(np.min(x_data)), float(np.max(x_data)), 200)
            y_line = (slope * x_line) + intercept
            scatter_ax.plot(
                x_line,
                y_line,
                color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                linestyle=lin_reg_style["linestyle"],
                linewidth=lin_reg_style["linewidth"],
                alpha=lin_reg_style["alpha"],
                label=f"{item_name} lin reg: coef={slope:.3e}, R²={r_squared:.3e}",
            )

    if (not ITEM_DISTINCTION) and DISPLAY_LIN_REG and len(df) >= 2:
        x_global = df["delta_g"].to_numpy(dtype=float)
        y_global = df["normalized_read_pool_coverage"].to_numpy(dtype=float)
        if np.std(x_global) > 0:
            slope, intercept = np.polyfit(x_global, y_global, deg=1)
            y_pred = (slope * x_global) + intercept
            ss_res = float(np.sum((y_global - y_pred) ** 2))
            ss_tot = float(np.sum((y_global - float(np.mean(y_global))) ** 2))
            r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

            x_line = np.linspace(float(np.min(x_global)), float(np.max(x_global)), 200)
            y_line = (slope * x_line) + intercept
            scatter_ax.plot(
                x_line,
                y_line,
                color="black",
                linestyle=lin_reg_style["linestyle"],
                linewidth=lin_reg_style["linewidth"],
                alpha=lin_reg_style["alpha"],
                label=f"Global lin reg: coef={slope:.3e}, R²={r_squared:.3e}",
            )

    if not ITEM_DISTINCTION:
        x_global = df["delta_g"].to_numpy(dtype=float)
        y_global = df["normalized_read_pool_coverage"].to_numpy(dtype=float)
        rho_value, p_value, spearman_issue = _compute_spearman_stats(x_global, y_global)
        if spearman_issue is None:
            print(
                f"Spearman global: n={len(df)}, rho={rho_value:.6f}, p-value={p_value:.3e}"
            )
            spearman_lines.append(
                f"Global: rho={rho_value:.3f}, p={p_value:.3e}, n={len(df)}"
            )
        else:
            print(f"Spearman global: n={len(df)}, undefined ({spearman_issue})")
            spearman_lines.append(f"Global: undefined ({spearman_issue}), n={len(df)}")

    spearman_block_text = None
    if spearman_lines:
        spearman_block_text = "Spearman correlation\n" + "\n".join(spearman_lines)

    scatter_ax.set_title(PLOT_NAME_MAP["scatter"])
    scatter_ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
    scatter_ax.set_ylabel(Y_AXIS_NAME_MAP["normalized_read_pool_coverage"])
    scatter_ax.grid(True, alpha=0.25, linestyle="--")
    scatter_ax.legend()

    fig_scatter.suptitle(FIGURE_SCATTER_TITLE, fontsize=15, y=0.98)
    fig_scatter.text(0.5, 0.94, figure_description, ha="center", va="top", wrap=True, fontsize=10)
    if spearman_block_text is not None:
        fig_scatter.text(
            0.015,
            0.89,
            spearman_block_text,
            ha="left",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none", "pad": 4.0},
        )
    layout_top = 0.80 if spearman_block_text is not None else 0.9
    fig_scatter.tight_layout(rect=(0.03, 0.04, 0.97, layout_top))

    plt.show()

    output_path_fig_scatter_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_fig_scatter_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig_scatter.savefig(output_path_fig_scatter_png, dpi=250)
    fig_scatter.savefig(output_path_fig_scatter_pdf)
    plt.close(fig_scatter)

    # -------------------------
    # Figure 3: Binned means
    # -------------------------
    n_rows_binned = len(item_ids_to_render) if ITEM_DISTINCTION else 1
    fig_binned, axes_binned = plt.subplots(
        n_rows_binned,
        1,
        figsize=(11.5, 4.4 * n_rows_binned),
        sharex=False,
    )
    if n_rows_binned == 1:
        axes_binned = [axes_binned]

    binned_style = PLOT_STYLE_MAP["binned"]

    if ITEM_DISTINCTION:
        for ax, item_id in zip(axes_binned, item_ids_to_render):
            item_name = ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}")
            item_df = df[df["item_id"] == item_id].copy()

            if item_df.empty:
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{PLOT_NAME_MAP['binned_mean']} - {item_name}")
                ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
                ax.set_ylabel(Y_AXIS_NAME_MAP["relative_frequency_binned"])
                ax.grid(True, alpha=0.25, linestyle="--")
                continue

            item_df["bin_upper"] = np.ceil(item_df["delta_g"] / DELTA_G_PRECISION) * DELTA_G_PRECISION
            item_binned = (
                item_df.groupby("bin_upper", as_index=False)
                .agg(
                    n_points=("normalized_read_pool_coverage", "size"),
                    mean_normalized_coverage=("normalized_read_pool_coverage", "mean"),
                )
            )
            item_binned = item_binned[item_binned["n_points"] >= MIN_POINTS_PER_BIN].copy()

            if item_binned.empty:
                ax.text(
                    0.5,
                    0.5,
                    "No bins pass threshold",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
            else:
                item_binned["bin_center"] = item_binned["bin_upper"] - (DELTA_G_PRECISION / 2.0)
                bar_width = DELTA_G_PRECISION * float(binned_style["bar_width_ratio"])
                ax.bar(
                    item_binned["bin_center"].to_numpy(dtype=float),
                    item_binned["mean_normalized_coverage"].to_numpy(dtype=float),
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=binned_style["alpha"],
                    edgecolor=binned_style["edgecolor"],
                    linewidth=binned_style["linewidth"],
                    label=f"{item_name} bins>= {MIN_POINTS_PER_BIN}",
                    align="center",
                )
                ax.legend()

            ax.set_title(f"{PLOT_NAME_MAP['binned_mean']} - {item_name}")
            ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
            ax.set_ylabel(Y_AXIS_NAME_MAP["normalized_read_pool_coverage_binned"])
            ax.grid(True, alpha=0.25, linestyle="--")
    else:
        binned_ax = axes_binned[0]
        global_df = df.copy()
        global_df["bin_upper"] = np.ceil(global_df["delta_g"] / DELTA_G_PRECISION) * DELTA_G_PRECISION
        global_binned = (
            global_df.groupby("bin_upper", as_index=False)
            .agg(
                n_points=("normalized_read_pool_coverage", "size"),
                mean_normalized_coverage=("normalized_read_pool_coverage", "mean"),
            )
        )
        global_binned = global_binned[global_binned["n_points"] >= MIN_POINTS_PER_BIN].copy()

        if global_binned.empty:
            binned_ax.text(
                0.5,
                0.5,
                "No bins pass threshold",
                ha="center",
                va="center",
                transform=binned_ax.transAxes,
            )
        else:
            global_binned["bin_center"] = global_binned["bin_upper"] - (DELTA_G_PRECISION / 2.0)
            bar_width = DELTA_G_PRECISION * float(binned_style["bar_width_ratio"])
            binned_ax.bar(
                global_binned["bin_center"].to_numpy(dtype=float),
                global_binned["mean_normalized_coverage"].to_numpy(dtype=float),
                width=bar_width,
                color="#4c4c4c",
                alpha=binned_style["alpha"],
                edgecolor=binned_style["edgecolor"],
                linewidth=binned_style["linewidth"],
                label=f"Global bins>= {MIN_POINTS_PER_BIN}",
                align="center",
            )
            binned_ax.legend()

        binned_ax.set_title(PLOT_NAME_MAP["binned_global"])
        binned_ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
        binned_ax.set_ylabel(Y_AXIS_NAME_MAP["normalized_read_pool_coverage_binned"])
        binned_ax.grid(True, alpha=0.25, linestyle="--")

    fig_binned.suptitle(FIGURE_BINNED_TITLE, fontsize=15, y=0.985)
    fig_binned.text(0.5, 0.955, figure_description, ha="center", va="top", wrap=True, fontsize=10)
    fig_binned.tight_layout(rect=(0.03, 0.04, 0.97, 0.92))

    plt.show()

    output_path_fig_binned_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_fig_binned_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig_binned.savefig(output_path_fig_binned_png, dpi=250)
    fig_binned.savefig(output_path_fig_binned_pdf)
    plt.close(fig_binned)


if __name__ == "__main__":
    main()
