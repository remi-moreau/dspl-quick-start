################
#### CONFIG ####
################

# ---- DATABASE ----
from pathlib import Path

DB_PATH = "../../database/dspl.db"

EXP_ID = "synthesis_bench"
READ_POOL_ID = "barcode01_agilent"

ITEM_IDS_TO_PLOT = [0, 1]

DELTA_G_PRECISION = 1.0
MIN_POINTS_PER_BIN = 2

DISPLAY_LIN_REG = True

SQL_QUERY = """
WITH read_pool_scope AS (
    SELECT
            r.exp_id,
            r.read_pool_id,
            r.enc_run_id,
            r.item_id,
            r.region_id,
            r.position_id
    FROM read r
    WHERE r.exp_id = ?
      AND r.read_pool_id = ?
      AND r.item_id IN ({item_placeholders})
),
item_totals AS (
    SELECT
            exp_id,
            read_pool_id,
            item_id,
            COUNT(*) AS n_reads_item_total
    FROM read_pool_scope
    GROUP BY exp_id, read_pool_id, item_id
),
aligned_ref_counts AS (
    SELECT
            exp_id,
            read_pool_id,
            enc_run_id,
            item_id,
            region_id,
            position_id,
            COUNT(*) AS n_reads_on_ref
    FROM read_pool_scope
    WHERE enc_run_id IS NOT NULL
      AND region_id IS NOT NULL
      AND position_id IS NOT NULL
    GROUP BY exp_id, read_pool_id, enc_run_id, item_id, region_id, position_id
),
ref_with_delta_g AS (
    SELECT
            a.item_id,
            a.region_id,
            a.position_id,
            ep.delta_g,
            a.n_reads_on_ref,
            t.n_reads_item_total
    FROM aligned_ref_counts a
    JOIN item_totals t
      ON t.exp_id = a.exp_id
     AND t.read_pool_id = a.read_pool_id
     AND t.item_id = a.item_id
    JOIN encoded_payload ep
      ON ep.exp_id = a.exp_id
     AND ep.enc_run_id = a.enc_run_id
     AND ep.item_id = a.item_id
     AND ep.region_id = a.region_id
     AND ep.position_id = a.position_id
)
SELECT
        item_id,
        region_id,
        position_id,
        delta_g,
        n_reads_on_ref,
        n_reads_item_total,
        (1.0 * n_reads_on_ref) / NULLIF(n_reads_item_total, 0) AS relative_frequency
FROM ref_with_delta_g
WHERE delta_g IS NOT NULL
ORDER BY item_id ASC, region_id ASC, position_id ASC
"""

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/read_pool_relative_frequency_vs_delta_g.png"
OUTPUT_PATH_SVG = "../plots/read_pool_relative_frequency_vs_delta_g.svg"

# ---- NAMES ----

FIGURE_TITLE = "Relative reference frequency in read pool against delta G"

FIGURE_DESCRIPTION = (
    f"Read pool {READ_POOL_ID} (exp_id={EXP_ID}): per-reference relative frequency "
    "(occurrence count / total reads in pool for the item) against delta G."
)

PLOT_NAME_MAP = {
    "scatter": "Per-reference relative frequency against delta G",
    "binned_mean": "Binned mean relative frequency over delta G intervals",
}

X_AXIS_NAME_MAP = {
    "delta_g": "Delta G",
}

Y_AXIS_NAME_MAP = {
    "relative_frequency": "Relative frequency in read pool",
    "relative_frequency_binned": "Mean relative frequency in bin",
}


# ---- STYLES ----

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green (fallback)
}

PLOT_STYLE_MAP = {
    "scatter": {
        "marker": "o",
        "s": 18,
        "alpha": 0.5,
        "edgecolors": "none",
    },
    "lin_reg": {
        "linestyle": "-",
        "linewidth": 2.0,
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
    query_params = [EXP_ID, READ_POOL_ID, *ITEM_IDS_TO_PLOT]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=query_params)

    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()

    if df.empty:
        raise ValueError(
            "No aligned reference data with delta_g found. Check EXP_ID, READ_POOL_ID and ITEM_IDS_TO_PLOT."
        )

    df["delta_g"] = pd.to_numeric(df["delta_g"], errors="coerce")
    df["relative_frequency"] = pd.to_numeric(df["relative_frequency"], errors="coerce")
    df = df.dropna(subset=["delta_g", "relative_frequency"])

    item_ids_to_render = [item_id for item_id in ITEM_IDS_TO_PLOT if item_id in ITEM_ID_NAME_MAP]
    if not item_ids_to_render:
        raise ValueError("No configured item_id is present in ITEM_ID_NAME_MAP.")

    # Build one scatter + one binned subplot per item.
    n_rows = 1 + len(item_ids_to_render)
    fig, axes = plt.subplots(
        n_rows,
        1,
        figsize=(11.5, 4.4 * n_rows),
        sharex=False,
    )
    if n_rows == 1:
        axes = [axes]

    scatter_ax = axes[0]
    scatter_style = PLOT_STYLE_MAP["scatter"]
    lin_reg_style = PLOT_STYLE_MAP["lin_reg"]
    binned_style = PLOT_STYLE_MAP["binned"]

    for item_id in item_ids_to_render:
        item_name = ITEM_ID_NAME_MAP.get(item_id, f"item_id={item_id}")
        item_df = df[df["item_id"] == item_id].copy()
        if item_df.empty:
            continue

        x_data = item_df["delta_g"].to_numpy(dtype=float)
        y_data = item_df["relative_frequency"].to_numpy(dtype=float)

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

        if DISPLAY_LIN_REG and len(item_df) >= 2 and np.std(x_data) > 0:
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

    scatter_ax.set_title(PLOT_NAME_MAP["scatter"])
    scatter_ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
    scatter_ax.set_ylabel(Y_AXIS_NAME_MAP["relative_frequency"])
    scatter_ax.grid(True, alpha=0.25, linestyle="--")
    scatter_ax.legend()

    for ax, item_id in zip(axes[1:], item_ids_to_render):
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
                n_points=("relative_frequency", "size"),
                mean_relative_frequency=("relative_frequency", "mean"),
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
                item_binned["mean_relative_frequency"].to_numpy(dtype=float),
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
        ax.set_ylabel(Y_AXIS_NAME_MAP["relative_frequency_binned"])
        ax.grid(True, alpha=0.25, linestyle="--")

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.955, FIGURE_DESCRIPTION, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.04, 0.97, 0.92))

    plt.show()

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=250)
    fig.savefig(output_path_svg)
    plt.close(fig)


if __name__ == "__main__":
    main()
