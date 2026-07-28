################
#### CONFIG ####
################

# ---- DATA PATH ----
from pathlib import Path

DATA_PATH_PREFIX = "../../database/exports/"

DATA_PATH_SUFFIX = "barcode01_agilent/metrics_view_cluster/at_decoding.csv"

DATA_PATH = Path(DATA_PATH_PREFIX) / Path(DATA_PATH_SUFFIX)

ITEM_ID_NAME_MAP = {
    0: "JPEG DNA reference",
    1: "JPEG DNA delta G",
}


# ---- DATA SPECS ----

EXPECTED_COLUMN_NAMES_IN_THE_CSV = [
    "label",
    "item_id",
    "region_id",
    "position_id",
    "delta_g",
    "count_used_for_consensus",
]

LABELS_NAME_MAP = {
    "barcode01_agilent_alignment_decoding": "Read pool from supplier 1 (Agilent)"
}

# ---- OUTPUT PATH ----

OUTPUT_PATH_PNG = "../plots/delta_g_correlation.png"
OUTPUT_PATH_SVG = "../plots/delta_g_correlation.svg"

# ---- NAMES ----

FIGURE_TITLE = "Cluster size at decoding against delta G."

FIGURE_DESCRIPTION = f"Cluster size at decoding against delta G, for label {'barcode01_agilent_alignment_decoding'}."

PLOT_NAME_MAP = {
    "delta_g_correlation": "Coverage at decoding against delta G (mean over runs by reference)",
}

X_AXIS_NAME_MAP = {
    "delta_g": "Delta G",
}

Y_AXIS_NAME_MAP = {
    "count_used_for_consensus": "Coverage at decoding (mean count_used_for_consensus)",
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
    }
}


##############
#### CODE ####
##############

import matplotlib.pyplot as plt
import pandas as pd


def _resolve_path_from_script(relative_path: Path) -> Path:
    script_dir = Path(__file__).resolve().parent
    return (script_dir / relative_path).resolve()


def _validate_columns(df: pd.DataFrame, expected_columns: list[str]) -> None:
    missing_columns = [column for column in expected_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing expected columns in CSV: {missing_columns}")


def main() -> None:
    data_path = _resolve_path_from_script(DATA_PATH)
    output_path_png = _resolve_path_from_script(Path(OUTPUT_PATH_PNG))
    output_path_svg = _resolve_path_from_script(Path(OUTPUT_PATH_SVG))

    df = pd.read_csv(data_path)
    _validate_columns(df, EXPECTED_COLUMN_NAMES_IN_THE_CSV)

    df = df[df["label"].isin(LABELS_NAME_MAP.keys())]
    df = df[df["item_id"].isin(ITEM_ID_NAME_MAP.keys())].copy()

    if df.empty:
        raise ValueError("No data available after filtering by label and item_id.")

    df["delta_g"] = pd.to_numeric(df["delta_g"], errors="coerce")
    df["count_used_for_consensus"] = (
        pd.to_numeric(df["count_used_for_consensus"], errors="coerce").fillna(0)
    )
    df = df.dropna(subset=["delta_g"])

    # One point per reference and item: mean coverage across runs.
    reference_level_df = (
        df.groupby(["item_id", "region_id", "position_id"], as_index=False)
        .agg(
            delta_g=("delta_g", "first"),
            count_used_for_consensus=("count_used_for_consensus", "mean"),
        )
    )

    fig, ax = plt.subplots(figsize=(11.5, 7.5))
    style = PLOT_STYLE_MAP["delta_g_correlation"]

    for item_id, item_name in ITEM_ID_NAME_MAP.items():
        item_df = reference_level_df[reference_level_df["item_id"] == item_id]
        if item_df.empty:
            continue

        ax.scatter(
            item_df["delta_g"],
            item_df["count_used_for_consensus"],
            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
            label=f"{item_name} (n={len(item_df)})",
            marker=style["marker"],
            s=style["s"],
            alpha=style["alpha"],
            edgecolors=style["edgecolors"],
        )

    ax.set_title(PLOT_NAME_MAP["delta_g_correlation"])
    ax.set_xlabel(X_AXIS_NAME_MAP["delta_g"])
    ax.set_ylabel(Y_AXIS_NAME_MAP["count_used_for_consensus"])
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.legend()

    fig.suptitle(FIGURE_TITLE, fontsize=16, y=0.985)
    fig.text(0.5, 0.948, FIGURE_DESCRIPTION, ha="center", va="top", wrap=True, fontsize=11)
    fig.tight_layout(rect=(0.03, 0.05, 0.97, 0.9))

    output_path_png.parent.mkdir(parents=True, exist_ok=True)
    output_path_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path_png, dpi=250)
    fig.savefig(output_path_svg)
    plt.close(fig)


if __name__ == "__main__":
    main()
