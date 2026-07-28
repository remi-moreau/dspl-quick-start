
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
    2: "Motif PAIRCODE"
}


# ---- DATA STRUCTURE ----

EXPECTED_COLUMN_NAMES_IN_THE_CSV = [
    None
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
    0: None,
    1: None,
    2: None,
}

Y_AXIS_STYLE_MAP = {
    "coverage": None,
    "hamming_distance_normalized_at_decoding": None,
    "perfectly_decoded_payload_ratio_at_decoding": None,
    "estimated_sequencing_duration_at_decoding": None
}


##############
#### CODE ####
##############





