"""
Settings for the experience six_images (encoding and decoding of six images of the JPEG AIC 3 dataset).
This file contains all the settings required for the different steps of the experience, such as encoding, decoding and plotting.
The settings are packed in dictionaries, such as PACKED_ENCODING_SETTINGS, that can be easily imported and used in the different scripts.
"""



#################################################################################################
#################################################################################################
##################################       GLOBAL SETTINGS      ###################################
#################################################################################################
#################################################################################################

# Experience name and description
EXPERIENCE_NAME = "synthesis_bench"
EXPERIENCE_DESCRIPTION = "Encoding and decoding of three copies of the same image" \
                        " (the bird, i.e. image 5, of the JPEG AIC-3 dataset). \n· The" \
                        " first copy is the one of reference, it uses JPEGDNA classical " \
                        "encoding. \n· The second copy is encoded using JPEGDNA but with free " \
                        "energy constraints, to see their effects on decoding. \n· The third copy" \
                        " is encoded using Motif PAIRCODE."

# Python executable of the virtual environment specific to the use of the jpegdna VM
from pathlib import Path
VM_PYTHON_EXECUTABLE = str(
    Path.home() / "miniconda3" / "envs" / "jpegdna_vm_py311" / "bin" / "python"
)

# Input/output data directories
INPUT_DATA_PATH = f"input_data/{EXPERIENCE_NAME}"
OUTPUT_DATA_PATH = f"output_data/{EXPERIENCE_NAME}"






#################################################################################################
#################################################################################################
###############################      ENCODING SETTINGS   ########################################
#################################################################################################
#################################################################################################

# There are two types of encoding runs:
# - "outer_encoding" run which only encodes the payload map, without performing the inner 
#   encoding and formatted oligo map generation. 
# - "full_encoding" run which performs the full encoding process, from the payload map encoding 
#   to the formatted oligo map generation. It allows to output the formatted oligo map which is 
#   the one that will be synthesized in the wetlab and sequenced.

###############################
# COMMON ENCODING RUNS SETTINGS
###############################

# The following settings are common for both "outer_encoding" and "full_encoding" runs. 

# Default outer and inner encoding settings
# -----------------------------------------
DEFAULT_OUTER_ENCODING_SETTINGS = {
	"name": "JPEGDNA",
	"jxl_quality": 78,
	"jxl_skip": False,
	"jxl_binary_path": None,
	"raptor_symbol_size": 29,
	"raptor_max_source_symbols_per_block": 2600,
	"raptor_num_source_blocks": None,
	"raptor_overhead": 0.1,
	"nucleotide_decoder": "A",
	"ext_header_rs_bytes": 17,
	"ext_main_section_rs_bytes": 6,
	"jpegdna_engine_root": None,
	"libjxl_tools_path": None,
	"libjxl_ld_library_path": None,
    "jpegdna_python_executable": VM_PYTHON_EXECUTABLE,
    # "jpegdna_python_executable": None,
}
DEFAULT_INNER_ENCODING_SETTINGS = {
    "name": "REGION_FORMATTING",
}


# Per-item specs and outer/inner encoding settings
# ------------------------------------------------
LEFT_PRIMER = "AGAGCGGCGTATTGTATTCG" 
RIGHT_PRIMER = "AGTTCACGTCCAGTCAGAGT"
HEADER_PAYLOAD_LENGTH_JPEGDNA = 155 
BODY_PAYLOAD_LENGTH_JPEGDNA = 156 # NOTE: Will be used by the clustering-consensus instance
                                  # to return a consensus of this exact length
DATA_PAYLOAD_LENGTH_MOTIF_PAIRCODE = 162
MOTIF_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE = 157
BYTE_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE = 157

ITEMS = {
    "type": "image",
    "folder": f"{INPUT_DATA_PATH}/images/",
    "items": [ # a LIST of the items we are encoding/decoding. The ORDER of the list must be KEPT during 
            # the ENTIRE encoding-decoding process, thus index i in the list of an item must be seen as 
            # its UNIQUE id.
        {
            "name":"Bird_ref_JPEG_DNA",
            "input_folder":None, # Populated by the orchestrator with the above "folder" key when is None. 
            "input_path":"00005_560x888.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 78,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "ATCCACA", "payload_len": HEADER_PAYLOAD_LENGTH_JPEGDNA},
                    {"name": "body", "region_len": None, "index": "TACCAC", "payload_len": BODY_PAYLOAD_LENGTH_JPEGDNA},
                ]
            }
        },
        {
            "name":"Bird_delta_G_JPEG_DNA",
            "input_folder":None,
            "input_path":"00005_560x888.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 78,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "AATGATA", "payload_len": HEADER_PAYLOAD_LENGTH_JPEGDNA},
                    {"name": "body", "region_len": None, "index": "TGTGAT", "payload_len": BODY_PAYLOAD_LENGTH_JPEGDNA},
                ]
            }
        },
        {
            "name":"Bird_MOTIF_PAIRCODE",
            "input_folder":None,
            "input_path":"00005_560x888.png",
            "outer_encoding_settings":{
                "name": "MOTIF_PAIRCODE",
                "per_region": [
                    {"name": "data", "oligo_len": DATA_PAYLOAD_LENGTH_MOTIF_PAIRCODE},
                    {"name": "motif_mapping", "oligo_len": MOTIF_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE},
                    {"name": "byte_mapping", "oligo_len": BYTE_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE},
                ],
                "motif_engine_root": None,
                "motif_libjxl_tools_path": None,
                "motif_libjxl_ld_library_path": None,
                "motif_python_executable": VM_PYTHON_EXECUTABLE,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "data", "region_len": None, "index": None, "payload_len": DATA_PAYLOAD_LENGTH_MOTIF_PAIRCODE},
                    {"name": "motif_mapping", "region_len": None, "index": "ATAAG", "payload_len": MOTIF_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE},
                    {"name": "byte_mapping", "region_len": None, "index": "CCATT", "payload_len": BYTE_MAPPING_PAYLOAD_LENGTH_MOTIF_PAIRCODE}
                ]
            }
        },
    ]
}



#############################
# OUTER ENCODING RUN SETTINGS
#############################

# Packed settings (what the encoding orchestrator will see)
# ---------------
PACKED_OUTER_ENCODING_SETTINGS = {
    "experience_name": EXPERIENCE_NAME,
    "description": f"Outer encoding run of {EXPERIENCE_NAME} experience.\nEncode {str(len(ITEMS['items']))} items into the same number of encoded payload map. Item names: {[item['name'] for item in ITEMS['items']]}. Protocol used: JPEGDNA for the first two items, MOTIF_PAIRCODE for the third item.",
    "output_data_path": f"{OUTPUT_DATA_PATH}/encoding/",
    "items": ITEMS,
    "default_outer_encoding_settings": DEFAULT_OUTER_ENCODING_SETTINGS,
    "save_payload_map": True,
    "save_run_metadata": True,
}


############################
# FULL ENCODING RUN SETTINGS
############################

# Packed settings (what the encoding orchestrator will see)
# ---------------
PACKED_ENCODING_SETTINGS = {
    "experience_name": EXPERIENCE_NAME,
    "description": f"Full encoding run of {EXPERIENCE_NAME} experience.\nEncode {str(len(ITEMS['items']))} items into the same number of encoded payload map AND formatted oligo map. Item names: {[item['name'] for item in ITEMS['items']]}. Outer protocol used: JPEGDNA for the first two items, MOTIF_PAIRCODE for the third item.",
    "output_data_path": f"{OUTPUT_DATA_PATH}/encoding/",
    "items": ITEMS,
    "default_outer_encoding_settings": DEFAULT_OUTER_ENCODING_SETTINGS,
    "default_inner_encoding_settings": DEFAULT_INNER_ENCODING_SETTINGS,
    "save_formatted_oligo_map": True,
    "save_payload_map": True,
    "save_run_metadata": True,
}

