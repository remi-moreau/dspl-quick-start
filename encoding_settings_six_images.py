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
EXPERIENCE_NAME = "six_images"
EXPERIENCE_DESCRIPTION = "Encoding and decoding of six images of the JPEG AIC 3 dataset, " \
"using JPEGDNA as the outer coding protocol. " \

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
}
DEFAULT_INNER_ENCODING_SETTINGS = {
    "name": "REGION_FORMATTING",
}


# Per-item specs and outer/inner encoding settings
# ------------------------------------------------
LEFT_PRIMER = "AGAGCGGCGTATTGTATTCG"
RIGHT_PRIMER = "AGTTCACGTCCAGTCAGAGT"
HEADER_PAYLOAD_LENGTH = 155 
BODY_PAYLOAD_LENGTH = 156 # NOTE: Will be used by the clustering-consensus instance
                          # to return a consensus of this exact length

ITEMS = {
    "type": "image",
    "folder": f"{INPUT_DATA_PATH}/images/jxl_decoded_by_me/",
    "items": [ # a LIST of the items we are encoding/decoding. The ORDER of the list must be KEPT during 
            # the ENTIRE encoding-decoding process, thus index i in the list of an item must be seen as 
            # its UNIQUE id.
        {
            "name":"Chest",
            "input_folder":None, # Populated by the orchestrator with the above "folder" key when is None. 
            "input_path":"00001_1192x832.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 78,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "AATGATA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "TGTGAT", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            }
        },
        {
            "name":"Woman",
            "input_folder":None,
            "input_path":"00002_853x945.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 79,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "ATCAAGA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "TCCAAG", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            }
        },
        {
            "name":"Burger",
            "input_folder":None,
            "input_path":"00003_945x840.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 76,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "TTAAGTA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "GGAAGT", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            }
        },
        {
            "name":"Bird",
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
                    {"name": "header", "region_len": 1, "index": "ATCCACA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "TACCAC", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            }
        },
        {
            "name":"Night",
            "input_folder":None,
            "input_path":"00007_1600x1200.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 76,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "AGTGGTA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "GATGGT", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            }
        },
        {
            "name":"Day",
            "input_folder":None,
            "input_path":"00010_2592x1946.png",
            "outer_encoding_settings":{
                "name": "JPEGDNA",
                "jxl_quality": 79,
            },
            "inner_encoding_settings": {
                "name": "REGION_FORMATTING",
                "left_primer": LEFT_PRIMER,
                "right_primer": RIGHT_PRIMER,
                "regions": [
                    {"name": "header", "region_len": 1, "index": "AGTTTGA", "payload_len": HEADER_PAYLOAD_LENGTH},
                    {"name": "body", "region_len": None, "index": "TATTTG", "payload_len": BODY_PAYLOAD_LENGTH},
                ]
            },
        }
    ]
}


#############################
# OUTER ENCODING RUN SETTINGS
#############################

# Packed settings (what the encoding orchestrator will see)
# ---------------
PACKED_OUTER_ENCODING_SETTINGS = {
    "experience_name": EXPERIENCE_NAME,
    "description": f"Outer encoding run of {EXPERIENCE_NAME} experience.\nEncode {str(len(ITEMS['items']))} items into the same number of encoded payload maps. Item names: {[item['name'] for item in ITEMS['items']]}. Outer protocol used: JPEGDNA for all items.",
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
    "description": f"Full encoding run of {EXPERIENCE_NAME} experience.\nEncode {str(len(ITEMS['items']))} items into the same number of encoded payload maps and formatted oligo maps. Item names: {[item['name'] for item in ITEMS['items']]}. Outer protocol used: JPEGDNA for all items. Inner protocol used: REGION_FORMATTING for all items.",
    "output_data_path": f"{OUTPUT_DATA_PATH}/encoding/",
    "items": ITEMS,
    "default_outer_encoding_settings": DEFAULT_OUTER_ENCODING_SETTINGS,
    "default_inner_encoding_settings": DEFAULT_INNER_ENCODING_SETTINGS,
    "save_formatted_oligo_map": True,
    "save_payload_map": True,
    "save_run_metadata": True,
}

