# 1. USAGE

The YAML configuration path is required:

```bash
python plot_dspl_database/main.py plot_config.yml
```

# 2. FOR DEVELOPERS

## 2.1. `main.py`

Role: 
 - import necessary protocol(s) and pydantic data structure for .yml config validation
 - import all necessary script function (assuming they follow the right protocol in protocols.py) and assign them a canonical name which will be used to fail fast if the .yml config mentions an unknown script
 - parse CLI input and extract config from the given .yml 
 - fail fast if config is invalid (using pydantic validation)
 - run all required plotting scripts

**Import rules:** Import protocols, plotting scripts. DO NOT import `optional_script_utils.py`

## 2.2. `protocols.py`

Role: defining
 - protocol for scripts main function to follow
 - pydantic data structure for main .yml config validation

**Import rules:** Import NOBODY. This script is the business core.

## 2.3. Plotting scripts

A valid plotting script should provide a main class which follows the corresponding protocol defined in `protocols.py`.
In particular, the signature of its initialisation includes a dict which contains script-specific settings. This dict is not verified by the main.py, thus the script itself should verify it.


**Import rules:** optionally import `optional_script_utils.py`. Do NOT import `main.py`. Can import `protocols.py` to ensure settings structure.

## 2.4. `optional_script_utils.py`

Optional utility function and/or classes a script can optionally use. 

Plotting scripts export one multipage PDF containing one or more metadata pages followed by the configured figures, plus one PNG per figure. Figure and script settings are validated by strict Pydantic models owned by each plotting script.

Top-level `show_figures: true` displays each configured figure in YAML order before export. Close the figure window to continue; any interactive resize is preserved in the PNG and PDF. Metadata pages are not displayed.

Every figure accepts `same_x_scale_across_inputs` and `same_y_scale_across_inputs`. When enabled, the corresponding limits are harmonized across populated input/label axes after rendering.

Top-level `figure_width_per_input` and `figure_height` define the initial canvas size for every plotting script. The `read_pool_stats` distribution setting `read_count_values_per_bin: N` groups `N` consecutive integer read-count values in each bar class before converting the X axis to normalized coverage. Set it to `null` to retain Freedman-Diaconis binning on normalized coverage.

For `ref-coverage-at-ref-decoding_vs_delta-g_scatter`, metadata report Spearman rho and its p-value per input and item. Only finite decoded-reference pairs are included; graphical drop-outs at infinity are excluded.

## 2.5. Other rules

No subpackage. Every module must be at the package root.

Systematically use pydantic to define and validate any settings.

When building a new script, do not forget to use or complete `optional_script_utils.py`.

