# WETLAB DATA ANALYSIS WITH `dspl`

There are four barcodes in total, three for experience 1 (`synthesis_bench`) and one for experience 2 (`six_images`):
1. `synthesis_bench`:
     * `barcode01_agilent`: Agilent sureprint: Bird JPEGDNA DELTAG + Motif Paircode
     * `barcode02_dynegene`: Dynegene: Bird JPEGDNA DELTAG + Motif Paircode
     * `barcode03_genscript`: Genscript: Bird JPEGDNA DELTAG + Motif Paircode
2. `six_images`:
     * `barcode04_six_images`: Agilent Hifi - JPEGDNA Wetlab post encapsulation

We put their data in four databases, one per barcode.

## Setup

You must first install the `dspl` package (clone the repo https://github.com/remi-moreau/dna-storage-pipeliner and follow its own README.md instructions).


## Usage

As described in the `dspl` own README.md, clone this repository (https://github.com/remi-moreau/dspl-quick-start):
```bash
git clone https://github.com/remi-moreau/dspl-quick-start.git
cd dspl-quick-start
```

First, fill-in the `decoding_*.yml` config files with the pathes to the Nanopore reads and other inputs you want to use.

Then, review your decoding configuration `decoding_config_synthesis_bench.yml` file. 
* In particular, verify that the parameter `decoding_settings -> input_settings -> read_folder_path` is set to the read folder containing the fastq you want to decode. By default, it is set to `progressive_decoding_wetlab_reads/barcode01_agilent`.
* You can also adjust multi_pass settings and inner/outer decoding settings.


Now that you are in the `dspl-quick-start` folder and you have setup the config files, you can run the commands of this document. They assume you are in the `dspl-quick-start` folder, and that you have installed the `dspl` package.

BE CAREFUL TO ADAPT DATABASE AND OTHER FILE PATHS TO YOUR OWN PATHS BOTH IN CLI AND .YML CONFIG FILES BEFORE RUNNING COMMANDS!


## Table of contents
- [WETLAB DATA ANALYSIS WITH `dspl`](#wetlab-data-analysis-with-dspl)
  - [Setup](#setup)
  - [Usage](#usage)
  - [Table of contents](#table-of-contents)
  - [1. Decoding `synthesis_bench`](#1-decoding-synthesis_bench)
    - [1.1 `barcode01_agilent`](#11-barcode01_agilent)
      - [Encoding run](#encoding-run)
      - [Decoding with JPEGDNA ref and delta\_g only](#decoding-with-jpegdna-ref-and-delta_g-only)
      - [Decoding with motif paircode (one fastq per pass)](#decoding-with-motif-paircode-one-fastq-per-pass)
      - [Decoding with motif paircode (multiple fastqs per pass)](#decoding-with-motif-paircode-multiple-fastqs-per-pass)
    - [1.2 `barcode02_dynegene`](#12-barcode02_dynegene)
      - [Encoding run](#encoding-run-1)
      - [Decoding with JPEGDNA ref and delta\_g only](#decoding-with-jpegdna-ref-and-delta_g-only-1)
    - [1.3 `barcode03_genscript`](#13-barcode03_genscript)
      - [Encoding run](#encoding-run-2)
      - [Decoding with JPEGDNA ref and delta\_g only](#decoding-with-jpegdna-ref-and-delta_g-only-2)
  - [2. Decoding `six_images`](#2-decoding-six_images)
      - [Encoding run](#encoding-run-3)
      - [Decoding with all images EXCEPT burger](#decoding-with-all-images-except-burger)
      - [Decoding with all images including burger](#decoding-with-all-images-including-burger)
  - [3. Plotting](#3-plotting)
    - [3.1 Plot `delta_g` metrics](#31-plot-delta_g-metrics)
    - [3.2 Plot JPEGDNA before/after encapsulation metrics](#32-plot-jpegdna-beforeafter-encapsulation-metrics)
    - [3.3 Plot Motif PAIRCODE metrics](#33-plot-motif-paircode-metrics)
    - [3.4 Plot Progressive decoding metrics](#34-plot-progressive-decoding-metrics)




## 1. Decoding `synthesis_bench`

### 1.1 `barcode01_agilent`

#### Encoding run

Create database and add encoding run
```bash
dspl db add-full-encoding database/barcode01_agilent.db \
--run output_data/synthesis_bench/encoding/full_encoding_15/ \
--experience-id synthesis_bench \
--delta-g-map input_data/synthesis_bench/delta_g_files/delta-g_mapping.csv
```

#### Decoding with JPEGDNA ref and delta_g only

Decode JPEGDNA ref and delta_g with alignment from read fastq files (NOTE: before, you must set the path to the reads in the `decoding_barcode01_agilent.yml` config file).
```bash
dspl run decoding decoding_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_decoding \
&& dspl run metrics-computation metrics_computation_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_decoding
```

Variant with CDHIT instead of alignment, using the reads just ingested in the database:
```bash
dspl run decoding decoding_barcode01_agilent.yml \
--labels barcode01_agilent_cdhit \
--override config/synthesis_bench/barcode01_agilent/cdhit.yml config/synthesis_bench/barcode01_agilent/from_read_pool.yml \
&& dspl run metrics-computation metrics_computation_barcode01_agilent.yml \
--labels barcode01_agilent_cdhit
```

#### Decoding with motif paircode (one fastq per pass)

Decode all three images (including motif paircode) from the previously ingested read pool and compute metrics:
```bash
dspl run decoding decoding_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_with_motif_paircode \
--override config/synthesis_bench/barcode01_agilent/from_read_pool.yml config/synthesis_bench/decode_also_motif_paircode.yml config/no_loop.yml \
&& dspl run metrics-computation metrics_computation_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_with_motif_paircode
```

#### Decoding with motif paircode (multiple fastqs per pass)
Decode all three images (including motif paircode) from the previously ingested read pool and compute metrics, but this time with multiple fastqs per pass (see `config/multi_pass_n_fastq_per_pass.yml`):
```bash
dspl run decoding decoding_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_with_motif_paircode_n_fastq_per_pass \
--override config/synthesis_bench/barcode01_agilent/from_read_pool.yml config/synthesis_bench/decode_also_motif_paircode.yml config/no_loop.yml config/n_fastq_per_pass.yml \
&& dspl run metrics-computation metrics_computation_barcode01_agilent.yml \
--labels barcode01_agilent_alignment_with_motif_paircode_n_fastq_per_pass
```
### 1.2 `barcode02_dynegene`

#### Encoding run

Create database and add encoding run
```bash
dspl db add-full-encoding database/barcode02_dynegene.db \
--run output_data/synthesis_bench/encoding/full_encoding_15/ \
--experience-id synthesis_bench \
--delta-g-map input_data/synthesis_bench/delta_g_files/delta-g_mapping.csv
```


#### Decoding with JPEGDNA ref and delta_g only

```bash
dspl run decoding decoding_barcode02_dynegene.yml \
--labels barcode02_dynegene_alignment_decoding \
&& dspl run metrics-computation metrics_computation_barcode02_dynegene.yml \
--labels barcode02_dynegene_alignment_decoding
```

### 1.3 `barcode03_genscript`

#### Encoding run

Create database and add encoding run
```bash
dspl db add-full-encoding database/barcode03_genscript.db \
--run output_data/synthesis_bench/encoding/full_encoding_15/ \
--experience-id synthesis_bench \
--delta-g-map input_data/synthesis_bench/delta_g_files/delta-g_mapping.csv
```


#### Decoding with JPEGDNA ref and delta_g only

```bash
dspl run decoding decoding_barcode03_genscript.yml \
--labels barcode03_genscript_alignment_decoding \
&& dspl run metrics-computation metrics_computation_barcode03_genscript.yml \
--labels barcode03_genscript_alignment_decoding
```

## 2. Decoding `six_images`

#### Encoding run

Create database and add encoding run
```bash
dspl db add-full-encoding database/barcode04_six_images.db \
--run output_data/six_images/encoding/full_encoding_11/ \
--experience-id six_images
```

#### Decoding with all images EXCEPT burger

```bash
dspl run decoding decoding_barcode04_six_images.yml \
--labels barcode04_six_images_alignment_decoding \
&& dspl run metrics-computation metrics_computation_barcode04_six_images.yml \
--labels barcode04_six_images_alignment_decoding
```

#### Decoding with all images including burger

```bash
dspl run decoding decoding_barcode04_six_images.yml \
--labels barcode04_six_images_alignment_decoding_with_burger \
--override config/six_images/with_burger.yml \
&& dspl run metrics-computation metrics_computation_barcode04_six_images.yml \
--labels barcode04_six_images_alignment_decoding_with_burger
```

## 3. Plotting

The plottings scripts are in folder `plot_dspl_database`. The main script is `main.py`, which takes as input a YAML configuration file. The configuration file specifies which plotting scripts to run and their settings.

### 3.1 Plot `delta_g` metrics

```bash
python plot_dspl_database/main.py plot_config_delta_g.yml
```

### 3.2 Plot JPEGDNA before/after encapsulation metrics

```bash
python plot_dspl_database/main.py plot_config_jpegdna_before_after_encapsulation.yml
```


### 3.3 Plot Motif PAIRCODE metrics

```bash
python plot_dspl_database/main.py plot_config_motif.yml
```

### 3.4 Plot Progressive decoding metrics

```bash
python plot_dspl_database/main.py plot_config_progressive.yml
```