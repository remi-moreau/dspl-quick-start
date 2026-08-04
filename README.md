# WETLAB DATA ANALYSIS WITH `dspl`

There are four barcodes in total, three for experience 1 (`synthesis_bench`) and one for experience 2 (`six_images`):
1. `synthesis_bench`:
     * `barcode01_agilent`: Agilent sureprint: Bird JPEGDNA DELTAG + Motif Paircode
     * `barcode02_dynegene`: Dynegene: Bird JPEGDNA DELTAG + Motif Paircode
     * `barcode03_genscript`: Genscript: Bird JPEGDNA DELTAG + Motif Paircode
2. `six_images`:
     * `barcode04_six_images`: Agilent Hifi - JPEGDNA Wetlab post encapsulation

We put their data in four databases, one per barcode.

## 1. `synthesis_bench`

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

## 2. `six_images`

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