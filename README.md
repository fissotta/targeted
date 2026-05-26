# QIIME 2 Amplicon Comparative Pipeline (`v2026-04-02a`)

A resumable QIIME 2 comparative pipeline for paired-end amplicon sequencing data. The workflow processes the same input dataset through two parallel feature-generation strategies: **DADA2** for ASV inference and **VSEARCH** for *de novo* OTU clustering.

The pipeline supports multi-database taxonomic classification, checkpoint-based recovery, downstream statistical analysis, differential abundance testing, and publication-ready SVG plots.

---

## Overview

```text
Paired-end FASTQ files
        |
        v
      Import
        |
        +-----------------------------+
        |                             |
        v                             v
      DADA2                        VSEARCH
  ASV inference                OTU clustering
        |                             |
        v                             v
 Low-frequency filtering      Low-frequency filtering
        |                             |
        v                             v
 Taxonomy classification      Taxonomy classification
 GG2024 / SILVA / GTDB        GG2024 / SILVA / GTDB
        |                             |
        +-------------+---------------+
                      |
                      v
                    Stats
      Alpha diversity, beta diversity,
      PERMANOVA, ANOSIM, PERMDISP,
      DESeq2, edgeR
                      |
                      v
                    Plots
       SVG figures for diversity,
       taxonomy, and differential abundance
```

---

## Key Features

### Dual Processing Strategy

The same paired-end FASTQ input is processed using two independent approaches:

- **DADA2**, producing amplicon sequence variants.
- **VSEARCH**, producing *de novo* clustered OTUs.

This enables direct comparison between denoising-based and clustering-based feature tables.

### Resumable Execution

Each major step writes checkpoint files. If the pipeline stops or fails, execution can resume from a selected stage without recomputing previous steps.

Available checkpoint stages:

- `import`
- `dada2`
- `vsearch`
- `stats`
- `plots`

### Multi-Classifier Taxonomy

The pipeline supports taxonomic classification against:

- Greengenes2 / GG2024
- SILVA
- GTDB

Each classifier is applied separately to both DADA2 and VSEARCH outputs.

### Downstream Statistical Analysis

The R-based downstream layer supports:

- Alpha-diversity comparisons
- Beta-diversity ordination
- Distance-based statistical analyses
- PERMANOVA
- ANOSIM
- PERMDISP
- DESeq2 differential abundance analysis
- edgeR differential abundance analysis

Batch covariates can be included when available.

### Publication-Ready Outputs

Figures are exported as vector graphics, mainly in SVG format, suitable for manuscripts, reports, presentations, and GitHub documentation.

---

## Requirements

### Operating System

The pipeline is intended for:

- Linux
- macOS
- HPC environments using schedulers such as Slurm

### Required Command-Line Tools

The following executables must be available in `PATH`:

- `qiime`
- `biom`
- `python3`
- `Rscript`
- `awk`
- `sed`
- `grep`
- `sort`
- `tee`

### Required Python Packages

Python 3 must include:

- `pandas`

Install with:

```bash
python3 -m pip install pandas
```

### Required R Packages

Install the required CRAN packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "tibble", "purrr", "stringr", "forcats", "vegan", "ape", "svglite", "scales", "colorspace"))
```

Install the required Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

---

## Classifier Files

By default, the pipeline expects classifier files in a local `CLASSIFIERS/` directory.

Expected files:

```text
CLASSIFIERS/
├── GG2024.09.backbone.v4.nb.sklearn-1.4.2_2024.5.qza
├── silva-138-99-nb-classifier_2024.5.qza
└── gtdb_classifier_r220_2024.5.qza
```

Classifier paths can be provided explicitly through command-line arguments when supported by the script.

---

## Input Files

### FASTQ Directory

The input FASTQ directory must contain paired-end reads in Casava 1.8 format.

Example:

```text
fastq/
├── Sample1_S1_L001_R1_001.fastq.gz
├── Sample1_S1_L001_R2_001.fastq.gz
├── Sample2_S2_L001_R1_001.fastq.gz
└── Sample2_S2_L001_R2_001.fastq.gz
```

### Metadata File

The metadata file must be tab-delimited.

Minimum required columns:

```text
sample-id
group column selected with --group-column
```

Example:

```text
sample-id	Treatment	SequencingRun
Sample1	Control	Run1
Sample2	Treatment	Run1
Sample3	Control	Run2
Sample4	Treatment	Run2
```

---

## Command-Line Usage

### Basic Full Run

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --threads 16
```

### Full Run with Batch Covariate

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --batch-column SequencingRun --output-dir qiime2_amplicon_out --threads 16
```

### Resume from VSEARCH

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --start-from vsearch --threads 16
```

### Resume from Statistics

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --start-from stats --threads 16
```

### Run Only Selected Stages

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --run-stages stats,plots --threads 16
```

### Force Rerun of Selected Stages

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --force-stages stats,plots --threads 16
```

### Stop After DADA2

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --stop-after dada2 --threads 16
```

---

## Parameters

| Flag | Value | Description | Default |
|---|---:|---|---|
| `--input-fastq` | `DIR` | Required. Directory containing paired-end Casava 1.8 FASTQ files. | None |
| `--metadata` | `FILE` | Required. Tab-delimited metadata file. | None |
| `--group-column` | `STRING` | Required. Metadata column used for grouping in statistics and plots. | None |
| `--sample-id-column` | `STRING` | Metadata column containing unique sample identifiers. | `sample-id` |
| `--batch-column` | `STRING` | Optional batch/covariate column for adjusted analyses. | None |
| `--output-dir` | `DIR` | Output directory. | `qiime2_amplicon_out` |
| `--threads` | `INT` | Number of CPU threads. | `16` |
| `--trunc-len-f` | `INT` | Forward-read truncation length for DADA2. | `270` |
| `--trunc-len-r` | `INT` | Reverse-read truncation length for DADA2. | `270` |
| `--vsearch-perc-id` | `FLOAT` | Identity threshold for VSEARCH *de novo* clustering. | `0.97` |
| `--start-from` | `STAGE` | Start execution from a specific stage. | `import` |
| `--stop-after` | `STAGE` | Stop execution after a specific stage. | `plots` |
| `--run-stages` | `CSV` | Run only selected stages. Overrides `--start-from` and `--stop-after`. | None |
| `--force-stages` | `CSV` | Force rerun of selected stages even if checkpoints exist. | None |

---

## Pipeline Stages

| Stage | Description |
|---|---|
| `import` | Imports paired-end FASTQ files into QIIME 2 artifacts. |
| `dada2` | Runs DADA2 denoising, feature-table generation, filtering, and taxonomy classification. |
| `vsearch` | Runs trimming, read merging, dereplication, *de novo* clustering, filtering, and taxonomy classification. |
| `stats` | Runs diversity statistics, multivariate tests, and differential abundance analyses. |
| `plots` | Generates SVG plots for taxonomy, alpha diversity, beta diversity, and statistical summaries. |

---

## Output Structure

```text
qiime2_amplicon_out/
├── RUN_SUMMARY.txt
├── metadata.tsv
├── demux.qza
├── checkpoints/
│   ├── import.done
│   ├── dada2.done
│   ├── vsearch.done
│   ├── stats.done
│   ├── plots.done
│   └── current_stage.txt
├── logs/
│   ├── pipeline.log
│   └── stage_summary.tsv
├── dada2/
│   ├── table.qza
│   ├── rep_seqs.qza
│   ├── taxonomy/
│   ├── alpha/
│   ├── beta/
│   ├── exports/
│   │   ├── raw_table/
│   │   │   └── feature-table.tsv
│   │   └── consensus/
│   ├── stats/
│   │   ├── read_tracking.tsv
│   │   └── differential/
│   └── plots/
│       ├── alpha/
│       ├── beta/
│       └── taxonomy/
├── vsearch/
│   ├── table.qza
│   ├── rep_seqs.qza
│   ├── taxonomy/
│   ├── alpha/
│   ├── beta/
│   ├── exports/
│   │   ├── raw_table/
│   │   │   └── feature-table.tsv
│   │   └── consensus/
│   ├── stats/
│   │   ├── read_tracking.tsv
│   │   └── differential/
│   └── plots/
│       ├── alpha/
│       ├── beta/
│       └── taxonomy/
└── scripts/
```

---

## Main Outputs

### QIIME 2 Artifacts

```text
dada2/table.qza
dada2/rep_seqs.qza
vsearch/table.qza
vsearch/rep_seqs.qza
```

### Exported Feature Tables

```text
dada2/exports/raw_table/feature-table.tsv
vsearch/exports/raw_table/feature-table.tsv
```

### Read Tracking Tables

```text
dada2/stats/read_tracking.tsv
vsearch/stats/read_tracking.tsv
```

### Differential Abundance Results

```text
dada2/stats/differential/
vsearch/stats/differential/
```

### Plots

```text
dada2/plots/
vsearch/plots/
```

Expected plot categories:

- Alpha diversity
- Beta diversity
- Taxonomic composition
- Differential abundance
- Read tracking
- Classifier comparisons

---

## Recommended Execution with `nohup`

For long runs, use `nohup`:

```bash
nohup bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --threads 16 > qiime2_amplicon_out.nohup.log 2>&1 &
```

Monitor progress with:

```bash
tail -f qiime2_amplicon_out.nohup.log
```

---

## Checkpoint and Resume Logic

The pipeline writes `.done` files after successful completion of each stage.

Example:

```text
checkpoints/import.done
checkpoints/dada2.done
checkpoints/vsearch.done
checkpoints/stats.done
checkpoints/plots.done
```

If a stage has already completed, it is skipped unless forced with `--force-stages`.

Example:

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --output-dir qiime2_amplicon_out --force-stages stats,plots
```

---

## Typical Workflow

A recommended complete workflow is:

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --batch-column SequencingRun --output-dir qiime2_amplicon_out --threads 16
```

If the run stops during VSEARCH:

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --batch-column SequencingRun --output-dir qiime2_amplicon_out --start-from vsearch --threads 16
```

If only plots need to be regenerated:

```bash
bash qiime2_amplicon_comparative_pipeline_resumable_fixed_v8.sh --input-fastq fastq --metadata metadata.tsv --group-column Treatment --batch-column SequencingRun --output-dir qiime2_amplicon_out --run-stages plots --force-stages plots --threads 16
```

---

## Notes

- Input FASTQ files must match the expected paired-end Casava naming format.
- Metadata sample identifiers must match FASTQ sample names after QIIME 2 import parsing.
- The `--group-column` must exist in the metadata file.
- The `--batch-column` is optional but recommended when sequencing runs, extraction batches, or technical covariates are known.
- DADA2 truncation lengths should be adjusted according to read-quality profiles.
- VSEARCH clustering identity defaults to `0.97`, corresponding to 97% OTU clustering.
- Differential abundance methods require sufficient replication per group.
- PERMANOVA, ANOSIM, and PERMDISP require valid beta-diversity distance matrices and compatible metadata grouping.

---

## Citation

If this pipeline is used in a publication, cite this repository and the tools used by the workflow, including QIIME 2, DADA2, VSEARCH, the selected taxonomy databases, DESeq2, edgeR, and vegan.

---

## License
- MIT License
---
