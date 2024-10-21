# prediction of metagenome functions

## Qiime2 Preprocessing Pipeline

a Nextflow pipeline for preprocessing biom and sequence files using Qiime2. The goal is to convert, filter, and summarize microbiome feature tables and sequences, particularly for applications like microbiome analysis or machine learning feature extraction.</br>
This pipeline automates the preprocessing of biom and sequence files, with a focus on filtering and converting files for downstream analysis. It allows you to:
- Convert BIOM and FASTA files into Qiime2-compatible formats.
- Filter features and samples based on user-defined abundance thresholds.
- Export filtered data back to BIOM and FASTA formats.
- Summarize the final output to help analyze the distribution of features.

### Requirements
- Qiime2: A microbiome bioinformatics platform (specifically the 2024.5 version used here).
- Biom-format: Command-line utilities for BIOM table handling.

### Example Usage
Run the pipeline with the following command:
```
nextflow run preprocessing.nf \
    --biom_table /path/to/feature-table.biom \
    --seq /path/to/representative-sequences.fasta \
    --input_dir /path/to/input/ \
    --proj_dir /path/to/project/ \
    --min_feature_abundance 10 \
    --min_samples_abundance 1000
```
This command will:
- Convert the BIOM and FASTA files to Qiime2 format.
- Filter features with a minimum abundance of 10 and samples with a minimum abundance of 1000.
- Output the results into the specified project directory.

