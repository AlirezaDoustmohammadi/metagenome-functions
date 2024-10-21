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


## PICRUSt2 Pipeline
a Nextflow-based pipeline that wraps PICRUSt2 tools for metagenome prediction. The pipeline automates the process of placing reads on a reference tree, performing hidden-state prediction of gene families, and generating metagenome and pathway-level predictions based on 16S rRNA gene sequence data.

### Requirements
PICRUSt2: Used for metagenome prediction and pathway inference.

### Example Usage
To run the pipeline, use the following command:
```
nextflow run picrust2_pipeline.nf \
    --biom_table /path/to/feature-table.biom \
    --seq /path/to/representative-sequences.fasta \
    --proj_dir /path/to/output/ \
    --cpu 8
```
This command will:
- Place sequences into a reference tree.
- Predict hidden states for gene families (KOs, ECs).
- Generate metagenome predictions for KOs and ECs.
- Perform pathway-level inference based on predicted metagenome abundances.

## KO to KEGG Pathway Mapper
<a href='KEGG%20db/ko_to_kegg_pathway_mapper.py'>ko_to_kegg_pathway_mapper.py</a> maps KEGG Orthology (KO) functions to their corresponding KEGG pathways
### Example Usage:
Ensure your KO prediction file is in a `.tsv.gz` format.</br>
Edit the <a href='KEGG%20db/ko_to_kegg_pathway_mapper.py#L7'>script to point</a> to the correct path for your input file:
```
df = pd.read_csv('path_to_your_file/KO_pred_metagenome_unstrat.tsv.gz', sep='\t', compression='gzip')
```
