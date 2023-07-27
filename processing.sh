#!/bin/bash

# start with `conda activate picrust2`

# Declare variables for biom and sequence files
raw_biom='feature-table.biom'
raw_seq='representative-sequences.fasta'
input_dir='data_pre_processing'
num_threads=12

# Extract filenames without extensions
raw_biom_filename="${raw_biom%.biom}"
raw_seq_filename="${raw_seq%.fna}"
raw_seq_filename="${raw_seq_filename%.fasta}"

# Function for placing reads into reference tree
place_seqs() {
  place_seqs.py -s ../${input_dir}/${raw_seq_filename}.filtered.features.samples.fasta -o out.tre -p ${num_threads} \
                --intermediate intermediate/place_seqs &> place_seqs.log
}

# Function for hidden-state prediction of gene families
predict_gene_families() {
  hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p ${num_threads} -n
  hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p ${num_threads}
}

# Function for generating metagenome predictions
generate_metagenome_predictions() {
  metagenome_pipeline.py -i ../${input_dir}/${raw_biom_filename}.filtered.features.samples.biom\
					     -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
                         -o EC_metagenome_out --strat_out
}

# Function for pathway-level inference
infer_pathway() {
  pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                      -o pathways_out -p ${num_threads}
}

# Function for adding functional descriptions
add_functional_descriptions() {
  add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

  add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                      -o pathways_out/path_abun_unstrat_descrip.tsv.gz
}

# Create a new directory and navigate into it
mkdir process
cd process

# Run each function in order
place_seqs
predict_gene_families
generate_metagenome_predictions
infer_pathway
add_functional_descriptions

