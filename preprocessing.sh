#!/bin/bash

# start with `conda activate qiime2-2023.5`

# Declare variables for biom and sequence files, and minimum abundances
raw_biom='feature-table.biom'
raw_seq='representative-sequences.fasta'
min_feature_abundance=15
min_samples_abundance=1500

# Extract filenames without extensions
raw_biom_filename="${raw_biom%.biom}"
raw_seq_filename="${raw_seq%.fna}"
raw_seq_filename="${raw_seq_filename%.fasta}"

convert_biom_to_tsv_and_summarize() {
  mkdir raw.summary
  cd raw.summary
  biom convert -i ../${raw_biom} -o ${raw_biom_filename}.tsv --to-tsv --header-key taxonomy
  # summariese the biome file
  biom summarize-table -i ../${raw_biom} -o ${raw_biom_filename}.summary
  python ../biomSummary.py --biom_file ${raw_biom}
  cd ..
}

convert_files_to_qza_format() {
  # Convert feature tables from biom files to QZA-format for filtering feature table using qiime
  qiime tools import \
    --input-path ${raw_biom} \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path ${raw_biom_filename}.qza

  # Convert ASVs sequence file from fasta files to QZA-format for filtering them using qiime
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path ${raw_seq} \
    --output-path ${raw_seq_filename}.qza
}

filter_features_and_samples_total_frequency_based() {

  # filtering features based on how frequently they are represented in the feature table (biom file).
  mkdir data_pre_processing
  cd data_pre_processing

  # remove all features with a total abundance (summed across all samples) of less than threshold
  qiime feature-table filter-features \
    --i-table ../${raw_biom_filename}.qza \
    --p-min-frequency ${min_feature_abundance} \
    --o-filtered-table ${raw_biom_filename}.total.frequency.filtered.features.qza

  # filtering samples whose total frequency is less than threshold
  qiime feature-table filter-samples \
    --i-table ${raw_biom_filename}.filtered.features.qza \
    --p-min-frequency ${min_samples_abundance} \
    --o-filtered-table ${raw_biom_filename}.total.frequency.filtered.features.samples.qza
}

filter_sequences_and_export_to_formats() {

  # filter sequences to match the filtered feature table 
  qiime feature-table filter-seqs \
    --i-data ../${raw_seq_filename}.qza \
    --i-table ${raw_biom_filename}.total.frequency.filtered.features.samples.qza \
    --o-filtered-data ${raw_seq_filename}.filtered.features.samples.qza

  # convert filtered feature table & ASVs sequence to biom & fasta format
  # convert filtered ASVs sequence file to fasta format
  qiime tools export \
    --input-path ${raw_seq_filename}.filtered.features.samples.qza \
    --output-path ../data_pre_processing
  mv dna-sequences.fasta ${raw_seq_filename}.filtered.features.samples.fasta

  # convert feature table to biom format
  qiime tools export \
    --input-path ${raw_biom_filename}.total.frequency.filtered.features.samples.qza \
    --output-path ../data_pre_processing/${raw_biom_filename}.total.frequency.filtered.features.samples.biom \
    --output-format BIOMV210Format
}

convert_filtered_biom_to_tsv_and_summarize() {
  mkdir summary
  cd summary
  # convert filtered biom file to tsv
  biom convert -i ../${raw_biom_filename}.filtered.features.samples.biom -o ${raw_biom_filename}.total.frequency.filtered.features.samples.tsv \
				--to-tsv --header-key taxonomy
				
  # summariese the biome file
  biom summarize-table -i ../${raw_biom_filename}.total.frequency.filtered.features.samples.biom -o ${raw_biom_filename}.total.frequency.filtered.features.samples.summary
  python ../../biomSummary.py --biom_file ${raw_biom_filename}.total.frequency.filtered.features.samples.biom
  
}

# Call each function in order
convert_biom_to_tsv_and_summarize
convert_files_to_qza_format
filter_features_and_samples_total_frequency_based
filter_sequences_and_export_to_formats
convert_filtered_biom_to_tsv_and_summarize
