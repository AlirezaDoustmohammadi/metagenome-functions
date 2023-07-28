# Data pre-processing
Depending on the pipeline you used to process your raw 16S data there may be many <b>extremely rare ASVs</b> found across your samples. 
Such rare ASVs typically only add noise to analyses (especially singletons - ASVs found only by 1 read in 1 sample) and should be removed. 

Similarly, it is important to check for any <b>low-depth samples</b> that should be removed. 
These can be identified based on the output of <I>biom summarize-table</i>.

<b>Note: </b> The best minimum cut-offs for excluding ASVs and samples varies by dataset since these cut-offs will differ depending on the overall read depth of your dataset.

## How can I filter the feature table (Biom) file?
An easy way is to use <a href="https://docs.qiime2.org/2023.5/">QIIME 2</a></br>
QIIME 2 is a powerful, extensible, and decentralized microbiome analysis package with a focus on data and analysis transparency.

### Installing QIIME 2
```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-linux-conda.yml
```

### Filtering data
1. converting feature table file & ASVs sequence file to qza format (readable for QIIME 2)</br>
   a. converting BIOM to qza:
   ```
        qiime tools import \
              --input-path feature_table.biom \
              --type 'FeatureTable[Frequency]' \
              --input-format BIOMV210Format \
              --output-path feature_table.qza
   ```
   b. converting fna/fasta to qza:
   ```
       qiime tools import \
              --type 'FeatureData[Sequence]' \
              --input-path reads.fasta \
              --output-path reads.qza
   ```
2. Total-frequency-based filtering</br>
   a. filtering <b>features</b> based on how frequently they are represented in the feature table
   ```
         qiime feature-table filter-features \
                 --i-table feature_table.qza \
                 --p-min-frequency 10 \
                 --o-filtered-table feature_table.total.frequency.filtered.features.qza
   ```

   b. filtering <b>samples</b> whose total frequency is less than threshold
   ```
         qiime feature-table filter-samples \
                 --i-table feature_table.total.frequency.filtered.features.qza \
                 --p-min-frequency 1500 \
                 --o-filtered-table feature_table.total.frequency.filtered.features.samples.qza
   ```
3. Contingency-based filtering
   a. filtering features from a table contingent on the number of samples they’re observed in
   ```
         qiime feature-table filter-features \
                 --i-table feature_table.total.frequency.filtered.features.samples.qza \
                 --p-min-samples 2 \
                 --o-filtered-table feature_table.contingency.filtered.features.qza
   ```
   b. filtering samples from a table contingent on the number of features they contain
   ```
         qiime feature-table filter-samples \
                 --i-table feature_table.contingency.filtered.features.qza \
                 --p-min-features 10 \
                 --o-filtered-table feature_table.contingency.filtered.features.samples.qza
   ```

<b>Note: </b> All of these methods can also be applied to filter on the maximum number of features or samples, using the `--p-max-features` and `--p-max-samples` parameters</br>
Other types of filtration:
1. Identifier-based filtering: </br>
   Identifier-based filtering is used to retain only a user-specified list of samples or features
2. Metadata-based filtering: </br>
   Metadata-based filtering is similar to identifier-based filtering, except that the list of IDs to keep is determined based on metadata search criteria rather than being provided by the user directly.
3. Taxonomy-based filtering

   
