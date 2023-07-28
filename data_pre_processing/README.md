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
1. converting feature table file & ASVs sequence file to qza format (readable for QIIME 2)
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
              --input-path rwads.fasta \
              --output-path reads.qza
   ```
   
