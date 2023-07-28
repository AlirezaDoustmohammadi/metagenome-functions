# Data pre-processing
Depending on the pipeline you used to process your raw 16S data there may be many extremely rare ASVs found across your samples. 
Such rare ASVs typically only add noise to analyses (especially singletons - ASVs found only by 1 read in 1 sample) and should be removed. 

Similarly, it is important to check for any low-depth samples that should be removed. 
These can be identified based on the output of biom summarize-table.

<b>Note: </b> The best minimum cut-offs for excluding ASVs and samples varies by dataset since these cut-offs will differ depending on the overall read depth of your dataset.
