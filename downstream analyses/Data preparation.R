
# Read metadata
metadata <- read.csv("../grouping.csv")

# Remove rows where Group is NA
metadata <- metadata[!is.na(metadata$Group), ]

# EC abundance
ec_ab <- read.table(gzfile("../process/EC_metagenome_out//pred_metagenome_unstrat.tsv.gz"), 
                    header = TRUE, sep = "\t", row.names = NULL)


# KO abundance
ko_ab <- read.table(gzfile("../process/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"), 
                      header = TRUE, sep = "\t", row.names = NULL)
# MetaCyc pathway abundance
MetaCyc_ab <- read.table(gzfile("../process/pathways_out/path_abun_unstrat.tsv.gz"), 
                         header = TRUE, sep = "\t", row.names = NULL)


# removing samples not in metadata (samples with NA label)
keep_samples <- metadata$SampleID

# Subset data to only the columns you want to keep
ko_ab <- ko_ab[, colnames(ko_ab) %in% c("function.", keep_samples)]
ec_ab <- ec_ab[, colnames(ec_ab) %in% c("function.", keep_samples)]
MetaCyc_ab <- MetaCyc_ab[, colnames(MetaCyc_ab) %in% c("pathway", keep_samples)]

# write files in tsv format
# metadata
write.table(metadata, "filtered_metadata.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)
# KO abundance
write.table(ko_ab, "filtered_KO_pred_metagenome_unstrat.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)

# EC abundance
write.table(ec_ab, "filtered_EC_pred_metagenome_unstrat.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)

# MetaCyc pathway abundance
write.table(MetaCyc_ab, "filtered_MetaCyc_pred_metagenome_unstrat.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)

