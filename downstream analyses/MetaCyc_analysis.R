library(readr)
library(ggpicrust2)
library(ggplot2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)



# Read metadata
metadata <- read.table("filtered_metadata.tsv", header = TRUE, sep = "\t", 
                       row.names = NULL)

# ko abundance
metacyc_ab <- read.table("filtered_MetaCyc_pred_metagenome_unstrat.tsv", header = TRUE, 
                    sep = "\t", row.names = NULL)


############################ Perform pathway DAA ############################ 
# Perform pathway DAA using edgeR method
daa_results_edgeR_df <- pathway_daa(abundance = metacyc_ab %>% column_to_rownames("pathway"), 
                                     metadata = metadata, group = "Group", 
                                     daa_method = "edgeR", select = NULL, reference = NULL)


############################ Annotate pathway ############################ 
# Annotate pathway results without KO to KEGG conversion
daa_annotated_edgeR_method_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                 daa_results_df = daa_results_edgeR_df, 
                                                                 ko_to_kegg = FALSE)

# write to tsv
write.table(daa_annotated_edgeR_method_results_df, 
            "MetaCyc DA/annotated_da_edgeR.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)


############################ Generate pathway error bar plot ############################ 
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
# edgeR
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = metacyc_ab %>% column_to_rownames("pathway"), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('3-HYDROXYPHENYLACETATE-DEGRADATION-PWY', 
                                                            'AEROBACTINSYN-PWY', 'ALL-CHORISMATE-PWY', 
                                                            'ARGDEG-PWY', 'AST-PWY', 'CATECHOL-ORTHO-CLEAVAGE-PWY',
                                                            'CHLOROPHYLL-SYN', 'CRNFORCAT-PWY', 'ECASYN-PWY', 
                                                            'ENTBACSYN-PWY', 'GALLATE-DEGRADATION-I-PWY', 
                                                            'GALLATE-DEGRADATION-II-PWY', 
                                                            'GLYCOLYSIS-TCA-GLYOX-BYPASS', 'HCAMHPDEG-PWY', 
                                                            'KETOGLUCONMET-PWY', 'METHYLGALLATE-DEGRADATION-PWY', 
                                                            'NADSYN-PWY', 'ORNARGDEG-PWY', 'ORNDEG-PWY', 'P105-PWY', 'PROTOCATECHUATE-ORTHO-CLEAVAGE-PWY', 'PWY-1541', 'PWY-1622', 'PWY-181', 'PWY-5088', 'PWY-5180', 'PWY-5181', 'PWY-5182', 'PWY-5415', 'PWY-5417'), 
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("MetaCyc DA/Error Bar/errorbar.edgeR1.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = metacyc_ab %>% column_to_rownames("pathway"), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('PWY-5431', 'PWY-5531', 'PWY-5651', 'PWY-6071', 
                                                            'PWY-6143', 'PWY-6182', 'PWY-6185', 'PWY-6562', 
                                                            'PWY-6690', 'PWY-7159', 'PWY-7373', 'PWY-7446', 
                                                            'PWY-7644', 'PWY0-1277', 'PWY0-321', 
                                                            'TCA-GLYOX-BYPASS', 'THREOCAT-PWY', 'TYRFUMCAT-PWY', 
                                                            'VALDEG-PWY'), p_value_bar = TRUE, 
                                                 colors = NULL, x_lab = "description")

ggsave("MetaCyc DA/Error Bar/errorbar.edgeR2.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")


############################ Generate pathway heatmap ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_edgeR_with_p_0.05 <- daa_annotated_edgeR_method_results_df %>% filter(p_adjust < 0.05)

# edgeR
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = metacyc_ab %>% 
                                                 filter(pathway %in% feature_edgeR_with_p_0.05$feature) %>% 
                                                 column_to_rownames("pathway"),
                                               metadata = metadata, group = "Group")
ggsave("MetaCyc DA/pathway heatmap/heatmap.edgeR.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

############################ Generate pathway PCA plot ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset

# edgeR
pathway_pca <- ggpicrust2::pathway_pca(abundance = metacyc_ab %>% 
                                         filter(pathway %in% feature_edgeR_with_p_0.05$feature) %>% 
                                         column_to_rownames("pathway"), metadata = metadata, group = "Group")

ggsave("MetaCyc DA/PCA/pca.edgeR.svg", plot = pathway_pca, width = 20, height = 12, units = "in")
