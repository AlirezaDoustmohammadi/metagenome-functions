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
ec_ab <- read.table("filtered_EC_pred_metagenome_unstrat.tsv", header = TRUE, 
                    sep = "\t", row.names = NULL)


############################ Perform pathway DAA ############################ 
# Perform pathway DAA using edgeR method
daa_results_edgeR_df <- pathway_daa(abundance = ec_ab %>% column_to_rownames("function."), 
                                     metadata = metadata, group = "Group", 
                                     daa_method = "edgeR", select = NULL, reference = NULL)

# Perform pathway DAA using LinDA method
daa_results_LinDA_df <- pathway_daa(abundance = ec_ab %>% column_to_rownames("function."), 
                                    metadata = metadata, group = "Group", 
                                    daa_method = "LinDA", select = NULL, reference = NULL)

# Perform pathway DAA using limma voom method
daa_results_limma_voom_df <- pathway_daa(abundance = ec_ab %>% column_to_rownames("function."), 
                                    metadata = metadata, group = "Group", 
                                    daa_method = "limma voom", select = NULL, reference = NULL)


############################ Annotate pathway ############################ 
# Annotate pathway results without KO to KEGG conversion
daa_annotated_edgeR_method_results_df <- pathway_annotation(pathway = "EC", 
                                                                 daa_results_df = daa_results_edgeR_df, 
                                                                 ko_to_kegg = FALSE)
daa_annotated_LinDA_method_results_df <- pathway_annotation(pathway = "EC", 
                                                                 daa_results_df = daa_results_LinDA_df, 
                                                                 ko_to_kegg = FALSE)
daa_annotated_limma_voom_method_results_df <- pathway_annotation(pathway = "EC", 
                                                            daa_results_df = daa_results_limma_voom_df, 
                                                            ko_to_kegg = FALSE)

# write to tsv
write.table(daa_annotated_edgeR_method_results_df, 
            "EC DA/annotated_da_edgeR.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)

write.table(daa_annotated_LinDA_method_results_df, 
            "EC DA/annotated_da_LinDA.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)

write.table(daa_annotated_limma_voom_method_results_df, 
            "EC DA/annotated_da_limma_voom.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)


############################ Generate pathway error bar plot ############################ 
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
# edgeR
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ec_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('EC:1.1.1.11', 'EC:1.1.1.122', 'EC:1.1.1.140',
                                                            'EC:1.1.1.179', 'EC:1.1.1.215', 'EC:1.1.1.220',
                                                            'EC:1.1.1.251', 'EC:1.1.1.264', 'EC:1.1.1.291', 
                                                            'EC:1.1.1.310', 'EC:1.1.1.313', 'EC:1.1.1.333', 
                                                            'EC:1.1.1.339', 'EC:1.1.1.373', 'EC:1.1.1.39',
                                                            'EC:1.1.1.56', 'EC:1.1.2.3', 'EC:1.1.5.2', 
                                                            'EC:1.1.5.4', 'EC:1.1.5.6', 'EC:1.1.5.8', 
                                                            'EC:1.1.99.1', 'EC:1.1.99.3', 'EC:1.10.9.1',
                                                            'EC:1.11.1.1', 'EC:1.11.1.21', 'EC:1.11.2.4',
                                                            'EC:1.12.5.1', 'EC:1.13.11.1', 'EC:1.13.11.11'), 
                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("EC DA/Error Bar/errorbar.edgeR1.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ec_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('EC:1.13.11.15', 'EC:1.13.11.16', 'EC:1.13.11.27',
                                                            'EC:1.13.11.3', 'EC:1.13.11.5', 'EC:1.13.11.53',
                                                            'EC:1.13.11.54', 'EC:1.13.11.57', 'EC:1.13.11.75',
                                                            'EC:1.14.11.1', 'EC:1.14.11.17', 'EC:1.14.11.33',
                                                            'EC:1.14.11.47', 'EC:1.14.12.10', 'EC:1.14.12.19',
                                                            'EC:1.14.13.1', 'EC:1.14.13.113', 'EC:1.14.13.127',
                                                            'EC:1.14.13.129', 'EC:1.14.13.149', 'EC:1.14.13.2',
                                                            'EC:1.14.13.59'),
                                                 p_value_bar = TRUE, 
                                                 colors = NULL, x_lab = "description")
ggsave("EC DA/Error Bar/errorbar.edgeR2.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

# limma voom
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ec_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_limma_voom_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select =NULL, p_value_bar = TRUE, 
                                                 colors = NULL, x_lab = "description")
ggsave("EC DA/Error Bar/errorbar.limma_voom.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")


# LinDA
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ec_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_LinDA_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = NULL, 
                             p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("EC DA/Error Bar/errorbar.LinDA.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

############################ Generate pathway heatmap ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_edgeR_with_p_0.05 <- daa_annotated_edgeR_method_results_df %>% filter(p_adjust < 0.05)
selected_rows_edgeR <- feature_edgeR_with_p_0.05[feature_edgeR_with_p_0.05$feature %in% 
                                                   c('EC:1.11.1.21'),]

feature_limma_voom_with_p_0.05 <- daa_annotated_limma_voom_method_results_df %>% filter(p_adjust < 0.05)

feature_LinDA_with_p_0.05 <- daa_annotated_LinDA_method_results_df %>% filter(p_adjust < 0.05)


# common pathways
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ec_ab %>% 
                                                 filter(function. %in% selected_rows_edgeR$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("EC DA/pathway heatmap/heatmap.common.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# edgeR
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ec_ab %>% 
                                                 filter(function. %in% feature_edgeR_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("EC DA/pathway heatmap/heatmap.edgeR.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# limma_voom
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ec_ab %>% 
                                                 filter(function. %in% feature_limma_voom_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("EC DA/pathway heatmap/heatmap.limma_voom.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# LinDA
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ec_ab %>% 
                                                 filter(function. %in% feature_LinDA_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("EC DA/pathway heatmap/heatmap.LinDA.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")


############################ Generate pathway PCA plot ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset

# edgeR
pathway_pca <- ggpicrust2::pathway_pca(abundance = ec_ab %>% 
                                         filter(function. %in% feature_edgeR_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("EC DA/PCA/pca.edgeR.svg", plot = pathway_pca, width = 20, height = 12, units = "in")


# limma_voom
pathway_pca <- ggpicrust2::pathway_pca(abundance = ec_ab %>% 
                                         filter(function. %in% feature_limma_voom_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("EC DA/PCA/pca.limma_voom.svg", plot = pathway_pca, width = 20, height = 12, units = "in")


# LinDA
pathway_pca <- ggpicrust2::pathway_pca(abundance = ec_ab %>% 
                                         filter(function. %in% feature_LinDA_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("EC DA/PCA/pca.LinDA.svg", plot = pathway_pca, width = 20, height = 12, units = "in")
