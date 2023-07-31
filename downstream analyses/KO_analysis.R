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
ko_ab <- read.table("filtered_KO_pred_metagenome_unstrat.tsv", header = TRUE, 
                    sep = "\t", row.names = NULL)


############################ Perform pathway DAA ############################ 
# Perform pathway DAA using edgeR method
daa_results_edgeR_df <- pathway_daa(abundance = ko_ab %>% column_to_rownames("function."), 
                                     metadata = metadata, group = "Group", 
                                     daa_method = "edgeR", select = NULL, reference = NULL)

# Perform pathway DAA using LinDA method
daa_results_LinDA_df <- pathway_daa(abundance = ko_ab %>% column_to_rownames("function."), 
                                    metadata = metadata, group = "Group", 
                                    daa_method = "LinDA", select = NULL, reference = NULL)

# Perform pathway DAA using limma voom method
daa_results_limma_voom_df <- pathway_daa(abundance = ko_ab %>% column_to_rownames("function."), 
                                    metadata = metadata, group = "Group", 
                                    daa_method = "limma voom", select = NULL, reference = NULL)


############################ Annotate pathway ############################ 
# Annotate pathway results without KO to KEGG conversion
daa_annotated_edgeR_method_results_df <- pathway_annotation(pathway = "KO", 
                                                                 daa_results_df = daa_results_edgeR_df, 
                                                                 ko_to_kegg = FALSE)
daa_annotated_LinDA_method_results_df <- pathway_annotation(pathway = "KO", 
                                                                 daa_results_df = daa_results_LinDA_df, 
                                                                 ko_to_kegg = FALSE)
daa_annotated_limma_voom_method_results_df <- pathway_annotation(pathway = "KO", 
                                                            daa_results_df = daa_results_limma_voom_df, 
                                                            ko_to_kegg = FALSE)

# write to tsv
write.table(daa_annotated_edgeR_method_results_df, 
            "KO DA/annotated_da_edgeR.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)

write.table(daa_annotated_LinDA_method_results_df, 
            "KO DA/annotated_da_LinDA.tsv", sep = "\t", row.names = FALSE, 
            quote = FALSE)

write.table(daa_annotated_limma_voom_method_results_df, 
            "KO DA/annotated_da_limma_voom.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)


############################ Generate pathway error bar plot ############################ 
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
# edgeR
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00007', 'K00028', 'K00039', 'K00064',
                                                            'K00078', 'K00090', 'K00094', 'K00098', 
                                                            'K00108', 'K00116', 'K00117', 'K00137', 
                                                            'K00138', 'K00141', 'K00146', 'K00151', 
                                                            'K00214', 'K00216', 'K00276', 'K00299', 
                                                            'K00316', 'K00322', 'K00363', 'K00367', 
                                                            'K00427', 'K00436', 'K00448', 'K00449', 
                                                            'K00453', 'K00455'), 
                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("KO DA/Error Bar/errorbar.edgeR1.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00457', 'K00464', 'K00480', 'K00481', 
                                                            'K00484', 'K00486', 'K00508', 'K00529', 
                                                            'K00673', 'K00689', 'K00692', 'K00839', 
                                                            'K00840', 'K00892', 'K00906', 'K00932', 
                                                            'K01031', 'K01032', 'K01055', 'K01083', 
                                                            'K01120', 'K01136', 'K01146', 'K01169', 
                                                            'K01227', 'K01252', 'K01390', 'K01407', 
                                                            'K01457', 'K01484'), p_value_bar = TRUE, 
                                                 colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.edgeR2.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K01521', 'K01545', 'K01602', 'K01631', 
                                                            'K01637', 'K01671', 'K01690', 'K01782', 
                                                            'K01825', 'K01826', 'K01856', 'K01857', 
                                                            'K01941', 'K02021', 'K02092', 'K02093',
                                                            'K02094', 'K02095', 'K02096', 'K02097', 
                                                            'K02100', 'K02167', 'K02255', 'K02284'),
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.edgeR3.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")



# limma voom
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_limma_voom_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select =c('K00058', 'K00133', 'K00134', 'K00151', 
                                                           'K00215', 'K00228', 'K00266', 'K00276', 
                                                           'K00384', 'K00455', 'K00457', 'K00484', 
                                                           'K00558', 'K00604', 'K00615', 'K00648', 
                                                           'K00655', 'K00688', 'K00758', 'K00789', 
                                                           'K00790', 'K00791', 'K00806', 'K00826', 
                                                           'K00845', 'K00873', 'K00876', 'K00919', 
                                                           'K00928', 'K00945'), p_value_bar = TRUE, 
                                                 colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.limma_voom1.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_limma_voom_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00951', 'K00962', 'K00981', 'K01091',
                                                            'K01104', 'K01358', 'K01409', 'K01448',
                                                            'K01462', 'K01479', 'K01537', 'K01591',
                                                            'K01624', 'K01652', 'K01710', 'K01733',
                                                            'K01738', 'K01740', 'K01752', 'K01756', 
                                                            'K01775', 'K01826', 'K01867', 'K01868', 
                                                            'K01872', 'K01875', 'K01876', 'K01889',
                                                            'K01939', 'K01955'), 
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.limma_voom2.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_limma_voom_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select =c('K01956', 'K01992', 'K02003', 'K02004',
                                                           'K02013', 'K02015', 'K02016', 'K02028',
                                                           'K02049', 'K02050', 'K02057', 'K02313',
                                                           'K02314', 'K02338', 'K02355', 'K02356',
                                                           'K02358', 'K02428', 'K02469', 'K02470',
                                                           'K02495', 'K02508', 'K02511', 'K02529'), 
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.limma_voom3.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")


# LinDA
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_LinDA_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00090', 'K00116', 'K00117', 'K00137',
                                                            'K00138', 'K00146', 'K00151', 'K00216',
                                                            'K00228', 'K00276', 'K00299', 'K00322',
                                                            'K00427', 'K00455', 'K00457', 'K00484',
                                                            'K00673', 'K00758', 'K00840', 'K00892',
                                                            'K00906', 'K00932', 'K01146', 'K01169',
                                                            'K01252', 'K01407', 'K01479', 'K01484',
                                                            'K01521', 'K01637'), 
                             p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("KO DA/Error Bar/errorbar.LinDA1.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_LinDA_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K01669', 'K01690', 'K01782', 'K01825',
                                                            'K01826', 'K02062', 'K02063', 'K02064',
                                                            'K02167', 'K02255', 'K02299', 'K02317',
                                                            'K02336', 'K02345', 'K02364', 'K02466',
                                                            'K02467', 'K02468', 'K02485', 'K02508',
                                                            'K02511', 'K02562', 'K02565', 'K02609',
                                                            'K02610', 'K02611', 'K02612', 'K02613',
                                                            'K02615', 'K02616'), 
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("KO DA/Error Bar/errorbar.LinDA2.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_LinDA_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K02617', 'K02618', 'K02846', 'K02853',
                                                            'K02972', 'K03078', 'K03087', 'K03112',
                                                            'K03181', 'K03184', 'K03207', 'K03208',
                                                            'K03214', 'K03304', 'K03425', 'K03435',
                                                            'K03468', 'K03472', 'K03477', 'K03485',
                                                            'K03535', 'K03580', 'K03668', 'K03712'), 
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("KO DA/Error Bar/errorbar.LinDA3.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")


# common pathways between edgeR, limma voom, LinDA
pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_edgeR_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00228', 'K00455', 'K00276', 'K00151',
                                                            'K00758', 'K00457', 'K00484'),
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.edgeR.common.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")

pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_limma_voom_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select =c('K00228', 'K00455', 'K00276', 'K00151',
                                                           'K00758', 'K00457', 'K00484'),
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")
ggsave("KO DA/Error Bar/errorbar.limma_voom.common.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")


pathway_errorbar <- ggpicrust2::pathway_errorbar(abundance = ko_ab %>% column_to_rownames("function."), 
                                                 daa_results_df = daa_annotated_LinDA_method_results_df, 
                                                 Group = metadata$Group, ko_to_kegg = FALSE, 
                                                 p_values_threshold = 0.05, order = "group", 
                                                 select = c('K00228', 'K00455', 'K00276', 'K00151',
                                                            'K00758', 'K00457', 'K00484'),
                                                 p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave("KO DA/Error Bar/errorbar.LinDA.common.svg", plot = pathway_errorbar, width = 20, height = 12, units = "in")



############################ Generate pathway heatmap ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_edgeR_with_p_0.05 <- daa_annotated_edgeR_method_results_df %>% filter(p_adjust < 0.05)
selected_rows_edgeR <- feature_edgeR_with_p_0.05[feature_edgeR_with_p_0.05$feature %in% 
                                                   c('K00228', 'K00455', 'K00276', 'K00151',
                                                     'K00758', 'K00457', 'K00484'),]

feature_limma_voom_with_p_0.05 <- daa_annotated_limma_voom_method_results_df %>% filter(p_adjust < 0.05)

feature_LinDA_with_p_0.05 <- daa_annotated_LinDA_method_results_df %>% filter(p_adjust < 0.05)


# common pathways
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ko_ab %>% 
                                                 filter(function. %in% selected_rows_edgeR$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("KO DA/pathway heatmap/heatmap.common.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# edgeR
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ko_ab %>% 
                                                 filter(function. %in% feature_edgeR_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("KO DA/pathway heatmap/heatmap.edgeR.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# limma_voom
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ko_ab %>% 
                                                 filter(function. %in% feature_limma_voom_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("KO DA/pathway heatmap/heatmap.limma_voom.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")

# LinDA
pathway_heatmap <- ggpicrust2::pathway_heatmap(abundance = ko_ab %>% 
                                                 filter(function. %in% feature_LinDA_with_p_0.05$feature) %>% 
                                                 column_to_rownames("function."),
                                               metadata = metadata, group = "Group")
ggsave("KO DA/pathway heatmap/heatmap.LinDA.svg", plot = pathway_heatmap, width = 20, height = 12, units = "in")


############################ Generate pathway PCA plot ############################ 
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset

# common pathways
pathway_pca <- ggpicrust2::pathway_pca(abundance = ko_ab %>% 
                                         filter(function. %in% selected_rows_edgeR$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("KO DA/PCA/pca.common.svg", plot = pathway_pca, width = 20, height = 12, units = "in")



# edgeR
pathway_pca <- ggpicrust2::pathway_pca(abundance = ko_ab %>% 
                                         filter(function. %in% feature_edgeR_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("KO DA/PCA/pca.edgeR.svg", plot = pathway_pca, width = 20, height = 12, units = "in")


# limma_voom
pathway_pca <- ggpicrust2::pathway_pca(abundance = ko_ab %>% 
                                         filter(function. %in% feature_limma_voom_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("KO DA/PCA/pca.limma_voom.svg", plot = pathway_pca, width = 20, height = 12, units = "in")


# LinDA
pathway_pca <- ggpicrust2::pathway_pca(abundance = ko_ab %>% 
                                         filter(function. %in% feature_LinDA_with_p_0.05$feature) %>% 
                                         column_to_rownames("function."), metadata = metadata, group = "Group")

ggsave("KO DA/PCA/pca.LinDA.svg", plot = pathway_pca, width = 20, height = 12, units = "in")
