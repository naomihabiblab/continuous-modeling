library(tibble)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(rhdf5)
library(future)
library(Matrix)
library(Seurat)
library(msigdbr)
library(SeuratDisk)
library(fastTopics)
library(clusterProfiler)
library(corrplot)
library(ggcorrplot)

source("R/utils.R")
source("R/pathways.analysis.R")
source("R/landscape.analysis.funs.R")
source("R/topic_space_functions.R")
source("R/figures.funs.R")
source("R/markers.R")
source("R/matrix_utils.R")
source("R/data_500_utils.R")
source("R/topics.R")
source("R/program_signature.R")
source("R/markers.funs.R")
source("validation/validation.functions.R")


cell_type <- 'inhibitory'
DATA_PATH <- file.path(get_result_dir(), "inhibitory_topics")

k <- Sys.getenv(x = 'SLURM_ARRAY_TASK_ID', unset = '17')


topic_output_name_prefix <- "all_inhibitory"

plots_dir <- file.path(get_result_dir(), 'InhibitoryAll')
dir.create(plots_dir, showWarnings = FALSE)

inhibitory <- LoadH5Seurat(file.path(SEURAT_OBJECT_DATA_PATH, "inhibitory.h5Seurat"))
inhibitory <- enrich_meta_data(inhibitory)
inhibitory <- fix_similarity_names(inhibitory)

pred_inh <- readRDS("../data/preds_inh.rds")
pred_inh_sublclass <- readRDS("../data/preds_inh_subclass.rds")
sublcass_hclust <- as.hclust(readRDS("../data/human_m1_10x_human_dendrogram.rds"))
sublcass_label_order <- c(sublcass_hclust$labels[sublcass_hclust$order], 'NA')

pred_inh$pruned.labels[is.na(pred_inh$pruned.labels)] <- 'NA'
pred_inh_sublclass$pruned.labels[is.na(pred_inh_sublclass$pruned.labels)] <- 'NA'
pred_inh_sublclass$pruned.labels <- factor(pred_inh_sublclass$pruned.labels, levels=sublcass_label_order)

pred_inh <- pred_inh[rownames(inhibitory@meta.data),]
pred_inh_sublclass <- pred_inh_sublclass[rownames(inhibitory@meta.data),]

inhibitory <- AddMetaData(inhibitory, pred_inh$labels, "predicted_subtype")
inhibitory <- AddMetaData(inhibitory, pred_inh$pruned.labels, "pruned_subtype")

inhibitory <- AddMetaData(inhibitory, pred_inh_sublclass$labels, "predicted_subclass")
inhibitory <- AddMetaData(inhibitory, pred_inh_sublclass$pruned.labels, "pruned_subclass")


marker_analysis(inhibitory, plots_dir, markers=neuronal_markers, group.by=c("state", "annotation"), markers_name="Neuronal")

fits <- read_topic_models(file.path(DATA_PATH, sprintf("%s_topic_fit", topic_output_name_prefix)))

print("Adding topics to metadata...")
topcis_and_cols <- add_topics_to_metadata(inhibitory, fits)
inhibitory <- topcis_and_cols$obj
topics_columns_list <- topcis_and_cols$columns_list
rm(topcis_and_cols)
gc()



k <- '17'

save_umap_densities(inhibitory, topics_columns_list[[k]])

print(paste0("Analyzing k=", k))
plots_dir_k <- file.path(plots_dir, k)
dir.create(plots_dir_k, showWarnings = FALSE)


# DE analysis  ------------------------------------------------------------
print("DE analysis")
dfas_vs_null <- read_de(prefix = file.path(DATA_PATH, sprintf("%s_dfa_vsnull", topic_output_name_prefix)))

fast_topics_de <- get_top_de_genes(dfas_vs_null[[k]])

filter_de <- function (de, fast_topics_de) {
  de %>%
    dplyr::left_join(fast_topics_de, by=c("topic", "gene")) %>%
    filter(!is.na(lfsr)) %>%
    filter(lfsr < 0.05) %>%
    filter(postmean > 0) %>%
    dplyr::select(-z)
}

upper_quantile <- 0.95
upper_quantile_str <- upper_quantile %>% sprintf("%0.2f", .) %>% gsub("\\.", "", .)

reweighted_f <- topic_reweight_f(fits[[k]]$F, upper_quantile=upper_quantile)

message("Computing de...")
df <- df_pos_des(obj=inhibitory, reweighted_f, fits[[k]], assay="SCT", organism='hsa')
de <- df %>% filter(z_score_log > 1)

filtered_de <- filter_de(de, fast_topics_de=fast_topics_de)

universe <- GeneIdMapping()$ids[rownames(fits[[k]]$F)]

pathways <- topic_gene_ontology(filtered_de, universe=universe)

dataset_path <- paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/de")
filtered_de$cell.type <- cell_type

attr(filtered_de, "upper_quantile") <- upper_quantile_str

print(glue::glue("Writing de of {cell_type} to h5 {dataset_path}"))
h5write(filtered_de, TOPICS_DATA_H5_PATH, dataset_path, write.attributes=TRUE)

dataset_path <- paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/pathways")
pathways$cell.type <- cell_type

attr(pathways, "upper_quantile") <- upper_quantile_str

print(glue::glue("Writing pathways of {cell_type} to h5 {dataset_path}"))
h5write(pathways, h5_path, dataset_path, write.attributes=TRUE)

##########
#   H5   #
##########
save_topics_score(cell_types = cell_type)
save_allen_annotation(obj = inhibitory, cell.type = cell_type)
save_abundance_stat(inhibitory, group = CELL_TYPE_TO_GROUP[[cell_type]], topics_columns=topics_columns_list[[k]])


##########
#   Within-cell permutation null (increasing programs k10 vs k8)
##########
# Self-contained: reads L from H5 (written above by save_topics_score), so this
# block can be moved/run independently of the rest of the pipeline.
source("R/topic_loading_permutation.R")

run_increasing_programs_permutation_null(
  cell.type = cell_type,
  output_dir = file.path(plots_dir, "permutation_null")
)


## Save subset fits to h5
subsets <- c("SST", "VIP", "PVALB")

for (subset_name in subsets) {
  fits <- read_topic_models(file.path(DATA_PATH, glue::glue("all_{subset_name}_inhibitory_topic_fit")), k=9)
  dataset_path <- glue::glue("{CELL_TYPE_TO_GROUP['inhibitory']}/subtype_{subset_name}")
  save_topic_fit(fit=fits[['9']], dataset_path = dataset_path)
}


### Senescence like genes signatures
senescence_like_signatures <- list(senescence_increasing=c("CHEK2", "ATM", "PIK3CA", "AKT3", "FOXO1", "BANP"),
                                   senescence_decreasing=c("BCL2", "BBC3"), 
                                   `NFKB pathway`=c("TRAF6", "NFKB1"),
                                   senescence_unknown=c("CDK2", "NBN"),
                                   cell_cycle=c("CENPC", "ANAPC16", "CEP290", "STAG2"),
                                   checkpoints=c("ATR", "RAD9A", "BRCA1"),
                                   translision=c("REV1", "POLK"),
                                   canonical_senescence=c("TP53", "CDKN1A", "CDKN2A"))

senescence_like_genes <- unlist(senescence_like_signatures, use.names = FALSE)


inhibitory <- AddModuleScore(inhibitory, senescence_like_signatures)
column_names <- paste0("Cluster", 1:length(senescence_like_signatures))
names(column_names) <- names(senescence_like_signatures)
inhibitory@meta.data <- inhibitory@meta.data %>% dplyr::rename(!!!column_names)

df <- FetchData(object = inhibitory, c(senescence_like_genes, 'projid', names(senescence_like_signatures)))
df <- df %>% rownames_to_column('cell')
h5delete(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/senescence_like_genes.signatures"))
h5write(df, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/senescence_like_genes.signatures"))



##########
#   Compute and save programs signature score 
##########
genes_mapping <- get_genes_by_topic(cell_type)

signature_score_prefix <- "Signature_"
inhibitory <- add_module_score_named(inhibitory, genes_mapping$genes_by_topic, prefix=signature_score_prefix)

signature_unique_score_prefix <- "SignatureUnique_"
inhibitory <- add_module_score_named(inhibitory, genes_mapping$unique_genes_by_topic, prefix=signature_unique_score_prefix)

h5write(inhibitory@meta.data %>% dplyr::select(starts_with("Signature")), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/Signature"))
h5write(inhibitory@meta.data %>% dplyr::select(starts_with("SignatureUnique")), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/SignatureUnique"))
