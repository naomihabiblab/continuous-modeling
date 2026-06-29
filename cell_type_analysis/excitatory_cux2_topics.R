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

DATA_PATH <- file.path(get_result_dir(), "excitatory_topics")

cux <- Sys.getenv(x = 'CUX', unset='plus')
if (!(cux %in% c("plus", "minus"))) {
  stop(paste("Wrong cux option", cux))
}

topic_output_name_prefix <- paste0("all_cux2_", cux, "_excitatory")

cux_file_suffix <- if(cux=="plus") '+' else '-'
neural_class <- paste("cux2", cux, sep = "_")
cell_type <- paste0('cux2', if(cux=="plus") 'p' else 'm')


plots_dir <- file.path(get_result_dir(), sprintf('ExcitatoryCux2%sAll', cux_file_suffix))
dir.create(plots_dir, showWarnings = FALSE)

object <- LoadH5Seurat(file.path(SEURAT_OBJECT_DATA_PATH, sprintf("cux2%s.h5Seurat", cux_file_suffix)))
object <- enrich_meta_data(object)
object <- fix_similarity_names(object)

pred_label <- readRDS(sprintf("../data/preds_cux%s.rds", cux_file_suffix))
pred_sublclass <- readRDS(sprintf("../data/preds_cux%s_subclass.rds", cux_file_suffix))
sublcass_hclust <- as.hclust(readRDS("../data/human_m1_10x_human_dendrogram.rds"))
sublcass_label_order <- c(sublcass_hclust$labels[sublcass_hclust$order], 'NA')

pred_label$pruned.labels[is.na(pred_label$pruned.labels)] <- 'NA'
pred_sublclass$pruned.labels[is.na(pred_sublclass$pruned.labels)] <- 'NA'
pred_sublclass$pruned.labels <- factor(pred_sublclass$pruned.labels, levels=sublcass_label_order)

pred_label <- pred_label[rownames(object@meta.data),]
pred_sublclass <- pred_sublclass[rownames(object@meta.data),]

object <- AddMetaData(object, pred_label$labels, "predicted_subtype")
object <- AddMetaData(object, pred_label$pruned.labels, "pruned_subtype")

object <- AddMetaData(object, pred_sublclass$labels, "predicted_subclass")
object <- AddMetaData(object, pred_sublclass$pruned.labels, "pruned_subclass")


# load fits
fits <- read_topic_models(file.path(DATA_PATH, sprintf("all_%s_excitatory_topic_fit", neural_class)))


print("Adding topics to metadata...")
topcis_and_cols <- add_topics_to_metadata(object, fits)
object <- topcis_and_cols$obj
topics_columns_list <- topcis_and_cols$columns_list
rm(topcis_and_cols)
gc()

# Specfic fit -------------------------------------------------------------

k <- if(cux=="plus") '14' else '16'
save_umap_densities(object, topics_columns_list[[k]])

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
df <- df_pos_des(obj=object, reweighted_f, fits[[k]], assay="SCT", organism='hsa')
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
save_allen_annotation(obj = object, cell.type = cell_type)
save_abundance_stat(inhibitory, group = CELL_TYPE_TO_GROUP[[cell_type]], topics_columns=topics_columns_list[[k]])


##########
#   Compute and save programs signature score 
##########
genes_mapping <- get_genes_by_topic(cell_type)

signature_score_prefix <- "Signature_"
object <- add_module_score_named(object, genes_mapping$genes_by_topic, prefix=signature_score_prefix)

signature_unique_score_prefix <- "SignatureUnique_"
object <- add_module_score_named(object, genes_mapping$unique_genes_by_topic, prefix=signature_unique_score_prefix)

h5write(object@meta.data %>% dplyr::select(starts_with("Signature")), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/Signature"))
h5write(object@meta.data %>% dplyr::select(starts_with("SignatureUnique")), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/SignatureUnique"))

