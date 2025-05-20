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

source("utils.R")
source("pathways.analysis.R")
source("landscape.analysis.funs.R")
source("topic_space_functions.R")
source("figures.funs.R")
source("markers.R")
source("matrix_utils.R")
source("data_500_utils.R")
source("topics.R")
source("new_method_topics_functions.R")

cell_type <- 'opcs'
DATA_PATH <- file.path(get_result_dir(), "opcs_topics")

k <- Sys.getenv(x = 'SLURM_ARRAY_TASK_ID', unset = '7')

obj <- LoadH5Seurat(file.path(SEURAT_OBJECT_DATA_PATH, "opcs.h5Seurat"))
obj <- enrich_meta_data(obj)
obj <- fix_similarity_names(obj)


plots_dir <- file.path(get_result_dir(), 'OPCAll')
dir.create(plots_dir, showWarnings = FALSE)

# load fits
fits <- read_topic_models(file.path(DATA_PATH, "all_opcs_topic_fit"))



print("Adding topics to metadata...")
topcis_and_cols <- add_topics_to_metadata(obj, fits)
obj <- topcis_and_cols$obj
topics_columns_list <- topcis_and_cols$columns_list
rm(topcis_and_cols)
gc()

k <- '7'
save_umap_densities(obj, topics_columns_list[[k]], cell.type = cell_type)


print(paste0("Analyzing k=", k))
plots_dir_k <- file.path(plots_dir, k)
dir.create(plots_dir_k, showWarnings = FALSE)



# DE analysis  ------------------------------------------------------------

print("DE analysis")
dfas_vs_null <- read_de(prefix = file.path(DATA_PATH, "all_opcs_dfa_vsnull"))


fast_topics_de <- get_top_de_genes(dfas_vs_null[[k]])

filter_de <- function (de, fast_topics_de) {
  de %>%
    dplyr::left_join(fast_topics_de, by=c("topic", "gene")) %>%
    filter(!is.na(lfsr)) %>%
    filter(lfsr < 0.05) %>%
    filter(postmean > 0) %>%
    dplyr::select(-z)
}

upper_quantile <- 0.8
upper_quantile_str <- upper_quantile %>% sprintf("%0.2f", .) %>% gsub("\\.", "", .)

reweighted_f <- topic_reweight_f(fits[[k]]$F, upper_quantile=upper_quantile)

message("Computing de...")
df <- df_pos_des(obj=obj, reweighted_f, fits[[k]], assay="SCT", organism='hsa')
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
save_abundance_stat(obj, group = CELL_TYPE_TO_GROUP[[cell_type]], topics_columns=topics_columns_list[[k]])
