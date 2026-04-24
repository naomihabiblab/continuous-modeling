#' @title Topic Model Analysis Functions
#' @description Comprehensive functions for fitting, analyzing, and visualizing
#'              topic models from single-cell data, including DE analysis,
#'              trait associations, and gene ontology enrichment.

library(rlang)
library(stringr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ltm)
library(ggpointdensity)


source("R/matrix_utils.R")
source("R/visualizations.R")
source("R/pathways.analysis.R")
source("R/landscape.analysis.funs.R")
source("R/figures.funs.R")
source("R/utils.R")
source("R/data_500_utils.R")

# Read topics -----------------------------------------------------------------

#' Extract model name from file path
#'
#' @param file_name File path string
#'
#' @return Model name extracted from filename (second-to-last dot-separated part)
extract_model_name_from_file_name <- function(file_name) {
  model_name <- file_name %>%
    strsplit('\\.') %>%
    dplyr::last() %>%
    tail(2) %>%
    head(1)
}

#' Read topic models from the disk by a pattern 
#' 
#' @return  a named list of topics models 
#' 
read_topic_models <- function(prefix, order_numeric = TRUE, k=NULL) {
  files <- Sys.glob(paste0(prefix, '\\.*', '.RDS'))
  topic_models <- list()
  print(files)
  for (file in files) {
    model_name <- extract_model_name_from_file_name(file)

    if (!is.null(k) & !(model_name %in% k)) {
      next
    }

    model <- readRDS(file)
    topic_models[[model_name]] <- model
  }
  
  models_names <- names(topic_models)
  if (order_numeric) {
    models_names = as.numeric(models_names)
  }
  topic_models <- topic_models[order(models_names)]
}

#' Read topic model likelihoods from disk
#'
#' @param prefix File path prefix (glob pattern)
#' @param order_numeric Sort models numerically by name (default: TRUE)
#'
#' @return Named list of likelihood objects
read_topic_likelihood <- function(prefix, order_numeric = TRUE) {
  files <- Sys.glob(paste0(prefix, '\\.*', '.RDS'))
  likelihood <- list()
  print(files)
  for (file in files) {
    model <- readRDS(file)
    model_name <- extract_model_name_from_file_name(file)
    likelihood[[model_name]] <- model
  }

  models_names <- names(likelihood)
  if (order_numeric) {
    models_names = as.numeric(models_names)
  }
  likelihood <- likelihood[order(models_names)]
}


#' Add topic loadings to Seurat metadata
#'
#' @param obj Seurat object
#' @param fits Named list of topic model fits
#'
#' @return List with modified obj and columns_list
add_topics_to_metadata <- function(obj, fits) {
  columns_list <- lapply(seq_along(fits), FUN = function(fit_index) {
    1:ncol(fits[[fit_index]]$L) %>% paste0(names(fits)[[fit_index]], ".k", .) %>% make.names()
  })
  names(columns_list) <- names(fits)
  for (fit_index in seq_along(fits)) {
    data <- data.frame(fits[[fit_index]]$L)
    cols <- columns_list[[fit_index]]
    obj <- AddMetaData(obj, metadata = data, col.name = cols)
  }

  list(obj = obj, columns_list = columns_list)
}


# DE analysis --------------------------------------------------------------

#' Read differential expression results from disk
#'
#' @param prefix File path prefix (glob pattern)
#' @param order_numeric Sort by name numerically (default: TRUE)
#'
#' @return Named list of DE analysis objects
read_de <- function(prefix, order_numeric = TRUE) {
  files <- Sys.glob(paste0(prefix, '\\.*', '.RDS'))
  des <- list()
  print(files)
  for (file in files) {
    de <- readRDS(file)
    de_name <- extract_model_name_from_file_name(file)
    des[[de_name]] <- de
  }
  
  de_names <- names(des)
  if (order_numeric) {
    de_names <- as.numeric(de_names)
  }
  des <- des[order(de_names)]
}

#' Run DE analysis on multiple topic fits
#'
#' @param fits List of topic model fits
#' @param counts_transpose Transposed count matrix
#' @param nc Number of cores (default: 6)
#' @param ... Additional arguments passed to de_analysis
#'
#' @return Named list of DE analysis results
topics_de_analysis <- function(fits, counts_transpose, nc = 6, ...) {
  sapply(fits, USE.NAMES = TRUE, simplify = FALSE, FUN = function(fit) {
    de_analysis(fit = fit, X = counts_transpose, control = list(nc = nc), ...)
  })
}

#' Get top differentially expressed genes from DE analysis
#'
#' @param dfa DE analysis object with postmean, lfsr, z matrices
#' @param number_of_genes Max genes per topic (-1 for all) (default: -1)
#' @param lfsr_threshold LFSR threshold for significance (default: 0.05)
#'
#' @return Data frame with gene, topic, postmean, lfsr, z columns
get_top_de_genes <- function(dfa, number_of_genes = -1, lfsr_threshold = 0.05) {
  postmean <- dfa$postmean %>% melt(value.name = 'postmean')
  lfsr <- dfa$lfsr %>% melt(value.name = 'lfsr')
  z <- dfa$z %>% melt(value.name = 'z')
  
  de_genes <- postmean %>% 
    dplyr::left_join(lfsr, by = c("Var1", "Var2")) %>%
    dplyr::left_join(z, by = c("Var1", "Var2")) %>%
    rename_with(~sub("Var1", "gene" , .x)) %>% 
    rename_with(~sub("Var2", "topic" , .x)) %>% 
    dplyr::filter(lfsr < lfsr_threshold) %>%
    dplyr::arrange(topic, desc(postmean)) %>% 
    dplyr::mutate(gene=as.character(gene))
  
  if (number_of_genes > 0) {
    de_genes <- de_genes %>% 
      group_by(topic) %>%
      top_n(n = number_of_genes, wt = postmean) %>% 
      ungroup()
  }
  
  de_genes
}

# correlations ------------------------------------------------------------

#' Internal: Compute correlation between two matrices
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix
#' @param reorder Hierarchically cluster result (default: TRUE)
#' @param method Correlation method (default: 'pearson')
#' @param ignore_missing_rows Use only common rows (default: TRUE)
#'
#' @return Correlation matrix between mat1 and mat2 columns
.correlation <- function(mat1, mat2, reorder = TRUE, method = 'pearson', ignore_missing_rows = TRUE) {
  # We we want to get a distance matrix in order to run hclust on it.
  # So we compare the correlation between all the topics even from the same fit
  colnames(mat1) <- paste("1", colnames(mat1), sep = "_")
  colnames(mat2) <- paste("2", colnames(mat2), sep = "_")

  if (ignore_missing_rows && nrow(mat1) != nrow(mat2)) {
    common_rownames <- intersect(rownames(mat1), rownames(mat2))
    message("Correlation between ", length(common_rownames), " rows")
    mat1 <- mat1[common_rownames,]
    mat2 <- mat2[common_rownames,]
  }

  mat <- cbind(mat1, mat2)
  correlation <- cor(mat, mat, method = method)
  
  if (isTRUE(reorder) && !is_square_matrix(correlation)) {
    warning("Reorder isn't supported on non square correlation matrix")
    reorder <- FALSE
  }

  if (isTRUE(reorder)) {
    hclustering <- hclust(as.dist((1-correlation) / 2) , method='ward.D2')
    correlation <- correlation[hclustering$order, hclustering$order]
  }
  
  # remove unwanted correlations between topics from the same fit
  correlation <- correlation[!colnames(correlation)  %in% colnames(mat1),
                             !colnames(correlation)  %in% colnames(mat2)]

  colnames(correlation) <- gsub("^1_", "", colnames(correlation))
  rownames(correlation) <- gsub("^2_", "", rownames(correlation))
  correlation
}

#' Compute correlation between topic model fits
#'
#' @param fit.1 First topic model fit
#' @param fit.2 Second topic model fit (default: same as fit.1)
#' @param fit_attribute Which matrix to correlate: 'L' or 'F' (default: 'L')
#' @param reorder Hierarchically cluster result (default: TRUE)
#' @param method Correlation method (default: 'pearson')
#'
#' @return Correlation matrix between topics
topics_correlation <- function(fit.1, fit.2 = NA, fit_attribute = 'L', reorder = TRUE, method = 'pearson') {
  if (missing(fit.2)) {
    fit.2 <- fit.1
  }

  .correlation(mat1 = fit.1[[fit_attribute]],
               mat2 = fit.2[[fit_attribute]],
               reorder = reorder,
               method = method)
}


# Gene ontology -----------------------------------------------------------

#' Perform gene ontology enrichment for topic DE genes
#'
#' @param dfa DE analysis object or data frame with gene/topic columns
#' @param lfsr_threshold Significance threshold (default: 0.01)
#' @param universe Background gene universe (default: NULL)
#' @param ... Additional arguments passed to EnrichmentAnalysis
#'
#' @return Data frame of GO enrichment results
topic_gene_ontology <- function(dfa, lfsr_threshold = 0.01, universe = NULL, ...) {
  if (is(dfa, 'topic_model_de_analysis')) {
    postmean <- dfa$postmean %>%
      reshape2::melt(value.name = 'postmean', variable.factor = FALSE) %>%
      rename_with(~sub("Var1", "gene" , .x)) %>%
      rename_with(~sub("Var2", "topic" , .x))

    lfsr <- dfa$lfsr %>%
      reshape2::melt(value.name = 'lfsr', variable.factor = FALSE) %>%
      rename_with(~sub("Var1", "gene" , .x)) %>%
      rename_with(~sub("Var2", "topic" , .x))

    de_genes <- postmean %>%
      dplyr::left_join(lfsr, by = c("gene", "topic"))

    if (is.null(universe)) {
      universe <- GeneIdMapping()$ids[rownames(dfa$z)]
    }
  } else {
    de_genes <- dfa
    stopifnot("Missing universe parameters"=!is.null(universe))
  }

  if ('postmean' %in% colnames(de_genes)) {
    de_genes <- de_genes %>% filter(postmean > 0 & lfsr < lfsr_threshold)
  }

  de_genes <- de_genes %>%
              dplyr::mutate(gene=as.character(gene)) %>%
              dplyr::mutate(id=GeneIdMapping()$ids[gene])

  formula <- as.formula("id~topic")
  de_go <- EnrichmentAnalysis(de_genes, formula = formula, universe=universe, ...)
  
  for(db_name in names(de_go)) {
    de_go[[db_name]]@compareClusterResult['db'] <- db_name
  }
  
  go_results <- do.call("rbind", lapply(de_go, attr, 'compareClusterResult')) %>%
    dplyr::arrange(topic, Description) %>% 
    dplyr::select(-c(desc_id, Cluster))
  
  gene_id_to_names <- GeneIdMapping()$names
  go_result_names <- lapply(str_split(go_results$geneID, "/"), function(x) unname(gene_id_to_names[x]))
  go_results$geneName <- go_result_names %>% sapply(str_c, collapse="/")
  go_results
}

# Per person quantification  ---------------------------------------------

#' Compute mean topic scores per person
#'
#' @param obj Seurat object
#' @param topics_column_names Character vector of topic column names
#' @param person_id_column Column with person IDs (default: 'projid')
#' @param score_column_prefix Prefix for output columns (default: '')
#'
#' @return Data frame with mean topic scores per person
compute_per_person_topic_score <- function(obj,
                                           topics_column_names, 
                                           person_id_column = 'projid',
                                           score_column_prefix = '') {
  person_id_and_topic_scores <- FetchData(obj, c(person_id_column, topics_column_names))
  mean_topic_per_person <- person_id_and_topic_scores %>%
    dplyr::group_by(across(all_of(person_id_column))) %>% 
    summarise(across(all_of(topics_column_names), list(mean = mean))) %>% 
    rename_with(~gsub("^X|_mean", "", .x)) %>% 
    rename_with(~paste0(score_column_prefix, .x), .cols = -all_of(person_id_column))
  
  mean_topic_per_person
}

#' Compute mean topic scores per person for specific annotation
#'
#' @param obj Seurat object
#' @param topics_column_names Character vector of topic column names
#' @param person_id_column Column with person IDs (default: 'projid')
#' @param score_column_prefix Prefix for output columns (default: '')
#' @param an Annotation value to filter on
#'
#' @return Data frame with mean topic scores per person for annotation subset
compute_per_person_topic_score_annotation <- function(obj,
                                           topics_column_names,
                                           person_id_column = 'projid',
                                           score_column_prefix = '', an) {
  person_id_and_topic_scores <- FetchData(obj, c(person_id_column, topics_column_names, 'annotation'))
  mean_topic_per_person <- person_id_and_topic_scores %>% filter(annotation == an) %>%
    dplyr::group_by(across(all_of(person_id_column))) %>%
    summarise(across(all_of(topics_column_names), list(mean = mean))) %>%
    rename_with(~gsub("^X|_mean", "", .x)) %>%
    rename_with(~paste0(score_column_prefix, .x), .cols = -all_of(person_id_column))

  mean_topic_per_person
}

#' Compute cell type proportions per person
#'
#' @param obj Seurat object with cell.type and projid columns
#' @param cell.types Character vector of cell types to include
#'
#' @return Data frame with cell type proportions per person
compute_cell_types_proportion <- function(obj, cell.types) {
  FetchData(obj, c("cell.type", "projid")) %>%
    dplyr::group_by(cell.type, projid) %>%
    count %>%
    filter(cell.type %in% cell.types) %>% 
    reshape2::dcast(., projid~cell.type, value.var = "n", fill = 0) %>% 
    drop_na(projid) %>% 
    column_to_rownames("projid") %>%
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    rownames_to_column(var = 'projid') %>% 
    mutate(projid = as.numeric(projid))
}


##########
# Dynamics
##########

#' Fit Green's state dynamics over program's pseudotime
#'
#' @description Fits GAM models for topic scores along pseudotime trajectory.
#'
#' @param data AnnData-like object with topic scores
#' @param data.da AnnData object with trajectory information
#' @param bootstrap.iterations Number of bootstrap iterations (default: 100)
#' @param bootstrap.proportion Proportion for bootstrap (default: 0.8)
#' @param bootstrap Enable bootstrapping (default: FALSE)
#'
#' @return Modified data object with fitted dynamics
fit_state_dynamics <- function(data, data.da, bootstrap.iterations=100, bootstrap.proportion=.8, bootstrap=FALSE) {
  trajectories <- data.da$uns$trajectories
  if ("palantir" %in% names(trajectories)) {
    trajectories <- trajectories$palantir
  }

  indecies <- trajectories$branch.probs %>% py_to_r %>% rownames %>% `%in%`(., rownames(data))
  pseudotime <- trajectories$pseudotime[indecies]
  branch.probs <- trajectories$branch.probs[indecies]

  data$uns$trajectories$da_fits_sqrt <- fit.dynamics(
    pseudotime = pseudotime,
    features = data$layers[['sqrt.prev']],
    trajectory.probs = branch.probs,
    trajectory.terminal.pseudotime = setNames(py_to_r(trajectories$terminals)[,1], trajectories$terminals$index),
    bootstrap=bootstrap,
    bootstrap.iterations=bootstrap.iterations,
    bootstrap.proportion=bootstrap.proportion)

  data
}
