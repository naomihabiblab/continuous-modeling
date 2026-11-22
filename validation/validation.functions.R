#' @title Validation Functions
#' @description Functions for validating program signatures across datasets,
#'              computing gene set coherence, and managing HDF5 data storage.
#' @notes go together with python implemntation in validation/validation_functions.py

library(rhdf5)
library(dplyr)

source("R/utils.R")
source("R/data_500_utils.R")
source("R/pathways.analysis.R")

#' Convert program name to standardized format
#'
#' @description Normalizes program names by removing dots, converting to lowercase,
#'              and removing cell type prefixes.
#'
#' @param name Character string of program name to convert
#'
#' @return Character string with standardized program name
convert_topic_name <- function(name) {
  # Remove dots, convert to lower, replace ".x." or ".X." with "_"
  name <- gsub("\\.", "_", name)
  name <- tolower(name)
  # If pattern like cux2m_x_k1, convert to cux2m_k1
  name <- gsub("_\\d+_", "_", name)
  name <- gsub("exc_", "", name)
  name <- gsub("inh_", "", name)
  name
}

#' Get genes associated with each topic for a cell type
#'
#' @description Loads DE genes from HDF5, groups by topic, and creates
#'              gene lists including union genes for neuronal annotations.
#'
#' @param cell.type Character string: "inhibitory", "excitatory", "microglia", etc.
#'
#' @return Named list with three elements:
#'         - genes_by_topic: All DE genes per topic
#'         - genes_by_topic_top_50: Top 50 genes by z_score_log per topic
#'         - unique_genes_by_topic: Genes appearing in only one topic
get_genes_by_topic <- function(cell.type) {
    if (cell.type == "excitatory") {
        des_cux2p <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["cux2p"], "/de"))
        des_cux2m <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["cux2m"], "/de"))

        des <- rbind(des_cux2p, des_cux2m)
        des <- des %>% mutate(topic = paste0(.data[['cell.type']], "_", topic))
    } else {
        des <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/de"))
    }

    genes_by_topic <- des %>%
        group_by(topic) %>%
        summarise(genes = list(gene)) %>%
        { setNames(.$genes, .$topic) }


    # add union of genes for topics in neuronal_topics_to_annotation
    if (cell.type %in% c("inhibitory", "excitatory")) {
        for (union_name in names(neuronal_topics_to_annotation)) {
            topic_names_raw <- neuronal_topics_to_annotation[[union_name]]
            topic_names <- vapply(names(topic_names_raw), convert_topic_name, character(1))
            # Only keep topic names that exist in genes_by_topic
            topic_names <- topic_names[topic_names %in% names(genes_by_topic)]
            if (length(topic_names) > 0) {
                union_genes <- unique(unlist(genes_by_topic[topic_names]))
                genes_by_topic[[union_name]] <- union_genes
            }
        }

        genes_by_topic <- genes_by_topic[names(neuronal_topics_to_annotation)]
    }

    genes_by_topic_top_50 <- des %>%
        group_by(topic) %>% 
        top_n(50, wt = z_score_log) %>% 
        summarise(genes = list(gene)) %>%
        { setNames(.$genes, .$topic)}

    topic_gene_df <- des %>% dplyr::select(topic, gene)

    # Count in how many topics each gene appears
    # Use genes_by_topic to compute unique genes and then filter non-unique genes from genes_by_topic

    # Flatten genes_by_topic to get a table of gene, topic
    gene_topic_list <- lapply(names(genes_by_topic), function(topic) {
        data.frame(topic = topic, gene = genes_by_topic[[topic]], stringsAsFactors = FALSE)
    })
    gene_topic_df <- do.call(rbind, gene_topic_list)

    # Count how many topics each gene appears in
    unique_genes <- gene_topic_df %>%
        group_by(gene) %>%
        summarise(n_topics = n_distinct(topic), .groups = "drop") %>%
        filter(n_topics == 1) %>%
        pull(gene)

    # For each topic, keep only genes that are unique
    unique_genes_by_topic <- lapply(genes_by_topic, function(glist) {
        intersect(glist, unique_genes)
    })

    return(list(genes_by_topic=genes_by_topic, genes_by_topic_top_50=genes_by_topic_top_50, unique_genes_by_topic=unique_genes_by_topic))
}

#' Get pathway gene signatures for a cell type
#'
#' @description Loads pathway-to-gene mappings, adding neuronal-specific
#'              pathways for inhibitory/excitatory cell types.
#'
#' @param cell.type Character string specifying the cell type
#'
#' @return Named list mapping pathway names to gene vectors
get_pathway_signatures <- function(cell.type) {
    # TODO: move this file
    pathways_to_genes_path <- file.path(get_result_dir(), "topic_space", glue::glue("Pathways_genes_{cell.type}.json"))
    
    if (file.exists(pathways_to_genes_path)) {
        pathway_to_genes <- jsonlite::read_json(pathways_to_genes_path, simplifyVector=TRUE)
    } else {
        pathway_to_genes <- list()
    }
    
    
    if(cell.type %in% c("inhibitory", "excitatory", "cux2p", "cux2m")) {
    
        pathways_to_name_and_topic_group <- list(
        "hsa04721" = c("Synaptic vesicle cycle", "increasing_early"),
        "R-HSA-5661231"=c("Metallothioneins bind metals", "increasing_early"),
        "GO:0051402" = c("neuron apoptotic process", "increasing_early"),
        "GO:0006979" = c("response to oxidative stress", "increasing_early"),
        "R-HSA-68886" = c("M Phase", "increasing_late"),
        "GO:0071103" = c("DNA conformation change", "increasing_late"),
        "R-HSA-73894" = c("DNA Repair", "increasing_late"),
        "R-HSA-3371453" = c("Regulation of HSF1-mediated heat shock response", "increasing_late"),
        # Decreasning pathways
        "GO:0005761" = c("mitochondrial ribosome", "decreasing_early"),
        "GO:0051028" = c("mRNA transport", "decreasing_late"),
        "GO:0099572" = c("postsynaptic specialization", "decreasing_late"),
        "antinsense" = c("antinsense", "decreasing_early")
    )


        neuornal_pathway_to_genes <- sapply(names(pathways_to_name_and_topic_group), function(pathway_id) {
            group_name <- pathways_to_name_and_topic_group[[pathway_id]][[2]]
            pathway_name <- pathways_to_name_and_topic_group[[pathway_id]][[1]]
            message(pathway_name)
            group <- neuronal_topics_to_annotation[[group_name]]
            de <- read.dataset.for.topics(names(group), dataset_path = "/de") 
            if (cell.type == "inhibitory") {
                de <- de %>% filter(.data[['cell.type']] == "inhibitory")
            } else if (cell.type == "excitatory") {
                de <- de %>% filter(.data[['cell.type']] %in% c("cux2p", "cux2m"))
            } else if (cell.type %in% c("cux2p", "cux2m")) {
                de <- de %>% filter(.data[['cell.type']] == .env[['cell.type']])
            }

            if (pathway_id == "antinsense") {
                de %>% filter(grepl("-AS\\d+", gene)) %>% pull(gene) %>% unique
            } else {
                pathway.genes <- get.pathway.genes(pathway_id)
                de %>% filter(gene %in% pathway.genes) %>% pull(gene) %>% unique
            }
        }, simplify = FALSE, USE.NAMES = TRUE)


        names(neuornal_pathway_to_genes) <- sapply(names(neuornal_pathway_to_genes), function(pathway_id) {
            pathways_to_name_and_topic_group[[pathway_id]][[1]]
        })

        pathway_to_genes <- c(pathway_to_genes, neuornal_pathway_to_genes)
    }


    return(pathway_to_genes)
}

#' Save correlation results to HDF5
#'
#' @description Writes correlation matrix and associated statistics to HDF5 file.
#'
#' @param correlation_object List with cor, pval, adj.pval, sig, names, and params
#' @param cell.type Optional cell type for path construction (default: NULL)
#' @param group_name HDF5 group name for storage
#' @param h5_path Path to HDF5 file (default: TOPICS_DATA_H5_PATH)
#'
#' @return NULL (called for side effects)
h5_save_correlation <- function(correlation_object, cell.type=NULL, group_name, h5_path=TOPICS_DATA_H5_PATH) {
    if (!is.null(cell.type)) {
        group_path <- paste0(CELL_TYPE_TO_GROUP[cell.type], "/", group_name)
    } else {
        group_path <- group_name
    }

    h5createGroup(h5_path, group_path)
    h5write(correlation_object$cor, h5_path, paste0(group_path, "/corr"))
    h5write(correlation_object$pval, h5_path, paste0(group_path, "/pval"))
    h5write(correlation_object$adj.pval, h5_path, paste0(group_path, "/adj.pval"))
    h5write(correlation_object$sig, h5_path, paste0(group_path, "/sig"))
    h5write(correlation_object$names, h5_path, paste0(group_path, "/names"))

    h5createGroup(h5_path, paste0(group_path, "/params"))
    h5write(correlation_object$params$cor.method, h5_path, paste0(group_path, "/params/cor.method"))
    h5write(correlation_object$params$cor.use, h5_path, paste0(group_path, "/params/cor.use"))
    h5write(correlation_object$params$p.adjust.method, h5_path, paste0(group_path, "/params/p.adjust.method"))
}

#' Read correlation results from HDF5
#'
#' @description Loads correlation matrix and statistics from HDF5 file.
#'
#' @param cell.type Optional cell type for path construction (default: NULL)
#' @param group_name HDF5 group name to read from
#' @param h5_path Path to HDF5 file (default: TOPICS_DATA_H5_PATH)
#'
#' @return List with cor, pval, adj.pval, sig, names, and params
h5_read_correlation <- function(cell.type=NULL, group_name, h5_path=TOPICS_DATA_H5_PATH) {
    if (!is.null(cell.type)) {
        group_path <- paste0(CELL_TYPE_TO_GROUP[cell.type], "/", group_name)
    } else {
        group_path <- group_name
    }

    correlation_object <- list(
        cor = h5read(h5_path, paste0(group_path, "/corr")),
        pval = h5read(h5_path, paste0(group_path, "/pval")),
        adj.pval = h5read(h5_path, paste0(group_path, "/adj.pval")),
        sig = h5read(h5_path, paste0(group_path, "/sig")),
        names = h5read(h5_path, paste0(group_path, "/names"))
    )
    correlation_object$params <- list(
        cor.method = h5read(h5_path, paste0(group_path, "/params/cor.method")),
        cor.use = h5read(h5_path, paste0(group_path, "/params/cor.use")),
        p.adjust.method = h5read(h5_path, paste0(group_path, "/params/p.adjust.method"))
    )
    return(correlation_object)
}

#' Filter correlation object by index
#'
#' @description
#' Takes a square correlation object (as loaded by \code{h5_read_correlation}) and subsets it by a set of row/column indices.
#'
#' @param c A correlation object list with elements 'cor', 'pval', 'adj.pval', 'sig', 'names', and 'params'.
#' @param index Integer or logical vector of indices to subset both dimensions.
#'
#' @return A correlation object subset to the specified indices.
#'
#' @examples
#' \dontrun{
#'   c <- h5_read_correlation(...)
#'   idx <- 1:10
#'   c_filt <- filter_correlation_by_index(c, idx)
#' }
filter_correlation_by_index <- function(c, index) {
  c_subset <- list(
    cor = c$cor[index, index],
    names = list(c$names[[1]][index], c$names[[2]][index]),
    pval = c$pval[index, index],
    adj.pval = c$adj.pval[index, index],
    sig = c$sig[index, index],
    params = c$params
  )

  return(c_subset)
}


#' Save mean signature scores to HDF5
#'
#' @description Aggregates signature scores by donor and optionally by subtype,
#'              then saves to HDF5 file.
#'
#' @param df Data frame with Signature_, SignatureUnique_, Pathway_ columns
#' @param cell.type Cell type for HDF5 path construction
#' @param validation_name Name of the validation dataset
#' @param group.by Column to group by for aggregation (default: 'projid')
#' @param h5_path Path to HDF5 file (default: TOPICS_DATA_H5_PATH)
#' @param subtype_column Optional column for subtype grouping (default: NULL)
#'
#' @return NULL (called for side effects)
h5_save_mean_scores <- function(df, cell.type, validation_name, group.by='projid', h5_path=TOPICS_DATA_H5_PATH, subtype_column=NULL) {
    group_path <- paste0(CELL_TYPE_TO_GROUP[cell.type], "/", validation_name)
    h5createGroup(h5_path, group_path)

    signatures_df <- df %>% dplyr::select(starts_with("Signature_"), starts_with("SignatureUnique_"), starts_with("Pathway_"), group.by, subtype_column)

    mean_df <- signatures_df %>%
        group_by(across(all_of(c(group.by)))) %>%
        summarise(across(starts_with("Signature_"), mean), across(starts_with("SignatureUnique_"), mean), across(starts_with("Pathway_"), mean))

    h5write(mean_df, h5_path, paste0(group_path, "/mean_scores"))

    if (!is.null(subtype_column)) {
        mean_by_subtype_df <- signatures_df %>%
            group_by(across(all_of(c(group.by, subtype_column)))) %>%
            summarise(across(starts_with("Signature_"), mean), across(starts_with("SignatureUnique_"), mean), across(starts_with("Pathway_"), mean)) %>%
            rename_with(~gsub("_mean", "", .x)) %>% ungroup()

        h5write(mean_by_subtype_df, h5_path, paste0(group_path, "/mean_scores_by_subtype"))
    }
}

#' Check if HDF5 path exists
#'
#' @param file Path to HDF5 file
#' @param path HDF5 internal path to check
#'
#' @return Logical indicating whether path exists
h5_exists <- function(file, path) {
  info <- h5ls(file, all = TRUE)
  any(paste0(info$group, "/", info$name) == path |
      (info$group == path & info$name == ""))
}

#' Compute cross-cell-type program correlations
#'
#' @description Loads mean signature scores for all cell types, merges by donor,
#'              computes correlation matrix, and saves to HDF5.
#'
#' @param project Project/validation dataset name
#' @param identity_column Column name for donor identity (default: 'projid')
#'
#' @return NULL (called for side effects - saves to HDF5)
compute_program_environment_correlation <- function(project, identity_column='projid') {
    cell_types <- c("inhibitory", "excitatory", "microglia", "astrocytes", "oligo", "opcs")

    df_list <- list()
    for (cell.type in cell_types) {
        group_path <- paste0("/", CELL_TYPE_TO_GROUP[cell.type], "/", project, "/mean_scores")
        dataset_name <- "/"
        if (!h5_exists(TOPICS_DATA_H5_PATH, group_path)) {
            warning(paste("Dataset missing for", group_path, "- skipping"))
            next  
        }
        df <- h5read(TOPICS_DATA_H5_PATH, group_path)

        # add cell type prefix to column names
        df <- df %>%
            filter(!is.na(get(identity_column))) %>%
            dplyr::select(starts_with("SignatureUnique_"), all_of(identity_column)) %>% 
            rename_with(~paste0(cell.type, "_", .x), .cols = -which(names(.) == identity_column))
    
        df_list[[cell.type]] <- df
    }

    if(length(df_list) == 0) {
        merged_df <- NULL
    } else {
        merged_df <- Reduce(function(x, y) {
            dplyr::full_join(x, y, by = identity_column)
        }, df_list)

        merged_df <- merged_df %>% 
            column_to_rownames(identity_column)
    }
    
    # count number of rows with missing values
    missing_rows <- merged_df %>% 
        filter_all(any_vars(is.na(.))) %>%
        nrow()

    if (missing_rows > 0) {
        warning(paste("Missing values in", missing_rows, "rows"))
        merged_df <- merged_df %>%
            filter_all(any_vars(!is.na(.)))
    }

    colnames(merged_df) <- sapply(colnames(merged_df), gsub, pattern = "SignatureUnique_", replacement = "")

   # filter out not in community
    cor_mat <- correlate_dataframes(merged_df, use = 'complete.obs', method = 'pearson')

    group_name <- "/validation"
    h5createGroup(TOPICS_DATA_H5_PATH, group_name)
    group_name <- paste0(group_name, "/", project)
    h5createGroup(TOPICS_DATA_H5_PATH, group_name)
    h5_save_correlation(cor_mat, group_name=paste0(group_name, "/programs_correlation"))
}


# Test coherence of signatures

#' Compute mean pairwise correlation for a gene set and compare to random background
#'
#' @param expr_matrix Expression matrix (genes x cells), e.g., from GetAssayData()
#' @param genes Character vector of gene names to test
#' @param n_permutations Number of random background permutations (default: 100)
#' @param method Correlation method: "spearman" or "pearson" (default: "pearson")
#' @param seed Random seed for reproducibility (default: 1)
#' @param show_progress Show progress bar (default: TRUE, ignored when n_cores > 1)
#' @param n_cores Number of cores for parallel processing (default: future::availableCores())
#'               Set to 1 for sequential processing with progress bar
#' @return List with observed mean correlation, background distribution, p-value, and z-score
compute_gene_set_coherence <- function(expr_matrix, genes, n_permutations = 100, 
                                        method = "pearson", seed = 1,
                                        show_progress = TRUE,
                                        n_cores = future::availableCores()) {
  if (!is.null(seed)) set.seed(seed)
  
  # Filter genes to those present in expression matrix
  valid_genes <- genes[genes %in% rownames(expr_matrix)]
  n_genes <- length(valid_genes)
  
  if (n_genes < 2) {
    warning("Less than 2 valid genes found. Cannot compute correlation.")
    return(list(observed_mean = NA, bg_means = NA, p_value = NA, z_score = NA, n_genes = n_genes))
  }
  
  all_genes <- rownames(expr_matrix)
  
  # Helper function to compute mean of upper triangle correlation
  # Memory-efficient: computes cor matrix and extracts mean in one step
  compute_mean_cor <- function(gene_subset) {
    expr_sub <- expr_matrix[gene_subset, , drop = FALSE]
    # Use cor() with use="pairwise.complete.obs" for robustness
    cor_mat <- cor(t(as.matrix(expr_sub)), method = method, use = "pairwise.complete.obs")
    # Extract upper triangle (excluding diagonal) and compute mean
    mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  }
  
  # Compute observed mean correlation for the gene set and print the time it took
  tictoc::tic(msg = glue::glue("Computing observed mean correlation for gene set of size {length(valid_genes)} using #{n_cores} cores..."), log = TRUE)
  observed_mean <- compute_mean_cor(valid_genes)
  tictoc::toc()
  
  # Compute background distribution using random gene sets of the same size
  message(glue::glue("Computing background distribution for gene set of size {length(valid_genes)} ({n_permutations} permutations, using #{n_cores} cores)..."))
  
  # Use parallel processing if n_cores > 1
  if (n_cores > 1) {
    # Pre-generate all random gene sets for reproducibility
    if (!is.null(seed)) set.seed(seed + 1)  # Different seed for background
    random_gene_sets <- lapply(seq_len(n_permutations), function(i) {
      sample(all_genes, n_genes, replace = FALSE)
    })
    
    # Parallel computation
    bg_means <- unlist(parallel::mclapply(random_gene_sets, compute_mean_cor, mc.cores = n_cores))
    
  } else {
    # Sequential computation with progress bar
    bg_means <- numeric(n_permutations)
    
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = n_permutations, style = 3, width = 50)
    }
    
    for (i in seq_len(n_permutations)) {
      random_genes <- sample(all_genes, n_genes, replace = FALSE)
      bg_means[i] <- compute_mean_cor(random_genes)
      
      if (show_progress) setTxtProgressBar(pb, i)
    }
    
    if (show_progress) close(pb)
  }
  
  # Compute empirical p-value (one-tailed: is observed higher than random?)
  # Add 1 to numerator and denominator for conservative estimate (Phipson & Smyth, 2010)
  p_value <- (sum(bg_means >= observed_mean) + 1) / (n_permutations + 1)
  
  # Compute effect size metrics
  bg_mean <- mean(bg_means, na.rm = TRUE)
  bg_sd <- sd(bg_means, na.rm = TRUE)
  
  # Z-score (standardized difference)
  z_score <- if (bg_sd > 0) (observed_mean - bg_mean) / bg_sd else NA
  
  # Cohen's d (same as z-score in this context, but conventionally reported)
  cohens_d <- z_score
  
  # Percentile rank: proportion of background values below observed
  # Higher = observed is more extreme (better coherence)
  percentile_rank <- mean(bg_means < observed_mean, na.rm = TRUE) * 100
  
  # Fold enrichment: ratio of observed to background mean
  # >1 means gene set is more coherent than random
  fold_enrichment <- if (bg_mean != 0) observed_mean / bg_mean else NA
  
  # Cliff's delta: non-parametric effect size
  # Ranges from -1 to 1; >0 means observed tends to be higher than background
  # Formula: (# times obs > bg - # times obs < bg) / total comparisons
  cliffs_delta <- (sum(observed_mean > bg_means) - sum(observed_mean < bg_means)) / n_permutations
  
  list(
    observed_mean = observed_mean,
    bg_means = bg_means,
    bg_mean = bg_mean,
    bg_sd = bg_sd,
    p_value = p_value,
    z_score = z_score,
    cohens_d = cohens_d,
    percentile_rank = percentile_rank,
    fold_enrichment = fold_enrichment,
    cliffs_delta = cliffs_delta,
    n_genes = n_genes,
    n_permutations = n_permutations
  )
}

#' Compute coherence for multiple gene sets
#'
#' @param expr_matrix Expression matrix (genes x cells)
#' @param gene_sets Named list of character vectors (gene names per set)
#' @param n_permutations Number of random background permutations (default: 50)
#' @param method Correlation method (default: "pearson")
#' @param seed Random seed (default: 1)
#' @param show_progress Show progress bar for each gene set (ignored when n_cores > 1)
#' @param n_cores Number of cores for parallel processing (default: future::availableCores())
#' @return Data frame with coherence statistics for each gene set
compute_coherence_multiple <- function(expr_matrix, gene_sets, n_permutations = 50,
                                       method = "pearson", seed = 1,
                                       show_progress = TRUE,
                                       n_cores = future::availableCores()) {
  set.seed(seed)
  
  message(glue::glue("Using {n_cores} cores for parallel processing"))
  
  results <- lapply(names(gene_sets), function(name) {
    message(glue::glue("\nComputing coherence for {name} of size {length(gene_sets[[name]])}..."))
    res <- compute_gene_set_coherence(expr_matrix, gene_sets[[name]], 
                                       n_permutations = n_permutations, 
                                       method = method, seed = NULL,
                                       show_progress = show_progress,
                                       n_cores = n_cores)
    data.frame(
      gene_set = name,
      n_genes = res$n_genes,
      observed_mean_cor = res$observed_mean,
      bg_means_cor = I(list(res$bg_means)),
      bg_mean_cor = res$bg_mean,
      bg_sd_cor = res$bg_sd,
      z_score = res$z_score,
      cohens_d = res$cohens_d,
      percentile_rank = res$percentile_rank,
      fold_enrichment = res$fold_enrichment,
      cliffs_delta = res$cliffs_delta,
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
  })
  
  result_df <- bind_rows(results)
  # Adjust p-values for multiple testing

  result_df$p_adj <- p.adjust(result_df$p_value, method = "BH")
  result_df
}
