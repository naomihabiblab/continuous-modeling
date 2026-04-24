#' @title Pathway Dynamics Analysis
#' @description Functions for computing pathway signature dynamics along disease trajectories.

source("R/utils.R")
source("R/pseudobulk.funs.R")
source("R/landscape.analysis.funs.R")


pseudobulks <- load.pseudobulks.of.projid()

#' Fit pathway signature dynamics along trajectories
#'
#' @description Computes pathway gene signatures from pseudobulk data and fits their
#'              dynamics along disease trajectories using GAM models.
#'
#' @param data AnnData object with trajectory information in data$uns$trajectories$palantir
#' @param genes Named list of gene vectors per cell type, or single vector for all types
#' @param pseudobulks Named list of pseudobulk expression matrices per cell type
#' @param scale Scale expression values (default: TRUE)
#' @param signature If TRUE, compute mean signature score; if FALSE, return individual genes (default: TRUE)
#' @param cell.types Character vector of cell types to include (default: all in pseudobulks)
#'
#' @return Dynamics fit object from fit.dynamics() with pathway scores over pseudotime
#'
#' @details For each cell type, extracts specified genes from pseudobulk data, optionally
#'          scales and aggregates to signature scores, then fits dynamics using fit.dynamics().
pathways.dynamics <- function (data, genes, pseudobulks, scale=T, signature=T, cell.types=names(pseudobulks)) {
  features <- sapply(cell.types, function (ct) {
    
    if (ct %in% names(genes)) {
      ct.genes <- genes[[ct]]
    } else {
      ct.genes <- genes
    }
    
    if (ct == 'cux2p') {
      ct <- 'cux2+'
    } else if (ct == 'cux2m') {
      ct <- 'cux2-'
    }
    
    features <- pseudobulks[[ct]] %>% dplyr::select(ct.genes) %>%
      filter(rownames(.) %in% rownames(data)) %>%
      arrange(match(rownames(.), rownames(data)))
    
    if (scale) {
      features <- features %>% scale
    }
    
    if (signature) {
      features <- data.frame(pathways.score=features %>% rowMeans())  
    }
    
    features
  }, simplify = FALSE) %>% as.data.frame()

  if (signature) {
    colnames(features) <- cell.types
  }

  dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                           features = features,
                           trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                           trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                           evaluate.fit = T,
                           bootstrap = F)
  
  dynamics
}
