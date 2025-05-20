source("utils.R")
source("pseudobulk.funs.R")
source("landscape.analysis.funs.R")





pseudobulks <- load.pseudobulks.of.projid()


is.nested <- function(x) {
  if (class(x) != "list") {
    stop("Expecting 'x' to be a list")
  }
  out <- any(sapply(x, is.list))
  return(out)
}

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

unify_dynamics <- function(fits) {
  has_duplicated_features <- lapply(fits, function(f) f$fitted.vals %>% pull(feature) %>% unique) %>% unlist %>% duplicated %>% any
  stopifnot(!has_duplicated_features)

  fit <- list(fitted.vals=NULL, pred.vals=NULL, evaluations=NULL, number_of_points=fits[[1]]$number_of_points)
  fit$fitted.vals <- lapply(fits, function(f) f$fitted.vals) %>% do.call(rbind, .)
  fit$pred.vals <- lapply(fits, function(f) f$pred.vals) %>% do.call(rbind, .)
  fit$evaluations <- lapply(fits, function(f) f$evaluations) %>% do.call(rbind, .)

  fit
}