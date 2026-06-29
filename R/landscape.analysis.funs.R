#' @title Landscape Analysis Functions
#' @description Functions for trajectory analysis, dynamics fitting, trait associations,
#'              and community detection in single-cell data using Palantir and Leiden clustering.

#####################################################################################################################
#                                           General settings and functions                                          #
#####################################################################################################################

# if(!"scanpy.env" %in% reticulate::conda_list()[,"name"]) {
#   # Create conda envoriment for running python analysis (mainly scanpy+palantir) through R
#   reticulate::conda_create("scanpy.env",  packages = c("numpy","seaborn", "scikit-learn", "statsmodels", "numba", "pytables"))
#   reticulate::conda_install("scanpy.env", packages = c("python-igraph","leidenalg"))
#   reticulate::conda_install("scanpy.env", packages = c("scanpy","palantir","phate","fa2"), pip = T)
# }

# reticulate::use_condaenv("scanpy.env")
invisible(lapply(c("dplyr","tibble","reticulate","reshape2","anndata","progress","rhdf5", "tidyr"), library, character.only = T))



#####################################################################################################################
#                                                 Analysis Functions                                                #
#####################################################################################################################

#' Compute local-similarities in embedding
#'
#' Pairwise adaptive RBF distance
#'
#' @param embedding dataframe-like of shape (n_samples, n_features)
#' @param knn integer used to calculate the local density around each point - mean distance to the k nearest neighbors
#'
#' @return matrix of shape (n_samples, n_samples) of pairwise similarities
embedding.similarity <- function(embedding, knn = 5) {
  sim <- rdist::cdist(embedding, embedding)^2
  s <- apply(sim, 1, function(.) sqrt(mean(sort(.)[2:(knn+1)])))
  sim <- exp(-sim / outer(s, s))
  sim[is.na(sim)] <- 0
  sim
}

#' Perform trait associations 
#'
#' Associate co-variates to traits will controlling for technical variables
#' 
#' @param traits dataframe of shape (n_samples, n_traits)
#' @param covariates dataframe of shape (n_samples, n_covariates)
#' @param controls dataframe of shape (n_samples, n_controls)
#'
#' @return dataframe with trait association results for every n_traits*n_covariates 
#' while correcting for multiple hypothesis testing within each of the given traits
associate.traits <- function(traits, covariates, controls, p.adjust.method="BH") {
  df <- data.frame(traits, covariates, controls, check.names = TRUE)
  
  # Handle covariate that include spaces or spacial characters  
  covariates_original_names <- colnames(covariates)
  names(covariates_original_names) <- make.names(covariates_original_names)
  
  
  df <- do.call(rbind, lapply(colnames(traits), function(trait)
    do.call(rbind, lapply(make.names(colnames(covariates)), function(covariate) {
      # Set controlled covariates
      control <- colnames(controls)
      
      if (length(control) > 0) {
        control <- case_when(trait == "cogng_demog_slope" ~ paste(setdiff(control, c("msex", "age_death")), collapse = " + "),
                      T ~ paste(control,collapse = " + "))
        formula <- stringr::str_interp("${trait} ~ `${covariate}` + ${control}")
      } else {     
        formula <- stringr::str_interp("${trait} ~ `${covariate}`")
      }

      m <- summary(lm(formula, df[!is.na(df[,trait]), ]))
      
      data.frame(trait = trait, 
                 covariate = covariates_orginal_names[[covariate]],
                 beta  = m$coefficients[covariate, 1],
                 se = m$coefficients[covariate, 2],
                 tstat = m$coefficients[covariate, 3],
                 pval  = m$coefficients[covariate, 4],
                 r.sq  = m$adj.r.squared,
                 n     =  sum(!is.na(df[,trait])),
                 formula = formula)
    }))))
  return(df %>% dplyr::group_by(trait) %>%
           dplyr::mutate(adj.pval = p.adjust(pval, method=p.adjust.method),
                  sig = cut(adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))) %>% 
           dplyr::ungroup())
}


#' Perform trait associations 
#'
#' Associate co-variates to traits will controlling for technical variables
#' 
#' @param traits dataframe of shape (n_samples, n_traits)
#' @param covariates dataframe of shape (n_samples, n_covariates)
#' @param controls dataframe of shape (n_samples, n_controls)
#' @param mixed_model_intercept_controls dataframe of shape (n_samples, n_mixed_model_intercept_controls)
#'
#' @return dataframe with trait association results for every n_traits*n_covariates 
#' while correcting for multiple hypothesis testing within each of the given traits
#'
#' @note Uses lmerTest::lmer for mixed-effects models with random intercepts.
associate.traits.lmer <- function(traits, covariates, controls, p.adjust.method="BH", 
                             mixed_model_intercept_controls) {
  df <- data.frame(traits, covariates, controls, mixed_model_intercept_controls, check.names = TRUE)
  
  # Handle covariate that include spaces or spacial characters  
  covariates_original_names <- colnames(covariates)
  names(covariates_original_names) <- make.names(covariates_original_names)
  
  
  df <- do.call(rbind, lapply(colnames(traits), function(trait)
    do.call(rbind, lapply(make.names(colnames(covariates)), function(covariate) {
      # Set controlled covariates
      control <- colnames(controls)
      control <- case_when(trait == "cogng_demog_slope" ~ paste(setdiff(control, c("msex", "age_death")), collapse = " + "),
                           T ~ paste(control,collapse = " + "))
      
      formula <- stringr::str_interp("${trait} ~ ${covariate} + ${control}")
      
      mixed_control_cols <- colnames(mixed_model_intercept_controls)
      formula <- paste0(formula, " + ",  paste(glue::glue("(1|{mixed_control_cols})"), collapse = '+'))
      
      m <- summary(lmerTest::lmer(formula, df[!is.na(df[,trait]), ]))
      data.frame(trait = trait, 
                 covariate = covariates_original_names[[covariate]],
                 beta  = m$coefficients[covariate, 1],
                 se = m$coefficients[covariate, 2],
                 tstat = m$coefficients[covariate, 4],
                 pval  = m$coefficients[covariate, 5],
                 n     =  sum(!is.na(df[,trait])),
                 formula = formula)
    }))))

  return(df %>% dplyr::group_by(trait) %>%
           dplyr::mutate(adj.pval = p.adjust(pval, method=p.adjust.method),
                         sig = cut(adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))) %>% 
           dplyr::ungroup())
}

#' Fit trajectories to cellular landscape using Palantir
#'
#' @param data AnnData object of cellular landscape
#' @param dm number of diffusion map components to compute. Passed to palantir$utils$run_diffusion_maps
#' @param dm.k number of neighbors used when constructing diffusion space. Passed to palantir$utils$run_diffusion_maps
#' @param palantir.k number of neighbors used when running palantir. Passed to palantir$run_palantir
#' @param scale logical. Should components be scaled or not. Passed to palantir$run_palantir
#'
#' @return Palantir result in the form of a named list:
#' -  pseudotime - named vector of donors' pseudotime
#' -  branch.prob - donors' trajectory probabilities
#' -  user.root - datapoint used as root for Palantir's run
#' -  params - named list of input parameters
fit.trajectories.palantir <- function(data, dm, dm.k, palantir.k, scale, root.clusters, exclude.clusters = c(), include=NULL,
                             use_early_cell_as_start=F, early_cell = NULL) {
  pa <- reticulate::import("palantir", convert = FALSE)
  pd <- import("pandas")
  
  if (is.null(early_cell)) {
    X <- t(data[data$obs$clusters %in% root.clusters]$X)
    early_cell <- names(which.min(Matrix::colSums((X - Matrix::rowMeans(X))**2)))
  }
  
  if (is.null(include)) {
    include <- rep(T, nrow(data$X))
  }


  ms.data <- pa$utils$run_diffusion_maps(data.frame(data[include & (!data$obs$clusters %in% exclude.clusters)]$X), n_components = as.integer(dm), knn = as.integer(dm.k)) %>%
    pa$utils$determine_multiscale_space(.)
  p.res <- pa$core$run_palantir(ms.data, early_cell, knn = as.integer(palantir.k), max_iterations = as.integer(100),
                                scale_components = scale,
                                n_jobs = as.integer(1),
                                use_early_cell_as_start = use_early_cell_as_start
                                )
  
  markov <- pa$core[["_construct_markov_chain"]](ms.data, as.integer(palantir.k), 
                                            reticulate::py_get_attr(p.res, "pseudotime"),
                                            as.integer(1))$todense() %>% py_to_r()
  
  res <-(merge(data.frame(row.names = data$obs_names),
               data.frame(pseudotime=as.vector((reticulate::py_get_attr(p.res, "pseudotime"))$values), py_to_r(p.res$branch_probs)),
               by="row.names", all.x=T) %>% 
           tibble::column_to_rownames("Row.names"))[data$obs_names,]

  params <- as.list(match.call())
  params$include <- NULL
  # return(res)
  return(list(
    pseudotime = setNames(res$pseudotime, rownames(res)),
    branch.probs = res[,-1],
    terminals = res[gsub("X","", colnames(res)[-1]),] %>% mutate(terminal = rownames(.)) %>% dplyr::select(-contains("X")),
    transitions = markov,
    user.root = early_cell,
    params = params[3:length(params)]
  ))
}

#' Fit dynamics over pseudotime
#'
#' @description
#' This function computes the dynamics of a feature over the pseudotime in different trajectories
#' For each given feature the function regresses the feature on the pseudotime while weighting samples
#' by the trajectory probabilities.
#'
#' Then, the function uses the fitted model to predict the dynamics over a equally spaced set of pseudotime
#' values ending at the pseudotime of the trajectories' terminal state.
#'
#' @param pseudotime numerical vector of length `n_samples`
#' @param features data.frame of shape (`n_samples`, `n_features`)
#' A data frame of features to fit dynamics for.
#' @param trajectory.probs data.frame of shape (`n_samples`, `n_trajectories`)
#' A data frame of trajectory probabilities used to weight samples by
#' @param min.prob.clip numeric
#' Exclude samples from trajectory specific dynamics whose probability in trajectory is below specified value.
#' Pseudotime values used for predicted dynamics in trajectory begin from the minimal pseudotime of samples
#' whose probability in trajectory is greater or equals to specified value.
#' @param trajectories.order Order of trajectories for factor levels (default: c("South", "prAD", "West", "ABA"))
#' @param evaluate.fit Evaluate model fit against null and unweighted models (default: TRUE)
#' @param weight_filter Minimum weight threshold (default: -Inf)
#' @param ... additional arguments used to `.fit.dynamics`
#'
#' @return Named list with:
#'         - fitted.vals: fitted values from model
#'         - pred.vals: predicted values over pseudotime
#'         - evaluations: model comparison results (if evaluate.fit=TRUE)
#'         - number_of_points: number of data points used
#'
#' @note This function assumes corresponding indices across `pseudotime`, `features` and `trajectory.probs`
#'
fit.dynamics <- function(
    pseudotime,
    features,
    trajectory.probs,
    trajectory.terminal.pseudotime = setNames(rep(1, ncol(trajectory.probs)), colnames(trajectory.probs)),
    min.prob.clip=0,
    trajectories.order = c("South", "prAD", "West", "ABA"),
    evaluate.fit = T, 
    weight_filter = -Inf,
    ...) {

  # Equally spaced pseudotime values to predict dynamics for
  xs.pred <- seq(min(pseudotime, na.rm = T),
                 max(pseudotime, na.rm = T),
                 length.out = 50)

  number_of_points <- pseudotime %>% as.vector %>% na.omit() %>% length

  pb <- progress_bar$new(format="Fitting dynamics along trajectories :current/:total [:bar] :percent in :elapsed. ETA :eta",
                         total = ncol(features)*ncol(trajectory.probs), clear=F, width=100, force = T)
  df <- data.frame(ps=pseudotime, trajectory.probs, features, check.names=FALSE)

  trajectories <- colnames(trajectory.probs)
  trajectories_levels <- trajectories[order(match(trajectories, trajectories.order))]

  res <- lapply(trajectories, function(t) {

    # Set trajectory range between pseudotime of first sample with weight in trajectory over `min.prob.clip`
    # and pseudotime specified by trajectory's terminal pseudotime
    beg <- min(df[df[,t] >= min.prob.clip, ]$ps, na.rm=T) %>% ifelse(is.infinite(.), 0, .)
    end <- trajectory.terminal.pseudotime[gsub("^X","",t)] %>% ifelse(is.infinite(.), 1, .)

    # Specify trajectory specific pseudotime values to predict dynamics for
    xs <- unique(c(beg, xs.pred[(beg <= xs.pred) & (xs.pred <= end)], end))
    
    # Fit- and predict dynamics for features along specific trajectory
    res <- lapply(colnames(features), function(p) {
      # TODO: replace df[,p] with way that don't support partial matching 
      
      # subset dataframe for specified trajectory and parameter to regress
      .df <- data.frame(x=df$ps, y=df[,p], w=df[,t], row.names = rownames(df)) %>%
          tidyr::drop_na() %>%
        filter(beg <= x & x <= end)

      res <- list(fit=data.frame(x=NA, fit = NA, se = NA, fit_sd = NA, se.fit_sd = NA, feature=p, trajectory=t, n=nrow(.df))[0,],
                  prd=data.frame(x=xs, fit = NA, se = NA, fit_sd = NA, se.fit_sd = NA, feature=p, trajectory=t))
      tryCatch({
        res <- .fit.dynamics(.df, xs, evaluate.fit=evaluate.fit, ...)
        res <- lapply(res, function(.df) .df %>% mutate(feature=p, trajectory=factor(t, levels=trajectories_levels)))
      }, error = function(e)
        warning("Failed computing dynamics for feature [", p, "] in trajectory [", t, "]. Error: ", e))

      pb$tick()
      return(res)
    })

    # Merge results of fit- and predict from different features
    return(lapply(1:length(res[[1]]), function(i) do.call(rbind, lapply(res, "[[", i))))
  })

  # Merge results of fit- and predict from different trajectories
  res <- lapply(1:length(res[[1]]), function(i) do.call(rbind, lapply(res, "[[", i)))

  names <- c("fitted.vals", "pred.vals")
  if(evaluate.fit)
    names <- c(names, "evaluations")
  res <- res %>% `names<-`(names)

  res$number_of_points <- number_of_points

  res
}

#' Fit dynamics over specific trajectory and parameter of interest
#'
#' @param df data.frame with columns `x`,`y`,`w` for the covariate, response and sample weight respectively
#' @param pred.x numeric vector of values to predict model outcome for
#' @param bootstrap If TRUE perform bootstrap sampling of the data (passed to the `sample` function as the `replace` value)
#' @param bootstrap.proportion numeric in range [0,1] indicating sub-sampling proportion
#' @param bootstrap.iterations number of bootstrap iterations to perform
#' @param bootstrap.group Aggregate bootstrap results by grouping (default: TRUE)
#' @param evaluate.fit Evaluate model fit (default: TRUE)
#' @param ... Additional arguments passed to mgcv::gam
#'
#' @returns a named list of model fitted- and predicted values
#'
#' @note Uses mgcv::gam with smooth term s(x) for fitting. If bootstrap.group=TRUE,
#'       results are aggregated across bootstrap iterations.
.fit.dynamics <- function(df,
                          pred.x=seq(min(df$x, na.rm = T), max(df$x, na.rm = T), length.out=50),
                          bootstrap = F,
                          bootstrap.proportion=1,
                          bootstrap.iterations=100,
                          bootstrap.group = T,
                          evaluate.fit = T,
                          ...) {

  if(!bootstrap) {
    bootstrap.proportion = 1
    bootstrap.iterations = 1
  }

  res <- lapply(1:bootstrap.iterations, function(i) {
    # Set subset of ids to use: sample/sub-sample with/without repetition
    ids <- sample(rownames(df), size = ceiling(nrow(df)*bootstrap.proportion), replace = bootstrap)

    # Fit model [y~s(x), weighted by w]
    args <- modifyList(list(formula=y~s(x), data=df[ids, ], weights=df[ids,]$w) , list(...))
    fit  <- do.call(mgcv::gam, args)

    fit.vals <- data.frame(fit$model, predict(fit, df[ids,], se.fit=T), n=nrow(df)) %>% mutate(i=i)
    prd.vals <- data.frame(x=pred.x, i, predict(fit, data.frame(x=pred.x), se.fit=T))


    if(!evaluate.fit)
      return(list(fit = fit.vals, prd = prd.vals))

    # Test fitted model
    models <- list(
      # against null model y~1
      null = do.call(mgcv::gam, modifyList(args, list(formula=y~1, weights=NULL))),
      # against un-weighted y~s(x) (i.e. not trajectory specific)
      unweighted = do.call(mgcv::gam, modifyList(args, list(weights=NULL))))

    evaluations <- do.call(rbind, lapply(names(models), function(m)
      data.frame(comparison=m, anova(fit, models[[m]], test = "F")[2,c("Deviance", "F", "Pr(>F)")], i=i)))

    return(list(fit = fit.vals, prd = prd.vals, evaluations=evaluations))
  })

  # Merge results of different subsampling iterations
  res <- lapply(1:length(res[[1]]), function(i) do.call(rbind, lapply(res, "[[", i)))

  # Summarise subsampling results for fitted- and predicted values
  if (bootstrap.group) {
    for(j in 1:2)
      res[[j]] <- res[[j]] %>%
                    dplyr::select(-i) %>%
                    group_by_at(vars(-fit, -se.fit)) %>%
                    summarise_all(.funs = list(mean=mean, sd=sd)) %>%
      `colnames<-`(gsub("_mean", "", colnames(.)))
  }

  return(res)
}

#' Leiden multiplex clustering over a set of layers
#'
#' This function partitions the vertices of the graphs using python's
#' leidenalg.Optimizer.optimise_partition_multiplex function.
#'
#' It searches for a resolution parameter maximizing the modularity
#' of the CPMVertexPartition quality
#'
#' @param graphs List of graphs over the same set of named vertices
#' @param weights vector of length `graphs` specifying the weight to be used for every layer
#' @param partition.type Type of partition (default: c("CPMVertexPartition", "RBERVertexPartition","RBConfigurationVertexPartition"))
#' @param resolution.range Range for resolution parameter search (default: c(1e-2,10))
#' @param set.seed Random seed (default: 1)
#'
#' @returns Named list of:
#' - membership: named vector with vertices partition
#' - resolution.parameters: vector of chosen resolution parameters corresponding given graphs
#'
#' @note Uses python leidenalg package via reticulate.
optimize.partition <- function(graphs,
                               weights = rep(1, length(graphs)),
                               partition.type=c("CPMVertexPartition", "RBERVertexPartition","RBConfigurationVertexPartition"),
                               resolution.range = c(1e-2,10),
                               set.seed=as.integer(1)) {
  la <- reticulate::import("leidenalg")
  partition.type <- match.arg(partition.type)
  partition.type <- la[[partition.type]]

  opt <- la$Optimiser()
  opt$set_rng_seed(set.seed)

  # Scan for resolution parameters which maximize the modularity
  partitions  <- lapply(graphs, function(g) {
    profile <- opt$resolution_profile(g,
                                      partition.type,
                                      resolution_range = resolution.range,
                                      weights = "weight",
                                      number_iterations = -1)

    res <- profile[[which.max(lapply(profile, function(p) p$modularity))]]$resolution_parameter
    return(partition.type(g, weights="weight", resolution_parameter = res))
  })

  if(length(graphs) == 1)
    opt$optimise_partition(partitions[[1]], n_iterations = -1)
  else
    opt$optimise_partition_multiplex(partitions, layer_weights=weights, n_iterations = -1)

  return(list(
    membership            = partitions[[1]]$membership %>% `names<-`(partitions[[1]]$graph$vs[["name"]]),
    resolution.parameters = lapply(partitions, "[[", "resolution_parameter")))
}


#' Construct communities from dynamics and correlations
#'
#' @description Performs multiplex Leiden clustering on dynamics adjacency and correlation matrices
#'              to identify communities of features with similar trajectory patterns.
#'
#' @param dynamics Data frame with columns: feature, trajectory, x, fit (from fit.dynamics)
#' @param correlations Correlation matrix between features
#' @param k Number of nearest neighbors for local density estimation
#' @param weights Vector of weights for graph layers (default: c(1, -1, 1) for pos, neg, dynamics)
#' @param dynamics.adjacency.clip Threshold to clip small adjacency values (default: 1e-4)
#' @param partition.type Type of partition for Leiden algorithm (default: c("CPMVertexPartition", "RBERVertexPartition","RBConfigurationVertexPartition"))
#'
#' @return Named list with:
#'         - membership: named vector of community assignments
#'         - dynamics.mtx: scaled dynamics matrix
#'         - dynamics.adjacency: adjacency matrix from dynamics
#'         - correlation.mtx: input correlation matrix
#'         - graphs: list of igraph graph objects
#'         - resolution.parameters: chosen resolution parameters
construct.communities <- function(dynamics, correlations, k, weights,
                                  dynamics.adjacency.clip=1e-4,
                                  partition.type=c("CPMVertexPartition", "RBERVertexPartition","RBConfigurationVertexPartition")) {

  dynamics <- dynamics %>% reshape2::dcast(feature~trajectory+x, value.var = "fit") %>%
    column_to_rownames("feature") %>%
    as.matrix()
  # Ensure identical row order for subsequent creating of graph vertices
  correlations <- correlations[rownames(dynamics), rownames(dynamics)]

  # ------------------------------------------ #
  # ---- Create dynamics adjacency matrix ---- #
  # ------------------------------------------ #
  m <- dynamics
  # Scale dynamics matrix
  m <- (m - apply(m, 1, mean, na.rm=T)) / apply(m, 1, sd, na.rm=T)

  # Adjust trajectory weights such that trajectories are weighted equally regardless to their size
  W <- table(gsub("_.*","", colnames(m))) %>% lapply(., function(a) rep(1/a, a)) %>% unlist() %>% diag(.)
  W <- W / min(W[W!=0])

  # Compute Mahalanobis distance
  D <- outer(1:nrow(m), 1:nrow(m), Vectorize(function(i,j) (m[i,] - m[j,]) %*% W %*% (m[i,] - m[j,])))

  # Estimate local density
  s <- apply(D, 1, function(.) sqrt(mean(sort(., partial=k+1)[2:(k+1)])))

  # Compute adaptive-RBF adjacency matrix
  A <- exp(-D/outer(s, s)) %>% `dimnames<-`(list(rownames(m), rownames(m)))
  A[A<dynamics.adjacency.clip] <- 0
  rm(W, D, s)


  # ------------------------------------------ #
  # --- Create graph layers for clustering --- #
  # ------------------------------------------ #
  ig <- reticulate::import("igraph")

  # Create graph layers with edges for positive- and negative correlations
  g <- ig$Graph$Weighted_Adjacency(correlations, loops = F, mode="undirected")
  g$vs$set_attribute_values("name", rownames(correlations))

  g.pos <- g$subgraph_edges(g$es$select(weight_gt = 0), delete_vertices=F)
  g.neg <- g$subgraph_edges(g$es$select(weight_lt = 0), delete_vertices=F)
  g.neg$es$set_attribute_values("weight", -g.neg$es[["weight"]])

  # Create graph layer for dynamics adjacency matrix
  g.dyn <- ig$Graph$Weighted_Adjacency(A, loops = F, mode = "undirected")
  g.dyn$vs$set_attribute_values("name", rownames(A))

  graphs <- list(g.positive.correlation=g.pos,
                 g.negative.correlation=g.neg,
                 g.dynamics=g.dyn)

  # ------------------------------------------ #
  # -- Partition states based on all layers -- #
  # ------------------------------------------ #
  partitions <- optimize.partition(list(g.pos, g.neg, g.dyn), weights=weights, partition.type = partition.type)
  return(list(
    membership            = partitions$membership,
    dynamics.mtx          = m,
    dynamics.adjacency    = A,
    correlation.mtx       = correlations,
    graphs                = graphs,
    resolution.parameters = partitions$resolution.parameters))
}



#'
#' Compute donor-level community proportions
#'
#' @description Converts state-level prevalences to community proportions by weighting
#'              states by their community membership.
#'
#' @param prevalences Matrix with states as columns, donors as rows (state prevalences per donor)
#' @param membership Named vector or data frame mapping states to communities
#'
#' @return Data frame with community proportions per donor (rows sum to 1)
#'
#' @details States are weighted by 1/(community size) to normalize for community size.
donor.community.proportion <- function(prevalences, membership) {
  # state-community one-hot encoding matrix, including any states not found in the state.community named vector
  indic <- membership %>% rownames_to_column() %>%
    `colnames<-`(c("r","p")) %>%
    mutate(v=1) %>%
    reshape2::dcast(r~p, value.var="v", fill=0) %>%
    column_to_rownames("r")
  # state weight in community: 1/(size community)
  indic <- t(t(indic) / Matrix::colSums(indic))
  
  prevalences <- t(t(prevalences) / Matrix::colSums(prevalences))
  proportions <- prevalences %*% indic[colnames(prevalences),]
  return(data.frame(proportions / Matrix::rowSums(proportions)))
}


#' Associate Allen Brain Atlas subtypes with AD traits
#'
#' @description Computes subtype proportions per donor and associates them with AD traits
#'              using associate.traits, controlling for age, sex, and PMI.
#'
#' @param df Data frame with projid and pruned_subtype columns
#' @param exclude_list Character vector of subtype names to exclude (default: empty)
#'
#' @return Data frame of trait association results
#'
associate.allen.traits <- function(df, data, exclude_list=c()) {
  features <- df %>%
    group_by(projid, pruned_subtype) %>%
    summarise(n=n()) %>%
    filter(pruned_subtype != 'NA' & !(pruned_subtype %in% exclude_list)) %>%
    pivot_wider(id_cols = 'projid', names_from = 'pruned_subtype', values_from = n, values_fill = 0) %>%
    column_to_rownames('projid') %>%
    filter(rownames(.) %in% rownames(data)) %>%
    arrange(match(rownames(.), rownames(data)))

  features <- features / rowSums(features)


  trait.association.pruned_subtype <- associate.traits(
    traits = data$obsm$meta.data[,c("sqrt.amyloid_mf","sqrt.tangles_mf", "sqrt.amyloid","sqrt.tangles", "cogng_demog_slope")],
    covariates = sqrt(features),
    controls = data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")]))

  trait.association.pruned_subtype
}


#' Associate within-subtype program metrics with AD traits
#'
#' @description Pivots donor-level program metrics (one row per projid × pruned_subtype)
#'              to wide format (one column per subtype) and runs associate.traits.
#'
#' @param dynamics_features Data frame with projid, pruned_subtype, and program metric columns
#' @param data AnnData object with obsm$meta.data containing trait and control columns
#' @param program_col Name of the program metric column to associate
#' @param exclude_list Character vector of subtype names to exclude (default: empty)
#' @param transform Function applied to covariates before association (default: identity; use sqrt for scores)
#'
#' @return Data frame of trait association results from associate.traits, with a program column
associate.subtype.program.traits <- function(dynamics_features, data, program_col, exclude_list=c(), transform=identity) {
  features <- dynamics_features %>%
    filter(pruned_subtype != 'NA' & !(pruned_subtype %in% exclude_list)) %>%
    dplyr::select(projid, pruned_subtype, all_of(program_col)) %>%
    mutate(pruned_subtype = stringr::str_to_upper(pruned_subtype)) %>%
    pivot_wider(id_cols = 'projid', names_from = 'pruned_subtype', values_from = all_of(program_col), values_fill = NA) %>%
    column_to_rownames('projid') %>%
    filter(rownames(.) %in% rownames(data)) %>%
    arrange(match(rownames(.), rownames(data)))

  traits_cols <- c("sqrt.amyloid_mf", "sqrt.tangles_mf", "sqrt.amyloid", "sqrt.tangles", "cogng_demog_slope")
  traits <- data$obsm$meta.data[, traits_cols] %>% mutate(across(all_of(traits_cols), ~as.numeric(as.character(.))))
  controls <- data.frame(data$obsm$meta.data[, c("age_death", "msex", "pmi")])

  associate.traits(
    traits = traits,
    covariates = transform(features),
    controls = controls
  ) %>% mutate(program = program_col)
}
