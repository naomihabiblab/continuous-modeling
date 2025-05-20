# This files is base on fit_topic_model.R from commit 65e3ad8a8e6f11ec6e5839e6e76dd14f37fca16d with caching added
library(fastTopics)

CACHE_DIR <- "cache"
dir.create(CACHE_DIR, showWarnings = FALSE)

fit_topics <-
  function(X, k, name, numiter.main = 100, numiter.refine = 100, method.main = "em",
           method.refine = "scd", init.method = c("topicscore", "random"),
           control.init = list(), control.main = list(numiter = 4),
           control.refine = list(numiter = 4, extrapolate = TRUE),
           verbose = c("progressbar", "detailed", "none")) {

    name <- paste(name, k, method.main, method.refine, sep="_")
    message("Using cache prefix: ", name)
    # Check the input data matrix.
    fastTopics:::verify.count.matrix(X)

    # Check and process input argument "verbose".
    verbose <- match.arg(verbose)

    # If necessary, remove all-zero columns from the counts matrix.
    if (fastTopics:::any.allzero.cols(X)) {
      X <- fastTopics:::remove.allzero.cols(X)
      warning(sprintf(paste("One or more columns of X are all zero; after",
                            "removing all-zero columns, %d columns will be",
                            "used for model fitting"), ncol(X)))
    }

    init_fit_path <- file.path(CACHE_DIR, paste0(name, "_init.RDS"))
    # Initialize the Poisson NMF model fit.

    if (file.exists(init_fit_path)) {
      message("Reading init results from cache")
      fit <- readRDS(init_fit_path)
    }
    else {
      fit <- init_poisson_nmf(X, k = k, init.method = init.method,
                            control = control.init,
                            verbose = ifelse(verbose == "none",
                                             "none", "detailed"))
      saveRDS(fit, init_fit_path)
    }

    fitted_fit_path <- file.path(CACHE_DIR, paste0(name, "_fitted.RDS"))

    # Perform the main model fitting step.
    if (file.exists(fitted_fit_path)) {
      message("Reading fitted results from cache")
      fit <- readRDS(fitted_fit_path)
    }
    else {
      fit <- fit_poisson_nmf(X, fit0 = fit, numiter = numiter.main,
                           method = method.main, control = control.main,
                           verbose = verbose)
      saveRDS(fit, fitted_fit_path)
    }


    # Perform the model refinement step.
    if (numiter.refine > 0) {
      refined_fit_path <- file.path(CACHE_DIR, paste0(name, "_refined.RDS"))
      if (file.exists(refined_fit_path)) {
        message("Reading refined results from cache")
        fit <- readRDS(refined_fit_path)
      }
      else {
        if (verbose != "none")
          cat("Refining model fit.\n")
        fit <- fit_poisson_nmf(X, fit0 = fit, numiter = numiter.refine,
                             method = method.refine, control = control.refine,
                             verbose = verbose)
        saveRDS(fit, refined_fit_path)
      }
    }

    # Output the multinomial topic model fit.
    return(poisson2multinom(fit))
}
