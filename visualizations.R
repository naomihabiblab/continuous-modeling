library(data.table)
library(reshape2)
library(ggplot2)
library(viridis)


FeatureDensity <- function(object, features, misc.entry="sim.umap", categorical=F, assay=DefaultAssay(object),
                           slot="data", include.na=F, ncol = 2,
                           return_object=FALSE, ...) {
  # Get similarity and feature values data for samples in the intersection of object, features and similarities
  cells <- data.frame(row.names = colnames(object))[rownames(object@misc[[misc.entry]]),]
  if(is.data.frame(features) | is.matrix(features)) {
    cells <- cells[rownames(features),]
    features <- features[rownames(cells), ]
  } else {
    DefaultAssay(object) <- assay
    # features  <- FetchData(object, vars = features, rownames(cells), slot)
    features  <- FetchData(object, vars = features, colnames(object), slot)
  }

  # If features are categorical, represent values of each feature as an independent feature (dummy variables)
  if(categorical) {
    feature.names <- colnames(features)
    features <- do.call(cbind, lapply(feature.names, function(c) features %>%
                                        dplyr::filter(include.na | !is.na(base::get(c))) %>%
                                        rownames_to_column() %>%
                                        dplyr::rename_at(vars(c), funs(paste0("col"))) %>%
                                        dplyr::mutate(col=paste0(c, ".", col), value=1) %>%
                                        spread(col, value, fill = 0) %>%
                                        dplyr::select_if(! colnames(.) %in% feature.names) %>%
                                        column_to_rownames("rowname")))
  }
  features[is.na(features)] <- 0

  cells <- intersect(rownames(cells), rownames(features))

  sim <- object@misc[[misc.entry]][cells, cells]
  dens <- sim %*% as.matrix(features) / Matrix::rowSums(sim)
  dens[dens < 1e-5] <- 0

  colnames(dens) <- paste0(gsub("-", ".", colnames(dens)), ".", misc.entry)
  object <- AddMetaData(object, dens)
  if (return_object) {
    return(object)
  }
  plots <- FeaturePlot(object, colnames(dens), combine = FALSE, ...)
  for(i in 1:length(plots)) {
    plots[[i]] <- plots[[i]] + NoLegend() + NoAxes()
  }
  return(plot_grid(plotlist = plots, ncol = ncol))
}
