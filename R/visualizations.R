#' @title Visualization Functions
#' @description General visualization functions including dot plots,
#'              density plots, and feature visualizations.

library(data.table)
library(reshape2)
library(ggplot2)
library(viridis)



#' Create dot plot for expression data
#'
#' @description Creates a dot plot showing expression level (color) and
#'              percentage of expressing cells (size) for features by group.
#'
#' @param data Data frame with feature columns
#' @param columns Character vector of column names to plot
#' @param group.vector Grouping vector for samples
#' @param order Hierarchically cluster columns (default: TRUE)
#' @param group.levels Optional factor levels for groups
#' @param do.return.order Return column order instead of plot (default: FALSE)
#' @param alpha.threshold Threshold for "expressing" (default: 0.1)
#' @param mean_threshold Threshold for mean calculation (default: 0)
#' @param scale Scale mean values per variable (default: FALSE)
#' @param column.group Optional grouping for columns (facet strips)
#'
#' @return ggplot object or column order vector if do.return.order=TRUE
dot_plot <- function(data, 
                     columns,
                     group.vector,
                     order=TRUE,
                     group.levels=NULL,
                     do.return.order=F,
                     alpha.threshold = 0.1,
                     mean_threshold = 0,
                     scale = FALSE,
                     column.group = NULL # New optional argument for grouping columns (strip in x axis)
                     ){

  alpha <- function(x) {mean(x > alpha.threshold) * 100}
  group.mean <- function(x) {mean(x[x > mean_threshold])}

  columns <- columns[columns %in% colnames(data)]

  # If column.group is supplied, ensure it has one entry per column
  if (!is.null(column.group)) {
    if (length(column.group) != length(columns)) {
      stop("When provided, column.group must have length equal to length(columns)")
    }
    # Allow column.group to be factor or character
    column.group.factor <- as.factor(column.group)
    column.group.df <- data.frame(variable = make.names(columns), column.group = column.group.factor, stringsAsFactors=FALSE)
  }

  t1 <- reshape2::melt(data.frame(data[,columns], groups = as.factor(group.vector)) %>%
                         group_by(groups) %>%
                         summarise_all(~ alpha(.)))
  t2 <- reshape2::melt(data.frame(data[,columns], groups = as.factor(group.vector)) %>%
                         group_by(groups) %>%
                         summarise_all(~ group.mean(.)))

  if (scale) {
    # For each variable (column), scale value.group.mean in t2 to have mean 0 and sd 1
    t2 <- t2 %>%
      group_by(variable) %>%
      mutate(value = as.numeric(scale(value))) %>%
      ungroup()
  }

  data.summarized <- merge(t1,
                           t2,
                           by=c('groups','variable'),
                           suffixes=c('.alpha','.group.mean'))
  data.summarized$variable <- make.names(data.summarized$variable)

  if (!is.null(column.group)) {
    # Attach the group info by variable (category assignment)
    data.summarized <- merge(
      data.summarized,
      column.group.df,
      by="variable",
      all.x=TRUE,
      sort=FALSE
    )
  }

  # Handle variable ordering and keep group/category associativity
  if (order) {
    hclustering <- hclust(as.dist(1 - cor(data[, columns])), method = 'ward.D2')
    ordered_variables <- make.names(columns[hclustering$order])
    data.summarized$variable <- factor(data.summarized$variable, levels = ordered_variables)
    if (!is.null(column.group)) {
      # Identify for each variable what its group is (many variables per group/category)
      col_ordered <- column.group.df[match(ordered_variables, column.group.df$variable), , drop=FALSE]
      data.summarized$column.group <- factor(data.summarized$column.group, levels = unique(col_ordered$column.group))
    }
  } else {
    ordered_variables <- make.names(columns)
    data.summarized$variable <- factor(data.summarized$variable, levels = ordered_variables)
    if (!is.null(column.group)) {
      # "natural" order of categories as in input, preserving multiple variables per group
      data.summarized$column.group <- factor(data.summarized$column.group, levels = unique(column.group))
    }
  }

  data.summarized <- data.summarized[order(data.summarized$variable),]

  if (!is.null(group.levels)) {
    levels(data.summarized$groups) <- group.levels
  }

  base_plot <- ggplot(data.summarized, aes(variable, groups)) +
    geom_point(aes(size = value.alpha, color = value.group.mean), stroke = 0) +
    scale_color_gradientn(colors = c("#00FFFF", "#FF9900")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = paste0('% cells have more than ', alpha.threshold),
         color = paste0('mean x > ', mean_threshold, ifelse(scale, ' (scaled)', '')))

  if (!is.null(column.group)) {
    # Add facet strip by column group/category, allowing several columns in a single group
    P <- base_plot +
      facet_grid(. ~ column.group, scales = "free_x", space = "free_x")
  } else {
    P <- base_plot
  }

  if (do.return.order) {
    return(columns[hclustering$order])
  }
  return(P)
}



#' Plot smoothed feature density on embedding
#'
#' @description Uses a similarity matrix to smooth feature values and
#'              create density-style feature plots.
#'
#' @param object Seurat object
#' @param features Features to plot (vector or data frame)
#' @param misc.entry Name of similarity matrix in object@misc (default: "sim.umap")
#' @param categorical Treat features as categorical (default: FALSE)
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use (default: "data")
#' @param include.na Include NA values in categorical (default: FALSE)
#' @param ncol Number of columns in output (default: 2)
#' @param return_object Return modified object instead of plot (default: FALSE)
#' @param ... Additional arguments passed to FeaturePlot
#'
#' @return Cowplot grid or modified Seurat object
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
