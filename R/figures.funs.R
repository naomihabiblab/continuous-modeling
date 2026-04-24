#####################################################################################################################
#                                           General settings and functions                                          #
#####################################################################################################################
#' @title Figure Functions for Program Analysis Visualization
invisible(lapply(c("reshape2","ggplot2","cowplot","dplyr","tibble","ComplexHeatmap","reticulate","anndata",
                   "patchwork","colorspace","ggrepel","rhdf5", "ggnewscale", "gg3D", "stringr"), library, character.only = T))

no.labs <- labs(x=NULL, y=NULL, color=NULL, fill=NULL)
no.axes <- theme(axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank())
theme_embedding <- theme_classic() + no.axes

embed.width <- embed.height <- 5
embed.width.small <- embed.height.small <- 2.5

blue2red <- colorRampPalette(c("midnightblue","#003BC0","#6689d9","white","#eb8694","#D80D29","firebrick4"))
green2purple <- colorRampPalette(c("darkgreen","white","darkorchid4"))
green2purple.less.white <- colorRampPalette(c("darkgreen","#5F9E5F","white","#A074B6","darkorchid4"))
cyan2orange <- colorRampPalette(c("#4db6ac", "white", "#ef6c00"))
cyan2orange.correlation <- circlize::colorRamp2(c(-1,0,1), c("#4db6ac","white","#ef6c00"))
window2fox <- colorRampPalette(c("#008585", "white", "#c7522a"))

white2orange <- colorRampPalette(c("white", "#ef6c00"))

clay_and_sea_pallet <- colorRampPalette(c(
  "#74a892", "#b8cdab", "#e5c185","#e0a278", "#db836b", "#c7522a", "#ad4623"
), bias=0.4)


#####################################################################################################################
#                                                 Plotting Functions                                                #
#####################################################################################################################

#' Create a customizable heatmap or dot plot from a data frame
#'
#' @description Converts a long-format data frame into a heatmap using ComplexHeatmap.
#'              Supports optional dot plot mode when size.by is specified.
#'
#' @param df Data frame in long format with at least 3 columns (rows, columns, values)
#' @param rows.by Column index, name, or expression for row identifiers (default: 1)
#' @param columns.by Column index, name, or expression for column identifiers (default: 2)
#' @param color.by Column index, name, or expression for color values (default: 3)
#' @param size.by Optional column index, name, or expression for dot sizes (creates dot plot)
#' @param scale.by Scaling method: "" (none), "row", or "col" (default: "")
#' @param row.order Named list for custom row ordering and splitting
#' @param column.order Named list for custom column ordering and splitting
#' @param cols Color palette function (default: cyan2orange)
#' @param size.factor Scaling factor for dot sizes (default: 1.0)
#' @param size.by.round Decimal places for size legend labels (default: 2)
#' @param annotations Function(rows, cols) returning annotation list
#' @param row_labels Function to transform row labels
#' @param column_labels Function to transform column labels
#' @param max_col Maximum value for color scale (default: NULL, auto-computed)
#' @param min_col Minimum value for color scale (default: NULL, auto-computed)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return Heatmap object, or list(hm, legend_list) if size.by is specified
gheatmap <- function(df, rows.by=1, columns.by=2, color.by=3, size.by=NULL,
                     scale.by=c("","row","col"),
                     row.order=NULL, column.order=NULL,
                     cols=cyan2orange,
                     size.factor=1.,
                     size.by.round=2,
                     annotations=function(rows, cols) list(),
                     row_labels=function(rows) rows,
                     column_labels=function(cols) cols,
                     max_col=NULL,
                     min_col=NULL,
                     ...) {
  scale.by = match.arg(scale.by)
  #TODO - support multiple gene lists and/or repeated gene names

  # Add columns as specified by function variables being a column index, name or expression
  for(p in c("rows.by", "columns.by")) {
    if (is.numeric(get(p)))
      assign(p, colnames(df)[get(p)])
    df <- dplyr::rename(df, {{p}}:=get(p))
  }

  # For color/size parameters support the use of expressions over one or many of the dataframes columns
  for(p in c("color.by","size.by")) {
    if(!is.null(get(p))) {
      if (is.numeric(get(p)))
        assign(p, colnames(df)[get(p)])

      used.cols <- names(base::Filter(isTRUE, sapply(colnames(df), grepl, get(p))))
      df <- df %>% dplyr::rowwise() %>% dplyr::mutate_at(used.cols, ~eval(parse(text=.)))
      df <- dplyr::mutate(df, {{p}} := eval(parse(text=get(p))))
    }
  }
  rm(p)

  # Create values matrix
  vals <- tidyr::pivot_wider(df, id_cols = "rows.by", names_from = "columns.by", values_from = "color.by", values_fill = 0) %>% tibble::column_to_rownames("rows.by") %>% as.matrix()

  # Apply row/column ordering
  if(!is.null(row.order)) {
    row.order <- sapply(row.order, function(v) intersect(v, rownames(vals)))
    vals <- vals[unlist(row.order),]
  }
  if(!is.null(column.order)) {
    column.order <- sapply(column.order, function(v) intersect(v, colnames(vals)))
    vals <- vals[,unlist(column.order)]
  }

  # Scaling of matrices
  if(scale.by == "col") vals <- scale(vals)
  if(scale.by == "row") vals <- t(scale(t(vals)))

  # Specify color scale for plot
  if(scale.by == "")
    col.vals <- seq(min_col %||% min(vals, na.rm=T), max_col %||% max(vals, na.rm=T), length.out=21)
  else {
    v <- max(abs(vals), na.rm = T)
    col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])
  }

  col.fun  <- circlize::colorRamp2(col.vals, cols(length(col.vals)))

  # If `size` is specified create dotplot instead of heatmap
  if(!is.null(size.by)) {
    sizes <- tidyr::pivot_wider(df, id_cols = "rows.by", names_from = "columns.by", values_from = "size.by", values_fill = 0) %>% column_to_rownames("rows.by") %>% as.matrix()
    sizes <- sizes / max(sizes)

    if(!is.null(row.order)) sizes <- sizes[unlist(row.order),]
    if(!is.null(column.order)) sizes <- sizes[,unlist(column.order)]

    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x=x, y=y, width=w, height=h, gp=gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r=pindex(sizes, i, j)*size.factor * unit(2, "mm"), gp = gpar(fill = col.fun(pindex(vals, i, j)), col = NA))}
    rect_gp = gpar(type = "none")

    legend.list = list(
      Legend(labels = round(seq(min(df[,"size.by"], na.rm = T), max(df[,"size.by"], na.rm = T), length.out=5), size.by.round), title = size.by,
             graphics = lapply(c(.01, .25, .5, .75, 1), function(p)
               function(x, y, w, h) grid.circle(x = x, y = y, r = p * unit(2, "mm"), gp = gpar(fill = "black"))),
             direction = "horizontal"))
  } else {
    layer_fun = function(j, i, x, y, w, h, fill) {}
    rect_gp = gpar(col = NA)
    legend.list = list()
  }

  args <- modifyList(modifyList(list(
    cluster_rows = T,
    cluster_columns = T,
    row_labels=row_labels(rownames(vals)),
    column_labels=column_labels(colnames(vals))),
                                annotations(rownames(vals), colnames(vals))),
                     list(...))
  if(!is.null(names(row.order)) & !hasName(args, "row_split")) args$row_split = stack(row.order)[,2]
  if(!is.null(names(column.order)) & !hasName(args, "column_split")) args$column_split = stack(column.order)[,2]

  # if(!is.null(mark.rows)) {
  #   mark.rows <- intersect(mark.rows, rownames(vals))
  #   ids <- which(rownames(vals) %in% mark.rows)
  #   ann <- anno_mark(at = ids, labels = rownames(vals)[ids], which = "row")
  #   if(hasName(args, "right_annotation")) args$right_annotation$named_rows <- ann
  #   else                                  args$right_annotation <- rowAnnotation(named_rows=ann)
  # }
  #
  # if(!is.null(mark.cols)) {
  #   mark.cols <- intersect(mark.cols, colnames(vals))
  #   ids <- which(colnames(vals) %in% mark.cols)
  #   ann <- anno_mark(at = ids, labels = colnames(vals)[ids], which = "column")
  #   if(hasName(args, "bottom_annotation")) args$bottom_annotation$named_rows <- ann
  #   else                                   args$bottom_annotation <- HeatmapAnnotation(named_rows=ann)
  # }

  hm <- do.call(Heatmap,
                modifyList(list(matrix=vals, col=col.fun, name=color.by,
                                row_order=unlist(row.order),
                                column_order=unlist(column.order),
                                rect_gp=rect_gp, layer_fun=layer_fun), args))

  if(!is.null(size.by))
    return(list(hm=hm, legend_list=legend.list))
  return(hm)
}


#' Heatmap of trait associations with significance annotations
#'
#' @description Creates a heatmap showing trait associations with t-statistics as colors
#'              and significance levels as text annotations (*, **, ***, ****).
#'
#' @param df Data frame with trait association results
#' @param params Character vector of trait names to include (default: amyloid, tangles, cognitive slope)
#' @param row.by Column name for row identifiers (default: "state")
#' @param col.by Column name for column identifiers (default: "trait")
#' @param value.by Column name for heatmap values (default: "tstat")
#' @param pval.by Column name for p-values (default: "adj.pval")
#' @param cols Color palette function (default: window2fox)
#' @param show.only.significant If TRUE, filter to rows with at least one significant association (default: TRUE)
#' @param plot If TRUE, draw the heatmap; if FALSE, return heatmap object (default: TRUE)
#' @param show_row_dend Show row dendrogram (default: FALSE)
#' @param row_order Optional vector specifying row order
#' @param column_labels Optional named vector for custom column labels
#' @param max_color_value Maximum absolute value for color scale (default: NULL, auto-computed)
#' @param clustering_distance_rows Distance metric for row clustering (default: "euclidean")
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return If plot=TRUE, draws heatmap and returns invisibly. If plot=FALSE, returns list(hm, legend).
plot.trait.associations.heatmap <- function(df, params=c("sqrt.amyloid_mf","sqrt.tangles_mf","cogng_demog_slope"),
                                    row.by="state", col.by="trait", value.by="tstat", pval.by="adj.pval",
                                    cols=window2fox,
                                    show.only.significant=T, plot=T, 
                                    show_row_dend=F,
                                    row_order=NULL,
                                    column_labels=NULL,
                                    max_color_value=NULL,
                                    clustering_distance_rows="euclidean", ...) {

  if("pandas.core.frame.DataFrame" %in% class(df))
    df <- py_to_r(df)

  df <- df %>% filter(trait %in% params)

  if(show.only.significant)
    df <- df %>% group_by_at(c(row.by)) %>%
      filter(if_any(pval.by, ~sum(. < .05) > 0)) %>%
      ungroup() %>%
      data.frame(check.names=FALSE)

  df$sig <- cut(df[[pval.by]], c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", ""))

  if (nrow(df) == 0) {
    print(grid::grid.text("No trait associated"))
    return()
  }

  vals <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = all_of(value.by), values_fill = NA_real_) %>%
    tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params))
  sig <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = "sig", values_fill = "") %>%
    tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params))

  if (!is.null(row_order)) {
    vals <- vals %>% arrange(match(make.names(rownames(.)), make.names(row_order)))
    sig <- sig %>% arrange(match(make.names(rownames(.)), make.names(row_order)))
  }

  vals <- vals %>% as.matrix()
  sig <- sig %>% as.matrix()

  if (is.null(max_color_value)) {
    max_color_value <- max(abs(vals), na.rm = T)
  }
  col.vals <- c(seq(-max_color_value, 0, length.out=11), seq(0, max_color_value, length.out=11)[-1])

  if (!is.null(column_labels) & length(column_labels) != ncol(vals)) {
    column_labels <- column_labels[colnames(vals)]
  } else if (is.null(column_labels)) {
    column_labels <- colnames(vals)
  }

  hm <- Heatmap(vals,
                name = value.by,
                col = circlize::colorRamp2(col.vals, cols(length(col.vals))),
                cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                cluster_columns = F,
                clustering_distance_rows = clustering_distance_rows,
                show_row_dend = show_row_dend,
                column_labels = column_labels,
                ...)
  legend <- Legend(title = "FDR", pch = c("****","***","**","*"), type = "points", labels = c("<0.0001","<0.001","<0.01", "<0.05"))
  if(!plot)
    return(list(hm=hm, legend=legend))
  draw(hm, annotation_legend_list = list(legend), merge_legend=T)
}

#' Plot 2D embedding landscape with feature coloring
#'
#' @description Creates scatter plots of 2D embeddings (e.g., UMAP, PHATE) colored by
#'              one or more features. Supports smoothing via similarity matrix.
#'
#' @param features Character vector of feature names, or data frame of features, or NULL for uniform coloring
#' @param embedding Name of embedding in data$obsm (default: "X_phate") or a matrix
#' @param cols Color palette function (default: clay_and_sea_pallet)
#' @param shape.by Optional column name for point shapes
#' @param size Point size (default: 2)
#' @param smoothened If TRUE, smooth values using similarity matrix (default: TRUE)
#' @param enforce.same.color.scale If TRUE, use same color scale across all features (default: TRUE)
#' @param sort.direction Direction for sorting points (1 or -1) to control overplotting (default: 1)
#' @param legend.position Legend position (default: "none")
#' @param raster If TRUE, use rasterized points for smaller file size (default: FALSE)
#' @param ncol Number of columns in plot grid (default: NULL)
#' @param nrow Number of rows in plot grid (default: NULL)
#' @param data. Data object containing X, obsm, obs, and obsp (default: data)
#'
#' @return Combined plot from cowplot::plot_grid()
plot.landscape <- function(
    features,
    embedding = "X_phate",
    cols = clay_and_sea_pallet,
    shape.by = NULL,
    size = 2,
    smoothened = TRUE,
    enforce.same.color.scale = TRUE,
    sort.direction = 1,
    legend.position="none",
    raster=FALSE,
    ncol=NULL,
    nrow=NULL,
    data. = data) {

  if(is.null(features)) {
    features <- data.frame(c=rep(1,nrow(data.)))
    cols <- function(n) rep("black", n)
    smoothened <- FALSE
  } else if (class(features) == "character") {
    features <- data.frame(data.$X, data.$obsm$meta.data, data.$obs) %>% dplyr::select_at(all_of(features))
  } else if (is.null(ncol(features))) {
    features <- data.frame(features)
  }

  if(class(embedding) == "character") {
    embedding.str = embedding
    embedding <- data.$obsm[[embedding]] %>% `colnames<-`(c("x","y"))
    embedding.axes <- colnames(embedding)
  }

  df <- data.frame(embedding, features)
  if(!is.null(shape.by))
    df$shape.by = data.$obsm[,shape.by]

  df <- df %>% filter(!if_all(embedding.axes, is.na))
  if(smoothened) {
    msk <- rownames(df)
    sim <- data.$obsp[[paste0("similarity_", embedding.str)]] %>%
      `dimnames<-`(list(data.$obs_names, data.$obs_names)) %>%
      `[`(msk, msk)
    for(p in colnames(features)) {
      if (class(features[,p]) == "numeric")
        df[!is.na(df[, p]), p] <- (sim[!is.na(df[, p]),!is.na(df[, p])] %*% matrix(df[!is.na(df[, p]), p]))
    }
  }

  limits <- NULL
  if(enforce.same.color.scale & ncol(features) > 1) {
    cmin <- min(apply(df[,colnames(features)], 1, min, rm.na=T))
    cmax <- max(apply(df[,colnames(features)], 1, max, rm.na=T))
    limits <- c(cmin, cmax)
  }

  lapply(1:ncol(features), function(i) {
    if(class(features[,i]) == "numeric")
      colors <- scale_color_gradientn(colors = cols(21), na.value = "lightgrey", limits = limits)
    else
      colors <- scale_color_manual(values=cols(length(unique(features[,i]))))

    args <- lapply(list(x="x", y="y", color=colnames(features)[[i]], shape=shape.by), function(x) if (!is.null(x)) sym(x))
    df. <- df
    if (class(features[,i]) == "numeric")
      df. <- df. %>% arrange(!is.na(!!sym(args$color)), sort.direction*!!sym(args$color))

    plt  <- ggplot(df., aes(,,!!!args))

    if(raster) plt <- plt + ggrastr::geom_point_rast(size=size, raster.dpi=600)
    else       plt <- plt + geom_point(size=size)

    plt + colors + theme_embedding + no.labs + labs(title=args$color) +
      scale_y_continuous(expand=expansion(add=.004)) +
      scale_x_continuous(expand=expansion(add=.004)) +
      theme(legend.position = legend.position,
            panel.background = element_blank())
  }) %>% plot_grid(plotlist = ., nrow=nrow, ncol=ncol)
}



#' Plot 3D embedding landscape with feature coloring
#'
#' @description Creates 3D scatter plots using gg3D package for 3D embeddings
#'              colored by one or more features.
#'
#' @param features Character vector of feature names, data frame, or NULL
#' @param embedding Name of 3D embedding in data$obsm (default: "X_all_3d_phate")
#' @param cols Color palette function (default: scales::hue_pal())
#' @param shape.by Optional column for point shapes
#' @param size Point size (default: 2)
#' @param smoothened Smooth values using similarity matrix (default: TRUE)
#' @param enforce.same.color.scale Use same color scale across features (default: TRUE)
#' @param sort.direction Sorting direction for overplotting (default: 1)
#' @param legend.position Legend position (default: "none")
#' @param theta Viewing angle theta for 3D projection (default: 0)
#' @param phi Viewing angle phi for 3D projection (default: 0)
#' @param ncol Number of columns in plot grid
#' @param nrow Number of rows in plot grid
#' @param data. Data object with obsm and obsp (default: data)
#'
#' @return Combined plot from cowplot::plot_grid()
plot.landscape.3D <- function(
    features = NULL,
    embedding = "X_all_3d_phate",
    cols = scales::hue_pal(),
    shape.by = NULL,
    size = 2,
    smoothened = TRUE,
    enforce.same.color.scale = TRUE,
    sort.direction = 1,
    legend.position="none",
    theta=0,
    phi=0,
    ncol=NULL,
    nrow=NULL,
    data. = data) {

  if(is.null(features)) {
    features <- data.frame(c=rep(1,nrow(data.)))
    cols <- function(n) rep("black", n)
    smoothened <- FALSE
  } else { if(class(features) == "character")
    features <- data.frame(data.$X, data.$obsm$meta.data, data.$obs) %>% dplyr::select_at(all_of(features))
  }

  if(class(embedding) == "character") {
    embedding.str = embedding
    embedding <- data.$obsm[[embedding]] %>% `colnames<-`(c("x","y","z"))
    embedding.axes <- colnames(embedding)
  }

  df <- data.frame(embedding, features)
  if(!is.null(shape.by))
    df$shape.by = data.$obsm[,shape.by]

  df <- df %>% filter(!if_all(embedding.axes, is.na))
  if(smoothened) {
    msk <- rownames(df)
    sim <- data.$obsp[[paste0("similarity_", embedding.str)]] %>%
      `dimnames<-`(list(data.$obs_names, data.$obs_names)) %>%
      `[`(msk, msk)
    for(p in colnames(features)) {
      if (class(features[,p]) == "numeric")
        df[!is.na(df[, p]), p] <- (sim[!is.na(df[, p]),!is.na(df[, p])] %*% matrix(df[!is.na(df[, p]), p]))
    }
  }

  limits <- NULL
  if(enforce.same.color.scale & ncol(features) > 1) {
    cmin <- min(apply(df[,colnames(features)], 1, min, rm.na=T))
    cmax <- max(apply(df[,colnames(features)], 1, max, rm.na=T))
    limits <- c(cmin, cmax)
  }

  lapply(1:ncol(features), function(i) {

    if(class(features[,i]) == "numeric")
      colors <- scale_color_gradientn(colors = cols(21), na.value = "lightgrey", limits = limits)
    else
      colors <- scale_color_manual(values=cols(length(unique(features[,i]))))

    args <- lapply(list(x="x", y="y",z="z", color=colnames(features)[[i]], shape=shape.by), function(x) if (!is.null(x)) sym(x))
    ggplot(df %>% arrange(!is.na(!!sym(args$color)), sort.direction*!!sym(args$color)), aes(,,!!!args)) +
      axes_3D(theta=theta, phi=phi) +
      stat_3D(theta=theta, phi=phi, size=size) +
      labs_3D() +
      colors + labs(x=NULL, y=NULL, title=args$color) +
      theme_classic() +
      theme(legend.position = legend.position,
            panel.background = element_blank(),
            axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()
            )
  }) %>% plot_grid(plotlist = ., nrow=nrow, ncol=ncol)
}

#' Wrapper function for plot.dynamics with data preprocessing
#'
#' @description Filters and prepares dynamics data for plotting with plot.dynamics().
#'              Handles both R data frames and Python pandas DataFrames.
#'
#' @param dynamics List with $fitted.vals and $pred.vals data frames
#' @param features Character vector of feature names to plot
#' @param trajectories Character vector of trajectories to include (default: all unique)
#' @param ... Additional arguments passed to plot.dynamics()
#'
#' @return ggplot object from plot.dynamics()
plot.dynamics.wrapper <- function(dynamics, features, trajectories = dynamics$fitted.vals %>% pull(trajectory) %>% unique(), ...) {
  if("pandas.core.frame.DataFrame" %in% class(dynamics$fitted.vals))
    dynamics$fitted.vals <- py_to_r(dynamics$fitted.vals)

  if("pandas.core.frame.DataFrame" %in% class(dynamics$pred.vals))
    dynamics$pred.vals <- py_to_r(dynamics$pred.vals)

  stopifnot(length(features) > 0)

  fits <- dynamics$fitted.vals %>% filter(feature %in% features) %>% filter(trajectory %in% trajectories) %>% dplyr::arrange(match(feature, features))
  pred <- dynamics$pred.vals %>% filter(feature %in% features) %>% filter(trajectory %in% trajectories) %>% dplyr::arrange(match(feature, features))

  stopifnot(nrow(fits) > 0)
  
  n_found <- fits$feature %>% n_distinct()
  if (n_found != length(features)) {
    warning(sprintf("Only %d of %d requested features found in dynamics data", n_found, length(features)))
  }

  return(plot.dynamics(fits, pred, ...))
}


#' Plot dynamics curves with confidence intervals and optional raw points
#'
#' @description Creates publication-quality plots of fitted dynamics curves with
#'              confidence ribbons, optional raw data points, and customizable styling.
#'
#' @param fits Data frame with fitted values (columns: x, y, fit, feature, trajectory, X.weights.)
#' @param preds Data frame with predictions (columns: x, fit, se.fit, feature, trajectory)
#' @param cols Named vector of colors for each group (default: NULL, auto-generated)
#' @param facet.by Faceting variable: "trajectory" or "feature" (default: "trajectory")
#' @param overlap.pseudotime Optional pseudotime value to highlight overlap region
#' @param include.points If TRUE, show raw data points (default: TRUE)
#' @param min.point.alpha Minimum alpha for points (default: 0.1)
#' @param label If TRUE, add labels to curves (default: TRUE)
#' @param legend.position Legend position (default: "right")
#' @param raster If TRUE, rasterize points (default: FALSE)
#' @param fold.change If TRUE, show fold-change from baseline (default: FALSE)
#' @param show.label.legend Show legend for labels (default: FALSE)
#' @param ymin, ymax Y-axis limits (default: NA, auto-computed)
#' @param labels Named vector for label text replacement
#' @param label_size Size of label text (default: 4)
#' @param line_size Line width (default: 0.5)
#' @param ribbon_alpha Alpha for confidence ribbon (default: 0.2)
#' @param show_n If TRUE and 'n' column exists, show sample size (default: TRUE)
#' @param background.features Features to plot in background (default: empty)
#' @param ... Additional arguments passed to facet_wrap()
#'
#' @return ggplot object
plot.dynamics <- function(fits, preds, cols = NULL, facet.by = c("trajectory","feature"),
                          overlap.pseudotime=NULL, include.points=T, min.point.alpha=.1, label=T,
                          legend.position = "right", raster=F, fold.change=F,
                          show.label.legend = F, ymin=NA, ymax=NA, labels=NULL,
                          label_size=4, line_size=0.5, ribbon_alpha=0.2,
                          show_n=TRUE, background.features = c(),
                          ...) {
  facet.by = match.arg(facet.by)
  group.by = setdiff(c("trajectory","feature"), facet.by)

  groups <- unique(fits[[group.by]])
  if(is.null(cols)) {
    cols <- setNames(scales::hue_pal()(length(groups)), groups)
  } else if (0 == length(names(cols)) && length(cols) == length(groups)) {
    names(cols) <- groups
  } else {
    unset_cols <- setdiff(groups, names(cols))
    if (length(unset_cols) > 0) {
      cols[unset_cols] <- scales::hue_pal()(length(unset_cols))
    }
  }

  if(fold.change) {
    fits <- fits %>% group_by(trajectory, feature) %>% mutate(across(c(y,fit), ~./.[x == min(x)]))
    preds <- preds %>% group_by(trajectory, feature) %>% mutate(across(c(fit), ~./.[x == min(x)]))
  }

  if (is.na(ymax)) ymax <- max(preds$fit+2*preds$se.fit, na.rm = T)
  if (is.na(ymin)) ymin <- min(preds$fit-2*preds$se.fit, na.rm = T)

  if(include.points) {
    ymax <- max(c(fits$y, ymax), na.rm = T)*1.05
    ymin <- min(c(fits$y, ymin), na.rm = T)*.95
  }

    p <- ggplot(fits) +
    facet_wrap(paste0("~",facet.by), ...) +
    scale_x_continuous(expand = c(0,0), labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ymax), labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = legend.position,
          axis.line = element_line())

  if(!is.null(overlap.pseudotime)) {
    p <- p +
      geom_ribbon(aes(x,ymax=ymax,ymin=ymin),
                  data.frame(x=c(0, overlap.pseudotime), ymax=rep(ymax, 2), ymin=rep(ymin,2)),
                  fill="black", alpha=.04) +
      geom_vline(xintercept = overlap.pseudotime, linetype="dashed")
  }

  if(include.points) {
    for(g in groups) {
      if(raster)
        p <- p + ggrastr::geom_point_rast(aes(x,y, color=X.weights., alpha=X.weights.), fits[fits[,group.by] == g,], inherit.aes = F, size=1, raster.dpi = 600)
      else
        p <- p + geom_point(aes(x,y, color=X.weights., alpha=X.weights.), fits[fits[,group.by] == g,], inherit.aes = F, size=1)

      p <- p + scale_color_gradientn(colors=colorRampPalette(c("white", cols[[g]]))(10)[-1]) +
        scale_alpha_binned(range=c(min.point.alpha,1), breaks=c(.25, .5, .75, 1), limits=c(0,1)) +
        labs(color=g, alpha=g)  +
        new_scale("alpha") +
        new_scale("color") +
        guides(color_new=guide_legend(g), alpha_new=guide_legend(g))
    }
  }

  for(g in groups) {
    p <- p + geom_ribbon(aes_string("x", "fit", ymin="fit-2*se.fit", ymax="fit+2*se.fit",
                                    fill=group.by),
                         preds[preds[, group.by] == g,], alpha=ribbon_alpha, color=NA, inherit.aes = F) +
      geom_line(aes_string("x", "fit", color=group.by),
                preds[preds[, group.by] == g,], size=line_size, inherit.aes = F) +
      scale_color_manual(values = cols[[g]]) +
      scale_fill_manual(values = cols[[g]]) +
      new_scale("color") +
      new_scale("fill") +
      guides(color_new=guide_legend(g), fill_new=guide_legend(g))

    if ('n' %in% colnames(fits) & show_n) {
      p <- p + geom_text(aes(x=Inf, y=Inf, label=glue::glue("n={n}")),
                colour="black", inherit.aes=FALSE, parse=FALSE, vjust=1, hjust=1)
    }
  }

  if(label)
    p <- p + ggrepel::geom_label_repel(aes_string("x", "fit", label='label', color=group.by), size=label_size,
                                       fits %>%
                                         group_by(feature, trajectory) %>%
                                         slice_min(x, n=1, with_ties = FALSE) %>%
                                         mutate(label=do.call(recode, c(list(.data[[group.by]]), labels, `filler`='red',
                                                                        `.default`=as.character(.data[[group.by]])))),
                                       min.segment.length = unit(0,"pt"), show.legend = show.label.legend) +
    scale_color_manual(values = cols)

  return(p)
}



#' Plot UMAP embedding colored by density of a feature
#'
#' @description Creates a rasterized UMAP plot colored by the density values
#'              of a specified column using a turbo color scale.
#'
#' @param embedding Data frame with x, y (or UMAP_1, UMAP_2) coordinates and feature columns
#' @param column Name of the column to color by
#' @param title Plot title (default: NA, auto-generated from column name)
#' @param title_size Title font size (default: 11)
#'
#' @return ggplot object
#'
#' @note Automatically converts UMAP_1/UMAP_2 column names to x/y.
plot_umap_density <- function(embedding, column, title = NA, title_size=11) {
  column_name <- gsub("X.\\d+.k", "p", column) %>% gsub(".sim.umap", "", .)
  title <- title %||% 'program {column_name} density'
  if ('UMAP_1' %in% colnames(embedding)) {
    embedding <- embedding %>% rename(x=UMAP_1, y=UMAP_2)
  }

  color_lab <- glue::glue('program {column_name} density')
  if (is.na(title)) {
    title <- color_lab
  }


  ggplot(embedding, aes_string('x', 'y', color=column)) +
      ggrastr::geom_point_rast(size=0.01, raster.dpi = 512) +
      scale_color_viridis_c(option='turbo') +
      theme_embedding +
      theme_void() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=title_size)) +
      labs(color="density")
}

#' Plot multiple UMAP density features in a faceted grid
#'
#' @description Creates a faceted plot with multiple density features shown
#'              as separate panels on the same UMAP embedding.
#'
#' @param embedding Data frame with coordinates and feature columns
#' @param columns Character vector of column names to plot
#' @param title Overall plot title (default: "")
#' @param ncol Number of facet columns (default: 2)
#' @param labels Named vector for facet label replacement (default: identity mapping)
#'
#' @return ggplot object with faceted panels
plot_umap_density2 <- function(embedding, columns, title = "", ncol=2, labels=setNames(columns, columns)) {
  if ('UMAP_1' %in% colnames(embedding)) {
    embedding <- embedding %>% dplyr::rename(x=UMAP_1, y=UMAP_2)
  }

  embedding_long <- embedding %>%
    pivot_longer(
      cols = all_of(columns),
      names_to = "variable",
      values_to = "value"
    )

  p <- ggplot(embedding_long, aes(x = x, y = y, color = value)) +
    ggrastr::geom_point_rast(size = 0.01, raster.dpi = 800) +
    scale_color_viridis_c(option = 'turbo') +
    theme_void() +
    ggtitle(title) +
    labs(color = "density") +
    theme(legend.position = "bottom")

  p <- p + facet_wrap(~ variable, scales = "free_y", ncol=ncol, labeller = labeller(variable=labels))

  p
}

#' Plot UMAP with state labels
#'
#' @description Creates a UMAP plot colored by categorical state with
#'              centroid labels for each state.
#'
#' @param embedding Data frame with coordinates and state column
#' @param column Name of the state column (default: "state")
#' @param title Plot title (default: "states")
#'
#' @return ggplot object
#'
#' @note State labels are placed at the centroid of each state's points.
plot_umap_state <- function(embedding, column='state', title = "states") {
  if ('UMAP_1' %in% colnames(embedding)) {
    embedding <- embedding %>% rename(x=UMAP_1, y=UMAP_2)
  }

  ggplot(embedding, aes_string('x', 'y', color=column, label=column)) +
    ggrastr::geom_point_rast(size=.01, raster.dpi = 800) +
    geom_text(data = embedding %>% group_by(!!sym(column)) %>%
                summarise_at(vars(x, y), list(mean)), color="black", size=4) +
    no.labs +
    theme_embedding +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none",)
}

#' Create a transformation that compresses a range of an axis
#'
#' @description Creates a custom scale transformation that compresses (squashes)
#'              a specified range by a given factor. Useful for de-emphasizing
#'              certain value ranges while maintaining continuity.
#'
#' @param from Left end of the range to squash
#' @param to Right end of the range to squash
#' @param factor Compression factor (values > 1 compress the range)
#'
#' @return A scales::trans_new() transformation object
#'
#' @examples
#' # Squash the range 10-100 by factor of 10:
#' # ggplot(...) + scale_x_continuous(trans = squash_axis(10, 100, 10))
squash_axis <- function(from, to, factor) {
  # A transformation function that squashes the range of [from, to] by factor on a given axis

  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function

  trans <- function(x) {
    if (is.na(x)) {
      return(x)
    }
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to

    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)

    return(x)
  }

  inv <- function(x) {
    if (is.na(x)) {
      return(x)
    }

    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor

    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))

    return(x)
  }

  # return the transformation
  return(scales::trans_new("squash_axis", trans, inv))
}

#' Create heatmap of gene dynamics across pseudotime
#'
#' @description Generates a ComplexHeatmap showing gene expression dynamics
#'              across pseudotime for one or both trajectory branches (prAD, ABA).
#'
#' @param dynamics List with $pred.vals containing predicted values (output of fit.dynamics for example)
#' @param genes Character vector of gene names to include
#' @param title Plot title (default: "Gene dynamics")
#' @param scale If TRUE, scale each gene's values (default: TRUE)
#' @param both_trajectories If TRUE, show both prAD and ABA trajectories (default: FALSE)
#' @param cluster_rows If TRUE, cluster rows hierarchically (default: TRUE)
#' @param col Color palette vector (default: RdYlBu reversed)
#' @param show_lfc If TRUE, show log fold-change annotation (default: TRUE)
#' @param name Legend name (default: "gene dynamics scaled")
#' @param row_names_mapping_func Optional function to transform row names
#' @param clustering_distance_rows Distance metric for clustering (default: "spearman")
#' @param column_split Custom column split factor
#' @param row_order Optional custom row order
#' @param use_slanter_ordering If TRUE, use slanter package for optimal ordering (default: FALSE)
#' @param row_split Factor for splitting rows into groups
#' @param in_community_var Data frame with community and sub-community information
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return ComplexHeatmap object
#'
#' @note Requires slanter package if use_slanter_ordering=TRUE.
plot.genes.dynamics.heatmap <- function (dynamics, genes, title="Gene dynamics",
                                         scale=TRUE, both_trajectories=FALSE, cluster_rows=TRUE,
                                         col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7,
                                                                                             name = "RdYlBu")))(100),
                                         show_lfc=TRUE,
                                         name="gene dynamics scaled",
                                         row_names_mapping_func=NULL,
                                         clustering_distance_rows='spearman',
                                         column_split=NULL,
                                         row_order=NULL,
                                         use_slanter_ordering=FALSE,
                                         row_split=NULL,
                                         in_community_var=NULL,
                                         ...) {


  predicted.vals.to.matrix <- function(pred.vals, trajectory) {
    if("pandas.core.frame.DataFrame" %in% class(pred.vals))
      pred.vals <- py_to_r(pred.vals)

    pred.vals %>%
      filter(feature %in% genes) %>%
      filter(.data[['trajectory']] == .env[['trajectory']]) %>%
      dplyr::select(c("x","feature","fit")) %>%
      pivot_wider(names_from = x, values_from = fit) %>%
      as.data.frame %>%
      column_to_rownames(., var = "feature") %>%
      dplyr::select(order(as.numeric(colnames(.)))) %>%
      dplyr::rename_with(~ round(as.numeric(.x),3))
  }

  prAD.mat <- predicted.vals.to.matrix(dynamics$pred.vals, "prAD")
  aba.mat <- predicted.vals.to.matrix(dynamics$pred.vals, "ABA")

  colnames(prAD.mat)[2:(ncol(prAD.mat)-1)] <- ""
  colnames(aba.mat)[2:(ncol(aba.mat)-1)] <- ""

  prAD.time.points <- ncol(prAD.mat)
  aba.time.points <- ncol(aba.mat)

  if (both_trajectories) {
    mat <- cbind(prAD.mat, aba.mat)
    colnames(mat)[startsWith(colnames(mat), "Var.")] <- ""
  } else {
    mat <- prAD.mat
    aba.time.points <- 0
  }

  if (scale) {
    mat <- mat %>% t %>% scale() %>% t
  }

  if (show_lfc) {
    prAD.lfc <- log(apply(prAD.mat + 1e-6, 1, max) / apply(prAD.mat + 1e-6, 1, min), 2)
    lfc_annotation <- rowAnnotation(LFC=prAD.lfc, col=c(LFC=circlize::colorRamp2(c(0, max(prAD.lfc)), c("white", "red"))))
  } else {
    lfc_annotation = NULL
  }

  if (!is.null(row_names_mapping_func)) {
    rownames(mat) <- row_names_mapping_func(rownames(mat))
  }


  column_split <- factor(c(rep('prAD', prAD.time.points), rep('ABA', aba.time.points)),
                         levels = c('prAD', 'ABA'))
  max_value <- max(abs(mat))

  col <- circlize::colorRamp2(seq(-max_value, max_value, length.out = length(col)),
                              colors = col,
                              )

  if (use_slanter_ordering) {
    if (!is.null(row_split)) {
      row_order <- unlist(lapply(split(1:nrow(mat), in_community_var$community), function(idx) {
        if (length(idx) == 1) return(index)
        similarity <- cor(t(mat)[idx, idx, drop = FALSE], method =  'spearman')
        similarity[similarity < 0] <- 0
        row_order <- slanter::slanted_orders(similarity)$rows
        idx[row_order]
      }), use.names = FALSE)
    } else {
      similarity <- cor(mat %>% t, method =  'spearman')
      similarity[similarity < 0] <- 0
      row_order <- slanter::slanted_orders(similarity)$rows %>% unname
    }
    cluster_rows <- FALSE
  }

  ComplexHeatmap::Heatmap(mat,
                          cluster_columns = FALSE,
                          clustering_distance_rows='spearman',
                          cluster_rows=cluster_rows,
                          column_split = column_split,
                          show_row_dend = FALSE,
                          name = name,
                          column_title = title,
                          column_gap = unit(4, "mm"),
                          col = col,
                          left_annotation = lfc_annotation,
                          row_order = row_order,
                          row_split = row_split,
                          ...
  )
}

#' Plot the F matrix (gene-topic loadings) as a heatmap
#'
#' @description Creates a heatmap visualization of the fastTopics F matrix showing
#'              gene loadings for each topic/program.
#'
#' @param F Data frame with gene column and topic score columns
#' @param ordered_des Data frame with gene ordering from differential expression (columns: gene, gene_score, topic)
#' @param marked_genes Character vector of genes to mark with arrows (default: empty)
#' @param col Color palette (default: cyan2orange(10))
#' @param column_title Title for columns (default: "Programs")
#' @param name Legend name (default: "Scaled program score")
#' @param use_raster If TRUE, rasterize heatmap cells (default: TRUE)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return ComplexHeatmap object with optional gene markers
#'
#' @note Genes are ordered by their top topic assignment from ordered_des.
plot_F_matrix <- function(F,
                          ordered_des,
                          marked_genes = c(),
                          col=cyan2orange(10),
                          column_title = "Programs",
                          name='Scaled program score',
                          use_raster=TRUE,
                          ...) {

  ordered_des <- ordered_des %>%
    group_by(gene) %>%
    slice_max(order_by = gene_score, n = 1) %>%
    ungroup %>%
    arrange(as.integer(gsub("k", "", topic)))

  df <- F %>%
    as.data.frame() %>%
    column_to_rownames("gene") %>%
    filter(rownames(.) %in% ordered_des$gene) %>%
    arrange(match(rownames(.), ordered_des$gene)) %>%
    t %>% scale %>% t

  hm <- Heatmap(df,
          show_row_dend = FALSE,
          show_row_names = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col=col,
          name=name,
          row_title = "Genes",
          column_title = column_title,
          use_raster = use_raster,
          raster_quality=5,
          ...)

  if (length(marked_genes) > 0) {
    hm <- hm + rowAnnotation(link = anno_mark(at = which(rownames(df) %in% marked_genes),
                             labels = rownames(df)[rownames(df) %in% marked_genes],
                             labels_gp = gpar(fontsize = 6), padding = unit(1, "mm")))
  }

  hm
}

#' Fit dynamics and create gene heatmap in one step
#'
#' @description Convenience function that fits dynamics model and creates
#'              a gene dynamics heatmap in a single call.
#'
#' @param data AnnData-like object with trajectory information in $uns
#' @param pseudobulks Named list of pseudobulk matrices by cell type take from "pseudobulks.R" function load.pseudobulks.of.projid()
#' @param ct Cell type name ("cux2p" -> "cux2+", "cux2m" -> "cux2-")
#' @param genes Character vector of genes to plot
#' @param title Plot title
#' @param both_trajectories Show both trajectories (default: FALSE)
#' @param scale Scale gene values (default: TRUE)
#' @param ... Additional arguments passed to plot.genes.dynamics.heatmap()
#'
#' @return ComplexHeatmap object
#'
#' @note Requires fit.dynamics() function to be available.
fit.plot.genes.dynamics.heatmap <- function(data, pseudobulks, ct, genes, title,
                                            both_trajectories = FALSE, scale = TRUE,
                                            ...) {
  if (ct == 'cux2p') {
    ct <- 'cux2+'
  } else if (ct == 'cux2m') {
    ct <- 'cux2-'
  }

  features <- pseudobulks[[ct]][, genes, drop=FALSE] %>%
    as.data.frame() %>%
    filter(rownames(.) %in% rownames(data)) %>%
    arrange(match(rownames(.), rownames(data)))

  dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                           features = features,
                           trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                           trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                           evaluate.fit = T,
                           bootstrap = F)

  plot.genes.dynamics.heatmap(dynamics, genes=genes,
                              both_trajectories=both_trajectories,
                              title, scale = scale, ...)
}

#' Heatmap of trait associations across multiple groups
#'
#' @description Creates a grouped heatmap showing trait associations with optional
#'              row grouping and custom labels. Similar to plot.trait.associations.heatmap
#'              but with additional row grouping support.
#'
#' @param df Data frame with trait association results
#' @param params Character vector of traits to include (default: amyloid, tangles, slope)
#' @param row.by Column for row identifiers (default: "covariate")
#' @param col.by Column for column identifiers (default: "trait")
#' @param value.by Column for heatmap values (default: "tstat")
#' @param pval.by Column for p-values (default: "adj.pval")
#' @param cols Color palette function (default: window2fox)
#' @param show.only.significant Filter to significant rows (default: TRUE)
#' @param plot If TRUE, draw heatmap (default: TRUE)
#' @param show_row_dend Show row dendrogram (default: FALSE)
#' @param row_order Optional row order vector
#' @param column_labels Named vector for column label replacement
#' @param row_labels Named vector for row label replacement
#' @param row_group_mapping Named vector mapping row labels to groups for splitting
#' @param clustering_distance_rows Distance metric (default: "euclidean")
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return If plot=TRUE, draws and returns invisibly. If FALSE, returns list(hm, legend).
plot.trait.assocation.multiple.groups <- function(df,
                                                  params=c("sqrt_amyloid_mf","sqrt_tangles_mf","cogng_random_slope"),
                                                  row.by="covariate", col.by="trait", value.by="tstat", pval.by="adj.pval",
                                                  cols=window2fox,
                                                  show.only.significant=T, plot=T,
                                                  show_row_dend=F,
                                                  row_order=NULL,
                                                  column_labels=NULL,
                                                  row_labels=NULL,
                                                  row_group_mapping=NULL,
                                                  clustering_distance_rows="euclidean", ...) {

    if("pandas.core.frame.DataFrame" %in% class(df))
      df <- py_to_r(df)

    df <- df %>% filter(trait %in% params)

    if(show.only.significant)
      df <- df %>% group_by_at(c(row.by)) %>%
      filter(if_any(pval.by, ~sum(. < .05) > 0)) %>%
      ungroup() %>%
      data.frame(check.names=FALSE)

    df$sig <- cut(df[[pval.by]], c(-.1, .001, .01, .05, Inf), c("***", "**", "*", ""))

    if (nrow(df) == 0) {
      print(grid::grid.text("No trait associated"))
      return()
    }

    vals <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = all_of(value.by), values_fill = NA_real_) %>%
      tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params))
    sig <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = "sig", values_fill = "") %>%
      tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params))

    if (!is.null(row_order)) {
      vals <- vals %>% arrange(match(make.names(rownames(.)), make.names(row_order)))
      sig <- sig %>% arrange(match(make.names(rownames(.)), make.names(row_order)))
    }

    vals <- vals %>% as.matrix()
    sig <- sig %>% as.matrix()

    v <- max(abs(vals), na.rm = T)
    col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])

    if (!is.null(column_labels) & length(column_labels) != ncol(vals)) {
      column_labels <- column_labels[colnames(vals)]
    } else if (is.null(column_labels)) {
      column_labels <- colnames(vals)
    }

    if (!is.null(row_labels) & length(row_labels) != nrow(vals)) {
      row_labels <- row_labels[make.names(rownames(vals))]
    } else if (is.null(row_labels)) {
      row_labels <- rownames(vals)
    }
    if (!is.null(row_group_mapping)) {
      row_split <- row_group_mapping[row_labels]
    } else {
      row_split <- NULL
    }


    hm <- Heatmap(vals,
                  name = value.by,
                  col = circlize::colorRamp2(col.vals, cols(length(col.vals))),
                  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                  cluster_columns = F,
                  clustering_distance_rows = clustering_distance_rows,
                  show_row_dend = show_row_dend,
                  column_labels = column_labels,
                  row_labels = row_labels,
                  row_split = row_split,
                  ...)
    legend <- Legend(title = "adj.pval", pch = c("***","**","*"), type = "points", labels = c("<0.001","<0.01", "<0.05"))
    if(!plot)
      return(list(hm=hm, legend=legend))
    draw(hm, annotation_legend_list = list(legend), merge_legend=T)
}

#' Calculate quantile thresholds from healthy/start group
#'
#' @description Computes quantiles of topic scores in the "start" group and
#'              calculates what percentile those thresholds correspond to in the "prAD" group.
#'
#' @param df Data frame with 'group' column and topic score column
#' @param topic Name of the topic/program score column
#' @param q Numeric vector of quantiles to compute (default: 0.5, 0.8, 0.9, 0.95, 0.99)
#'
#' @return Data frame with columns:
#'         - v: quantile threshold value
#'         - prAD_quantile: what percentile that threshold is in prAD group
#'
#' @note Assumes df has 'group' column with values "start" and "prAD".
topic_quantiles_healthy <- function(df, topic, q = c(0.5,0.8,0.9,0.95, 0.99)) {
  q_start <- quantile(df %>% filter(group == 'start') %>% pull(topic), q) %>% data.frame(v=.)
  q_start$prAD_quantile <- sapply(q_start %>% pull(v), function (v) {
    below_threshold <- mean((df %>% filter(group == 'prAD') %>% pull(topic)) < v)

    return(as.integer(below_threshold * 100))
  })
  q_start
}

#' Plot density distributions with quantile threshold lines
#'
#' @description Creates overlapping density plots for "start" and "prAD" groups
#'              with vertical lines showing quantile thresholds.
#'
#' @param df Data frame with 'group' column and topic scores
#' @param topic Name of the topic/program score column
#' @param quantiles_df Data frame with columns: threshold, start_quantile, prAD_quantile
#' @param title Plot title
#' @param qunatiles_to_plot Character vector of quantiles to show (default: "80%", "90%", "95%")
#'
#' @return ggplot object
plot_topic_quantiles <- function(df, topic, quantiles_df, title, quantiles_to_plot=c("80%", "90%", "95%")) {
  quantiles_df <- quantiles_df %>% filter(start_quantile %in% quantiles_to_plot)
  ggplot(df, aes_string(x=topic)) +
    geom_density(data=df %>% filter(group == 'prAD'), mapping=aes(fill='prAD'), alpha=0.4) +
    geom_density(data=df %>% filter(group == 'start'), mapping=aes(fill='start'), alpha=0.4) +
    geom_vline(data=quantiles_df,
               mapping = aes(xintercept = threshold),
               linetype = "dashed", color='black') +
    geom_text(data=quantiles_df, aes(x=threshold, y= Inf, label = start_quantile), colour = "#619CFF", vjust = 3,  size = 5) +
    geom_text(data=quantiles_df, aes(x=threshold, y= Inf, label = glue::glue("{prAD_quantile}%")), colour = "#F8766D", vjust = 5, size = 5) +
    theme_minimal() +
    labs(title = title, x="program score", y="density") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

#' Map neuronal topic names to colors
#'
#' @description Returns a named vector of colors based on neuronal cell type
#'              indicated by topic names.
#'
#' @param topics Character vector of topic names (e.g., "Inh.k1", "Exc.Cux2P.k2")
#'
#' @return Named character vector of hex colors:
#'         - Inhibitory (Inh): orange (#ef6c00)
#'         - Excitatory Cux2+ (Exc.Cux2P): dark red (#7c242c)
#'         - Excitatory Cux2- (Exc.Cux2M): teal (#4db6ac)
neuronal_topic_to_color <- function(topics) {
  sapply(topics, function(topic) case_when(str_detect(topic, "Inh") | topic == "inhibitory"  ~ "#ef6c00",
                                           str_detect(topic, "Exc.Cux2P") | topic == "cux2p" ~ "#7c242c",
                                           str_detect(topic, "Exc.Cux2M") | topic == "cux2m" ~ "#4db6ac"),
         USE.NAMES = T)
}

#' Create violin plots comparing groups with statistical significance
#'
#' @description Creates faceted violin plots with boxplot overlays and t-test
#'              significance annotations comparing two groups.
#'
#' @param df Data frame with numeric columns to plot
#' @param group Vector of group assignments (length must match nrow(df))
#' @param ylab Y-axis label (default: "")
#' @param xlab X-axis label (default: "Group")
#' @param title Plot title (default: "Violin plot")
#' @param fill_col Legend title for fill color (default: "Group")
#' @param test Statistical test for ggsignif (default: "t.test")
#'
#' @return ggplot object with faceted violin plots
#'
#' @note Uses ggsignif::geom_signif() for statistical annotations.
#' @note Assumes exactly 2 groups (comparisons = list(c(1, 2))).
#'       Will fail silently or incorrectly for more groups.
plot_violin_by_group <- function(df, group, ylab = "", xlab = "Group", title = "Violin plot", fill_col = "Group", test = "t.test") {
  # Ensure group is the same length
  stopifnot(nrow(df) == length(group))
  
  df_long <- df %>%
    mutate(Group = group) %>%
    pivot_longer(cols = -Group, names_to = "variable", values_to = "value")
  
  p <- ggplot(df_long, aes(x = as.factor(Group), y = value, fill = as.factor(Group))) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.2) +
    facet_wrap(~variable, scales = "free_y") +
    geom_signif(comparisons = list(c(1, 2)), test = test) +
    theme_bw() +
    labs(x = xlab, y = ylab, title = paste(title, " (n=", sum(!is.na(group)), ")"), fill = fill_col) +
    guides()
  
  return(p)
}


plot_violin_signature_per_group <- function(df, feature_name, ct) {
  legend_df <- data.frame(signif = c("*", "**", "***"), x = 1, y = 1)
      
  comparisons_list <- list(c("start neither", "prAD cell-cycle reentry"),
                            c("start neither", "prAD stress/apoptosis"),
                            c("prAD cell-cycle reentry", "prAD stress/apoptosis"))
  p_values <- sapply(comparisons_list, function(comp) {
    test_result <- wilcox.test(formula=as.formula(glue::glue("{feature_name} ~ cell_group")),
                                data = df %>% filter(cell_group %in% comp))
    test_result$p.value
  })
  
  p_adj <- p.adjust(p_values, method = "BH")
  
  
  signif_labels <- ifelse(p_adj < 0.001, "***",
                          ifelse(p_adj < 0.01, "**",
                                  ifelse(p_adj < 0.05, "*", "ns")))
  
  annotation_df <- data.frame(
    group1 = sapply(comparisons_list, `[[`, 1),
    group2 = sapply(comparisons_list, `[[`, 2),
    p_adj = p_adj,
    label = signif_labels
  )
  
  means_per_group <- df %>%
      group_by(cell_group) %>%
      summarise(mean_val = mean(.data[[feature_name]], na.rm=TRUE))

  ggviolin(df, x="cell_group", fill="cell_group", y=feature_name, legend.title='cell group') + 
      stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") + 
      geom_signif(
        comparisons = comparisons_list,
        annotations = annotation_df$label,
        y_position = max(df[[feature_name]]) + c(0.1, 0.25, 0.40),
        tip_length = 0.02
      ) +
      geom_point(data = legend_df, aes(x = x, y = y, shape = signif), color = "black", size = 4) +
      scale_shape_manual(
        name = "FDR",
        values = c("FDR" = 8),  # Use star shape for legend
        labels = c(glue::glue("* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001"))
      ) +
      geom_text(
        data = cell_group_counts,
        aes(x = cell_group, y = max(df[[feature_name]]) + 0.7, label = glue::glue("n={n}")),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = means_per_group,
        aes(x = cell_group, y = max(df[[feature_name]]) + 0.55, label = glue::glue("{round(mean_val, 2)}")),
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values=c('start neither'='#BCBEC0',
                                  'prAD neither'='#6D6E71',
                                  'prAD cell-cycle reentry'='#4db6ac',
                                  'prAD stress/apoptosis'="#ef6c00")) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme(legend.position = "bottom", legend.box = "vertical") +
      labs(y=feature_name) +
      ggtitle(glue::glue("{feature_name} in prAD vs start neither cells  - {ct}"))
}
