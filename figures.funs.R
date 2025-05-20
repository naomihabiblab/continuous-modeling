
#####################################################################################################################
#                                           General settings and functions                                          #
#####################################################################################################################
invisible(lapply(c("reshape2","ggplot2","cowplot","dplyr","tibble","ComplexHeatmap","reticulate","anndata",
                   "patchwork","colorspace","ggrepel","rhdf5", "ggnewscale", "gg3D"), library, character.only = T))

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
window2fox <- colorRampPalette(c("#008585", "white", "#c7522a"))

white2orange <- colorRampPalette(c("white", "#ef6c00"))

clay_and_sea_pallet <- colorRampPalette(c(
  "#74a892", "#b8cdab", "#e5c185","#e0a278", "#db836b", "#c7522a", "#ad4623"
), bias=0.4)


#####################################################################################################################
#                                                 Plotting Functions                                                #
#####################################################################################################################

#'
#'
#'
#'
#'
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


#' Scatter plot of trait associations
#' For each given parameter, scatter plot the trait association results with the scaled regression coefficient on the
#' x-axis, -log10 of the adjusted pvalue of the y-axis, colored and sized by the R^2
#'
#' @param trait.analysis.df dataframe of results, with the structure as returned from the associate.traits function
#' @param params vector of traits to plot for
#' @param label.thr Adjusted p-value threshold for adding labels
#' @param sig.thr Threshold of significance - used as threshold above which
#'
# plot.trait.associations <- function(trait.analysis.df, params, label.thr=0.05, sig.thr=0.01) {
#   associations <- trait.analysis.df %>%
#     filter(trait %in% params) %>%
#     mutate(trait = factor(trait, levels=params, ordered=T),
#            label = case_when(adj.pval < label.thr ~ covariate),
#            color = r.sq,
#            size = r.sq)
#
#   return(lapply(levels(associations$trait), function(t)
#     ggplot(associations %>% filter(trait == t), aes(beta, -log10(adj.pval), label=label,color=color,size=size)) +
#       geom_point() +
#       geom_text_repel(size=4, color="black", force = 3) +
#       geom_hline(yintercept = -log10(sig.thr), linetype="dashed") +
#       scale_y_sqrt() +
#       scale_color_gradient(na.value = "lightgrey", low="lightgrey", high = "red3") +
#       labs(x=NULL, y=NULL, title=NULL, color=NULL, size=NULL)+
#       labs(x="Effect size", title=t)) %>%
#       ggpubr::ggarrange(plotlist = ., common.legend = T, legend = "right", nrow = 1))
# }

plot.trait.associations <- function(trait.analysis.df, params, label.thr=0.05, sig.thr=0.01, nrow=1) {
  associations <- trait.analysis.df %>%
    filter(trait %in% params) %>%
    dplyr::mutate(trait = factor(trait, levels=params, ordered=T),
           label = case_when(adj.pval < label.thr ~ covariate),
           color = r.sq,
           size = r.sq)

  return(lapply(levels(associations$trait), function(t)
    ggplot(associations %>% filter(trait == t), aes(beta, -log10(adj.pval), label=label,color=color,size=size)) +
      geom_hline(yintercept = -log10(sig.thr), linetype="dashed",size=2) +
      geom_hline(yintercept = -log10(0.05), linetype="dashed",size=1, color="grey40") +
      geom_point() +
      ggrepel::geom_label_repel(size=4,face="bold", color="black", force = 3) +

      scale_y_sqrt() +
      scale_color_gradient(na.value = "lightgrey", low="lightgrey", high = "red3") +
      labs(x="Effect size", y=NULL, title=t, color=NULL, size=NULL)+
      theme(axis.line = element_line(size=1),
            axis.ticks = element_line(size=1),
            axis.text = element_text(face="bold"),
            text = element_text(face="bold"))) %>%
      ggpubr::ggarrange(plotlist = ., common.legend = T, legend = "right", nrow = nrow))
}

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
  legend <- Legend(title = "adj.pval", pch = c("****","***","**","*"), type = "points", labels = c("<0.0001","<0.001","<0.01", "<0.05"))
  if(!plot)
    return(list(hm=hm, legend=legend))
  draw(hm, annotation_legend_list = list(legend), merge_legend=T)
}

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



plot.landscape.3D <- function(
    features = NULL,
    embedding = "X_all_3d_phate",
    cols = awesome_pallet,
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

#'
#'
plot.dynamics.old <- function(fits, predicts, cols="black", fold.change=F, branch=NULL, facet.by.branch=T, facet.scales="free", ncol=1, label=T,
                          features = NULL) {
  if (!is.null(features)) {
    fits <- fits %>% filter(feature %in% features)
    predicts <- predicts %>% filter(feature %in% features)
  }

  if(fold.change) {
    fits <- fits %>% group_by(feature, trajectory) %>% mutate(fit = fit/fit[x == min(x)])
    predicts <- predicts %>% group_by(feature, trajectory) %>% mutate(fit = fit/fit[x == min(x)])
  }

  if(!is.null(branch)) {
    fits <- fits[fits$trajectory == branch,]
    predicts <- predicts[predicts$trajectory == branch,]
  }

  group.by = ifelse(facet.by.branch, "feature", "trajectory")
  facet.by = ifelse(facet.by.branch, "~trajectory", "~feature")

  if(is.character(cols) & length(cols) == 1)
    cols <- rep(cols, length(unique(fits[, group.by])))
  
  names <- unique(fits[, group.by])

  tryCatch({
    names <- unique(fits[, group.by])[[group.by]]
  }, error=function(e){} )

  cols <- setNames(cols, names)

  p <- ggplot(fits, aes_string(x="x", color=group.by)) +
    geom_line(aes(y=fit), size=.5) +
    geom_ribbon(aes_string("x", "fit", ymin="fit-se.fit", ymax="fit+se.fit", fill=group.by), data=predicts, alpha=.1, linetype="dashed") +
    geom_line(aes(x, fit), data=predicts, linetype="dashed", size=.5) +
    facet_wrap(facet.by, scales=facet.scales, ncol = ncol) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x=NULL, y=NULL, color=NULL) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = "none")

  if(label)
    p <- p + ggrepel::geom_label_repel(aes_string("x", "fit", label=group.by), fits %>% group_by(feature, trajectory) %>% slice_max(x, n=1))
  return(p)
}

plot.dynamics.wrapper <- function(dynamics, features, trajectories = dynamics$fitted.vals %>% pull(trajectory) %>% unique(), ...) {
  if("pandas.core.frame.DataFrame" %in% class(dynamics$fitted.vals))
    dynamics$fitted.vals <- py_to_r(dynamics$fitted.vals)

  if("pandas.core.frame.DataFrame" %in% class(dynamics$pred.vals))
    dynamics$pred.vals <- py_to_r(dynamics$pred.vals)
  
  stopifnot(length(features) > 0)
  
  fits <- dynamics$fitted.vals %>% filter(feature %in% features) %>% filter(trajectory %in% trajectories) %>% dplyr::arrange(match(feature, features)) 
  pred <- dynamics$pred.vals %>% filter(feature %in% features) %>% filter(trajectory %in% trajectories) %>% dplyr::arrange(match(feature, features))
  
  stopifnot(nrow(fits) > 0)
  warningCondition((fits$feature %>% n_distinct) != length(features))
  
  return(plot.dynamics(fits, pred, ...))
}


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


plot.dynamics_topic_cluster.wrapper <- function(dynamics, features, cols = NULL, facet.by = c("trajectory","feature"),
                          overlap.pseudotime=NULL, include.points=T, min.point.alpha=.1, label=T,
                          legend.position = "right", raster=F, ...) {
  if("pandas.core.frame.DataFrame" %in% class(dynamics$fitted.vals))
    dynamics$fitted.vals <- py_to_r(dynamics$fitted.vals)

  if("pandas.core.frame.DataFrame" %in% class(dynamics$pred.vals))
    dynamics$pred.vals <- py_to_r(dynamics$pred.vals)

  fits <- dynamics$fitted.vals %>% filter(feature %in% features) # %>% mutate(feature = gsub("\\d+\\.k", "T", feature))
  preds <- dynamics$pred.vals %>% filter(feature %in% features) # %>% mutate(feature = gsub("\\d+\\.k", "T", feature))
  # names(cols) <- gsub("\\d+\\.k", "T", names(cols))
  
  facet.by = match.arg(facet.by)
  group.by = setdiff(c("trajectory","feature"), facet.by)

  groups <- unique(fits[[group.by]])
  # groups <- gsub("\\d+\\.k", "T", features)
  if(is.null(cols)) {
    cols <- setNames(scales::hue_pal()(length(groups)), groups)
  }

  if (length(groups) < 2) {
    stop("Call with one topic and one cluster feature")
  }

  first_group <- groups[1:(length(groups)-1)]
  second_group <- groups[[length(groups)]]

  preds_first_group <- preds %>% filter(feature%in%first_group)
  preds_second_group <- preds %>% filter(feature == second_group)

  ymax <- max(preds_first_group$fit+2*preds_first_group$se.fit, na.rm = T)
  ymin <- min(preds_first_group$fit-2*preds_first_group$se.fit, na.rm = T)

  second_ymax <- max(preds_second_group$fit+2*preds_second_group$se.fit, na.rm = T)
  second_ymin <- min(preds_second_group$fit-2*preds_second_group$se.fit, na.rm = T)
  
  a.diff <- ymax - ymin
  b.diff <- second_ymax - second_ymin
  
  preds[preds[, group.by] == second_group,]$fit <- (preds[preds[, group.by] == second_group,]$fit - second_ymin) / b.diff * a.diff + ymin
  preds[preds[, group.by] == second_group,]$se.fit <- (preds[preds[, group.by] == second_group,]$se.fit ) / b.diff * a.diff
  
  
  ymax <- max(preds$fit+2*preds$se.fit, na.rm = T)
  ymin <- min(preds$fit-2*preds$se.fit, na.rm = T)
  
  p <- ggplot(fits) +
    facet_wrap(paste0("~",facet.by), ...) +
    scale_x_continuous(expand = c(0,0), labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ymax), labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x)),
                       sec.axis = sec_axis(~((. -ymin) * b.diff / a.diff) + second_ymin,)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = legend.position,
          axis.line = element_line())

  if(!is.null(overlap.pseudotime)) {
    p <- p +
      geom_ribbon(aes(x,ymax=ymax,ymin=ymin),
                  data.frame(x=c(0, overlap.pseudotime), ymax=rep(ymax, 2), ymin=rep(ymin,2)),
                  fill="black", alpha=.05) +
      geom_vline(xintercept = overlap.pseudotime, linetype="dashed")
  }

  for(g in groups) {
    p <- p + geom_ribbon(aes_string("x", "fit", ymin="fit-2*se.fit", ymax="fit+2*se.fit",
                                    fill=group.by),
                         preds[preds[, group.by] == g,], alpha=.2, color=NA, inherit.aes = F) +
      geom_line(aes_string("x", "fit", color=group.by),
                preds[preds[, group.by] == g,], size=.5, inherit.aes = F) +
      scale_color_manual(values = cols[[g]]) +
      scale_fill_manual(values = cols[[g]]) +

      new_scale("color") +
      new_scale("fill") +
      guides(color_new=guide_legend(g), fill_new=guide_legend(g))
  }

  if(label)
    p <- p + ggrepel::geom_label_repel(aes_string("x", "fit", label=group.by, color=group.by),
                                       preds %>% group_by(feature, trajectory) %>% slice_max(x, n=1),
                                       min.segment.length = unit(0,"pt"), show.legend = F) +
    scale_color_manual(values = cols)
  return(p)
}



#'
#'
#'
#'
#'
#'
#'
#'
#'
plot.dynamics.clustered <- function(dynamics,
                                    bottom.annotation.df,
                                    cols = colorspace::sequential_hcl(n = 20, palette = "Reds", rev = T),
                                    #bot.ann.cols = function(vals) circlize::colorRamp2(quantile(vals, probs = seq(0,1, length.out=20)), cols),
                                    bot.ann.cols = function(vals) circlize::colorRamp2(seq(min(vals), max(vals), length.out=20), cols),
                                    branch.order = NULL,
                                    mark.rows = NULL,
                                    ...) {

  # Bottom annotations of traits' dynamics
  bot.ann <- bottom.annotation.df %>% reshape2::dcast(param~branch+x, value.var = "fit") %>%
    column_to_rownames("param") %>%
    t() %>% data.frame()

  if(class(bot.ann.cols) == "function") {
    bot.ann.cols = sapply(colnames(bot.ann), function(c) bot.ann.cols(bot.ann[,c]), simplify = F, USE.NAMES = T)
  }

  bot.ann <-  HeatmapAnnotation(df=bot.ann,
                                col = bot.ann.cols,
                                annotation_name_side = "left",
                                annotation_legend_param = list(direction="horizontal"))

  # Top annotations - pseudotime in branches
  top.ann <- data.frame(name=colnames(dynamics)) %>%
    tidyr::separate(name, c("branch","ps"), "_", remove = F) %>%
    group_by(branch) %>%
    dplyr::mutate(ps = as.numeric(ps),
           ps.text = case_when(row_number() %in% c(1,n()) ~ as.character(round(ps,1)))) %>%
    column_to_rownames("name")

  top.ann <- HeatmapAnnotation(a=anno_simple(top.ann$ps, pch = top.ann$ps.text,
                                             col = circlize::colorRamp2(seq(0, 1, length.out=20), cols)),
                               annotation_name_side = "left",
                               annotation_label = c("Pseudotime"))

  column.split = gsub("_.*", "", colnames(dynamics))
  if(!is.null(branch.order))
    column.split = factor(column.split, levels = branch.order)

  right_annotation <- NULL
  if(!is.null(mark.rows)) {
    mark.rows <- intersect(mark.rows, rownames(dynamics))
    ids <- which(rownames(dynamics) %in% mark.rows)
    ann <- anno_mark(at = ids, labels = rownames(dynamics)[ids], which = "row")
    right_annotation <- rowAnnotation(named_rows=ann, gp=gpar(fontsize=10, fontface="bold"))
  }

  dynamics.heatmap <- Heatmap(
    dynamics,
    cluster_columns = F,
    show_column_names = F,

    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize=8),

    column_split = column.split,
    cluster_column_slices = F,

    #col = colorRampPalette(c("blue","white","red"))(20),#colorspace::diverge_hcl(20, palette = "Blue-Red 3"),
    top_annotation = top.ann,
    bottom_annotation = bot.ann,
    right_annotation = right_annotation,

    border = T,
    heatmap_legend_param = list(title="scaled dynamics", direction = "horizontal", legend_width = unit(80, "pt")),
    ...)

  return(prepare(dynamics.heatmap))
}


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
    vals <<- vals
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

topic_quantiles_healthy <- function(df, topic, q = c(0.5,0.8,0.9,0.95, 0.99)) {
  q_start <- quantile(df %>% filter(group == 'start') %>% pull(topic), q) %>% data.frame(v=.)
  q_start$prAD_quantile <- sapply(q_start %>% pull(v), function (v) {
    below_threshold <- mean((df %>% filter(group == 'prAD') %>% pull(topic)) < v)

    return(as.integer(below_threshold * 100))
  })
  q_start
}

plot_topic_quantiles <- function(df, topic, quantiles_df, title, qunatiles_to_plot=c("80%", "90%", "95%")) {
  quantiles_df <- quantiles_df %>% filter(start_quantile %in% qunatiles_to_plot)
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

neuronal_topic_to_color <- function(topics) {
  sapply(topics, function(topic) case_when(str_detect(topic, "Inh") | topic == "inhibitory"  ~ "#ef6c00",
                                           str_detect(topic, "Exc.Cux2P") | topic == "cux2p" ~ "#7c242c",
                                           str_detect(topic, "Exc.Cux2M") | topic == "cux2m" ~ "#4db6ac"),
         USE.NAMES = T)
}