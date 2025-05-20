source("utils.R")


SEA_AD_ROOT_PATH <- file.path(get_shared_dir(), "SeaAD/sc_data_files/sea_ad_single_cell_profiling")
SEA_AD_SUPPLEMENTARY_PATH <- file.path(SEA_AD_ROOT_PATH, "MTG/RNAseq/Supplementary Information")


sea_ad_read_de <- function(subset=c('all', 'early', 'late')) {
  subset <- match.arg(subset)

  if (subset == 'all') {
    subset <- ''
  } else {
    subset <- paste0('_', subset)
  }

  effect_size <- anndata::read_h5ad(file.path(SEA_AD_SUPPLEMENTARY_PATH, "Nebula Results", glue::glue("effect_sizes{subset}.h5ad")))
  pvalue <- anndata::read_h5ad(file.path(SEA_AD_SUPPLEMENTARY_PATH, "Nebula Results", glue::glue("pvalues{subset}.h5ad")))

  return(list(effect_size=effect_size, pvalue=pvalue))
}

sea_ad_read_de <- memoize_first(sea_ad_read_de)

anndata_get_var <- function(adata, keys, layer=NULL) {
  if (is.null(layer)) {
    x <- adata$X
  } else {
    x <- adata$layers[[layer]]
  }

  keys <- unique(keys)

  obs_keys <- keys[keys %in% rownames(x)]
  var_keys <- keys[keys %in% colnames(adata$var)]
  missing_keys <- setdiff(keys, union(obs_keys, var_keys))
  if (length(missing_keys) > 0) {
    stop(paste0("Some keys are missing: ", paste(missing_keys, collapse = ", ")))
  }

  df <- x[obs_keys, , drop=FALSE] %>% t

  df <- cbind(df, adata$var[, var_keys, drop=FALSE])

  df
}


sea_ad_de_summary <- function(stage, genes, p_value_threshold = 0.05) {
  ret <- sea_ad_read_de(stage)
  pvalue <- ret$pvalue
  effect_size <- ret$effect_size

  effect_df <- anndata_get_var(effect_size, c(genes, "Class"), layer='effect_sizes')
  p_value_df <- anndata_get_var(pvalue, genes, layer='pvalues')

  significant_df <- (effect_df %>% dplyr::select(-Class)) * (abs(p_value_df) < p_value_threshold)
  significant_df <- cbind(significant_df, effect_df %>% dplyr::select(Class))
  df <- significant_df %>%
    rownames_to_column('Subtype') %>%
    tidyr::pivot_longer(-c(Class, Subtype), names_to='gene', values_to='effect_size')


  res <- df %>%
    group_by(Class, gene) %>%
    summarise(mean=mean(effect_size),
              positive=sum(effect_size > 0),
              non_sig= sum(effect_size == 0),
              negative= sum(effect_size < 0)
    ) %>%
    mutate(diff = positive - negative) %>%
    filter(Class != 'Non-neuronal and non-neural') %>%
    mutate(stage=stage)

  res
}

sea_ad_de_plot_heatmap <- function(df,
                                   show_numbers=TRUE,
                                   height = 2*unit(15, "mm"),
                                   column_title = glue::glue("Number of subtype with de"),
                                   ...) {
  max_diff <- max(abs(df$diff))
  cols <- circlize::colorRamp2(c(-max_diff,0, max_diff/ 3, max_diff), c("#4db6ac","white","#ef6c00", "#c7522a"))

  stage <- df$stage[[1]]

  mat <- df %>%
    tidyr::pivot_wider(id_cols = "Class", names_from = "gene", values_from = diff) %>%
    mutate(Class=gsub("Neuronal: ", "", Class)) %>%
    column_to_rownames('Class') %>% t

  if (show_numbers) {
    cell_fun <- function(j, i, x, y, w, h, fill) grid.text(mat[i,j], x,y, gp = gpar(fontsize=8))
  } else {
    cell_fun <- NULL
  }

  Heatmap(mat,
          column_title = column_title,
          row_title = stage,
          show_column_dend = FALSE,
          cluster_rows = FALSE,
          col = cols,
          name = "Number of subtypes",
          cell_fun = cell_fun,
          ...
  )
}
