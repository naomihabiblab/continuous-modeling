#### Spatial script for paper
# Spatial transcriptomics analysis of AD tissue samples
# Analyzes co-localization of neuronal programs, astrocyte stress, and microglial activation

# =============================================================================
# PATH CONFIGURATION
# =============================================================================
# This section defines all paths used in the spatial analysis.
# Paths automatically detect if running on cluster or local mount.
#
# Required data:
# - Spatial h5ad object: {SPATIAL_BASE_DIR}/for_paper/final_obj/spatial_obj_for_paper.h5ad
# - Enrichment statistics:
#   - Programs: {SPATIAL_BASE_DIR}/for_paper/enrichment_data/programs/all_pvals_programs.csv
#   - Pathways: {SPATIAL_BASE_DIR}/for_paper/enrichment_data/pathways/all_pvals_pathways.csv
#
# Outputs:
# - Plots saved to: {PLOTS_DIR}/spatial/
# =============================================================================

# Libraries
library(reticulate)
library(ggplot2)
library(dplyr)
library(ggforce)
library(tidyverse)
library(ggrastr)

# Load project utilities
source("R/utils.R")  # Load get_result_dir() and get_shared_dir()
source("R/source_data.R")  # Source Data helpers (sd_add, sd_tee); writes to {SD_DIR}/staging only

# Detect environment and set base directories
PLOTS_DIR <- get_result_dir()           # Auto-detects cluster vs local
SHARED_DIR <- get_shared_dir()          # Auto-detects cluster vs local


SPATIAL_BASE_DIR <- file.path(get_result_dir(), "spatial")

# Specific input files
SPATIAL_OBJECT_PATH <- file.path(SPATIAL_BASE_DIR, "for_paper/final_obj/spatial_obj_for_paper.h5ad")
ENRICHMENT_BASE_DIR <- file.path(SPATIAL_BASE_DIR, "for_paper/enrichment_data")
ENRICHMENT_PROGRAMS_PATH <- file.path(ENRICHMENT_BASE_DIR, "programs/all_pvals_programs.csv")
ENRICHMENT_PROGRAMS__HEALTHY_PATH <- file.path(ENRICHMENT_BASE_DIR, "programs/all_pvals_programs_healthy.csv")
ENRICHMENT_PATHWAYS_PATH <- file.path(ENRICHMENT_BASE_DIR, "pathways/all_pvals_pathways.csv")

# Output directory for spatial plots
# SPATIAL_PLOTS_DIR env var lets a re-run write elsewhere (Source Data run: source_data/figures_rerun/spatial)
SPATIAL_PLOTS_DIR <- Sys.getenv("SPATIAL_PLOTS_DIR", file.path(PLOTS_DIR, "spatial"))
dir.create(SPATIAL_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Verify required files exist
if (!file.exists(SPATIAL_OBJECT_PATH)) {
  stop("Spatial object not found at: ", SPATIAL_OBJECT_PATH)
}
if (!file.exists(ENRICHMENT_PROGRAMS_PATH)) {
  warning("Programs enrichment data not found at: ", ENRICHMENT_PROGRAMS_PATH)
}
if (!file.exists(ENRICHMENT_PATHWAYS_PATH)) {
  warning("Pathways enrichment data not found at: ", ENRICHMENT_PATHWAYS_PATH)
}


#data:
anndata <- import("anndata", convert = FALSE)
adata <- anndata$read_h5ad(SPATIAL_OBJECT_PATH, backed = "r") 

# TODO: use obs and obsm and not the full adata?
h5ad_obs_to_df <- function(obs) {
  n <- length(obs$`_index`)
  df <- lapply(obs, function(field) {
    if (is.list(field) && !is.null(field$categories) && !is.null(field$codes)) {
      if (length(field$codes) == n) {
        field$categories[field$codes + 1L]
      } else {
        NULL
      }
    } else if (length(field) == n) {
      field
    } else {
      NULL
    }
  })
  df <- df[!sapply(df, is.null)]
  df <- as.data.frame(df)
  rownames(df) <- obs$`_index`
  df$`_index` <- NULL
  df
}

test_obs <- rhdf5::h5read(SPATIAL_OBJECT_PATH, "obs")
obs_df <- h5ad_obs_to_df(test_obs)

obs <- adata$obs
obsm <- adata$obsm

obs_df <- py_to_r(obs)

AD_samples <- obs_df%>%dplyr::select(c('sample_id','levels'))%>%unique()%>%filter(levels!='control')
AD_sample_id <- AD_samples$sample_id

non_AD_sample <- obs_df %>% dplyr::select(c('sample_id', 'levels')) %>% unique() %>% filter(levels == 'control')
non_AD_sample_id <- non_AD_sample$sample_id

#####################
#boxplots programs: #
#####################

plot_pairs_boxplot <- function(
    programs_data,
    pairs_df,
    y_offset = 0.15
) {
  
  sample_cols <- c(
    "A02092E1",
    "B01809C2",
    "B02008C6",
    "B02008D2",
    "B02009F6",
    "C02248B5"
  )
  
  sample_colors <- c(
    "A02092E1 (severe)"   = "#9E0B0F",  
    "B02008C6 (severe)"   = "#E31A1C",  
    
    "B02009F6 (moderate)" = "#FF8C1A",  
    "B02008D2 (moderate)" = "#FFC04D",  
    
    "C02248B5 (advanced)" = "#6B2E1F",  
    "B01809C2 (advanced)" = "#A0522D"  
  )
  sample_labels <- c(
    "A02092E1" = "A02092E1 (severe)",
    "B02008C6" = "B02008C6 (severe)",
    "B02009F6" = "B02009F6 (moderate)",
    "B02008D2" = "B02008D2 (moderate)",
    "C02248B5" = "C02248B5 (advanced)",
    "B01809C2" = "B01809C2 (advanced)"
  )
  
  ## keep only relevant pairs
  df_sub <- programs_data %>%
    semi_join(pairs_df, by = c("ref_topic", "neighbor_topic"))
  
  ## long format
  plot_long <- df_sub %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(
      pair_label = paste(ref_topic, "→", neighbor_topic),
      legend_label  = unname(sample_labels[as.character(sample)])
    )
  
  ## significance stars
  stars_df <- df_sub %>%
    mutate(
      pair_label = paste(ref_topic, "→", neighbor_topic),
      star = case_when(
        p_bh < 0.001 ~ "***",
        p_bh < 0.01  ~ "**",
        p_bh < 0.05  ~ "*",
        TRUE         ~ ""
      )
    ) %>%
    filter(star != "") %>%
    distinct(pair_label, star) %>%
    left_join(
      plot_long %>%
        group_by(pair_label) %>%
        summarise(y = max(value, na.rm = TRUE) + y_offset),
      by = "pair_label"
    )
  
  ## plot
  ggplot(plot_long, aes(x = pair_label, y = value, color = legend_label)) + 
    geom_boxplot(
      fill = NA,
      color = "black",
      outlier.shape = NA,
      width = 0.5
    ) +
    geom_jitter(width = 0.25, size = 3) +
    geom_text(
      data = stars_df,
      aes(x = pair_label, y = y, label = star),
      inherit.aes = FALSE,
      size = 6
    ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(values = sample_colors) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Reference → Neighbor Topic",
      y = "Mean Enrichment",
      color = "Sample"
    )
}

plot_pairs_boxplot_healthy <- function(
  programs_data_healthy,
  pairs_df,
  healthy_sample_cols,
  y_offset = 0.15
) {

  ## keep only relevant pairs
  df_sub <- programs_data_healthy %>%
    semi_join(pairs_df, by = c("ref_topic", "neighbor_topic"))

  ## derive ordered pair labels from pairs_df
  ordered_labels <- paste(pairs_df$ref_topic, "→", pairs_df$neighbor_topic)

  ## long format
  plot_long <- df_sub %>%
    pivot_longer(
      cols = all_of(healthy_sample_cols),
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(
      pair_label = factor(
        paste(ref_topic, "→", neighbor_topic),
        levels = ordered_labels
      )
    )

  ## significance stars
  stars_df <- df_sub %>%
    mutate(
      pair_label = factor(
        paste(ref_topic, "→", neighbor_topic),
        levels = ordered_labels
      ),
      star = case_when(
        p_bh < 0.001 ~ "***",
        p_bh < 0.01  ~ "**",
        p_bh < 0.05  ~ "*",
        TRUE         ~ ""
      )
    ) %>%
    filter(star != "") %>%
    distinct(pair_label, star) %>%
    left_join(
      plot_long %>%
        group_by(pair_label) %>%
        summarise(y = max(value, na.rm = TRUE) + y_offset),
      by = "pair_label"
    )

  ## plot
  ggplot(plot_long, aes(x = pair_label, y = value)) +
    geom_boxplot(
      fill = NA,
      color = "black",
      outlier.shape = NA,
      width = 0.5
    ) +
    geom_jitter(width = 0.25, size = 1) +
    geom_text(
      data = stars_df,
      aes(x = pair_label, y = y, label = star, size = 2),
      inherit.aes = FALSE,
      size = 6
    ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Reference → Neighbor Topic",
      y = "Mean Enrichment",
      color = "Sample"
    )
}

# Source Data for the co-localization boxplots: per-sample ratios and the stats behind the stars
sd_pairs <- function(data, pairs_df, sample_cols, fig, panel, what) {
  df <- data %>% semi_join(pairs_df, by = c("ref_topic", "neighbor_topic"))
  sd_add(df %>% dplyr::select(ref_topic, neighbor_topic, bin, all_of(sample_cols)), fig, panel, "per_sample",
         description = glue::glue("{what}: co-localization enrichment of the neighbour program around the reference program, one column per sample (1 = no enrichment; boxplot points)"))
  sd_add(df %>% dplyr::select(ref_topic, neighbor_topic, bin, mean_enrichment, p_raw, shapiro_W, shapiro_p, p_bh), fig, panel, "stats",
         description = glue::glue("{what}: one-sample t-test of the per-sample values against 1 (p_raw) and BH-adjusted p (p_bh, gives the stars; BH within each reference program and neighbour cell type, 2-7 tests, as in spatial_analysis.ipynb), from the stored enrichment CSV"))
}
AD_SAMPLE_COLS <- c("A02092E1", "B01809C2", "B02008C6", "B02008D2", "B02009F6", "C02248B5")

programs_data <- read.csv(ENRICHMENT_PROGRAMS_PATH)
programs_data_healthy <- read.csv(ENRICHMENT_PROGRAMS__HEALTHY_PATH)

pairs_df_80 <- tibble(
  ref_topic = c(
    "AST_k6_score_bin80",
    "Ex_k3_score_bin80",
    "Microglia_k8_score_bin80",
    "AST_k6_score_bin80",
    "Ex_k9_score_bin80",
    "AST_k6_score_bin80"
  ),
  neighbor_topic = c(
    "Ex_k3_score_bin80",
    "AST_k6_score_bin80",
    "AST_k6_score_bin80",
    "Microglia_k8_score_bin80",
    "AST_k6_score_bin80",
    "Ex_k9_score_bin80"
  )
)

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "pairs_boxplot_80.pdf"), width=6, height=6, { 
  plot_pairs_boxplot(
    programs_data = programs_data,
    pairs_df = pairs_df_80
  ) %>% print
})
sd_pairs(programs_data, pairs_df_80, AD_SAMPLE_COLS, "Fig6", "f", "AD samples, 80% quantile programs")

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "pairs_boxplot_80_healthy.pdf"), width=6, height=6, { 
  plot_pairs_boxplot_healthy(
    programs_data = programs_data_healthy,
    pairs_df = pairs_df_80,
    healthy_sample_cols = non_AD_sample_id
  ) %>% print
})  
sd_pairs(programs_data_healthy, pairs_df_80, non_AD_sample_id, "ED10", "e", "Control samples, 80% quantile programs")

pairs_df_85 <- tibble(
  ref_topic = c(
    "AST_k6_score_bin85",
    "Ex_k3_score_bin85",
    "Microglia_k8_score_bin85",
    "AST_k6_score_bin85",
    "Ex_k9_score_bin85",
    "AST_k6_score_bin85"
  ),
  neighbor_topic = c(
    "Ex_k3_score_bin85",
    "AST_k6_score_bin85",
    "AST_k6_score_bin85",
    "Microglia_k8_score_bin85",
    "AST_k6_score_bin85",
    "Ex_k9_score_bin85"
  )
)

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "pairs_boxplot_85.pdf"), width=6, height=6, { 
  plot_pairs_boxplot(
    programs_data = programs_data,
    pairs_df = pairs_df_85
  ) %>% print
})
sd_pairs(programs_data, pairs_df_85, AD_SAMPLE_COLS, "ED10", "b", "AD samples, 85% quantile programs")


#pathways: 
pathways_data=read.csv(ENRICHMENT_PATHWAYS_PATH)

pairs_df_pathways <- tibble(
  ref_topic = c(
    "AST_k6_score_bin80",
    "AST_k6_score_bin80",
    "AST_k6_score_bin80",
    "Microglia_k8_score_bin80",
    "Ex_cell_cycle_score_second_pick_bimodal_threshold",
    "Ex_DNA_Repair_score_second_pick_bimodal_threshold"
  ),
  neighbor_topic = c(
    "Ex_Metallothioneins_bind_metals_score_second_pick_bimodal_threshold",
    "Ex_MT_genes_score_second_pick_bimodal_threshold",
    "Ex_response_to_oxidative_stress_score_second_pick_bimodal_threshold",
    "AST_k6_score_bin80",
    "AST_k6_score_bin80",
    "AST_k6_score_bin80"
  )
)

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "pairs_boxplot_pathways.pdf"), width=6, height=6, { 
  plot_pairs_boxplot(
      programs_data = pathways_data,
      pairs_df = pairs_df_pathways
    ) %>% print
})
sd_pairs(pathways_data, pairs_df_pathways, AD_SAMPLE_COLS, "ED10", "c", "AD samples, pathways (bimodal threshold) and programs (80% quantile)")




###########################
##### slices (Fig.4) ######
###########################

plot_topic_xy_all_samples <- function(df, topic_col,color_topic) {
  sample_ids <- unique(df$sample_id)
  plot_list <- list()
  
  for (sid in sample_ids) {
    df_s <- df %>% filter(sample_id == sid)
    
    df_s <- df_s %>%
      mutate(
        topic_value = .data[[topic_col]],
        is_na = is.na(topic_value)
      )
  
    is_binary <- all(na.omit(df_s$topic_value) %in% c(0,1))
    
    p <- ggplot(df_s, aes(x = x, y = y)) +
      geom_point_rast(data = df_s %>% filter(is_na), color = "grey80", size = 0.5, alpha = 0.8, raster.dpi = 300)
    
    if (is_binary) {
      
      p <- p +
        geom_point_rast(
          data = df_s %>% filter(!is_na & topic_value == 1),
          color = color_topic,
          size = 0.8,
          alpha = 0.8,
          raster.dpi = 300
        )
      
    } else {
      p <- p +
        geom_point_rast(
          data = df_s %>% filter(!is_na),
          aes(color = topic_value),
          size = 0.8,
          alpha = 0.9,
          raster.dpi = 300
        ) +
        scale_color_viridis_c(option = "plasma")
    }
    
    p <- p + theme_minimal() +
      theme(panel.grid = element_blank(), aspect.ratio = 1) +
      labs(title = paste0("Sample ", sid," topic=",topic_col), x = "x", y = "y", color = topic_col)
    
    plot_list[[sid]] <- p
  }
  
  return(plot_list)
}



plots <- plot_topic_xy_all_samples(obs_df%>%filter(sample_id%in%AD_sample_id), "Ex_k3_score_bin80",color_topic='#c7522a')
ggsave(
  file.path(SPATIAL_PLOTS_DIR, "full_slice_stress_topic.pdf"),
  plots[['B02009F6']],
  width = 6,
  height = 6
)

plots <- plot_topic_xy_all_samples(obs_df%>%filter(sample_id%in%AD_sample_id), "Ex_k9_score_bin80",color_topic = '#008585')
ggsave(
  file.path(SPATIAL_PLOTS_DIR, "full_slice_cell_cycle_topic.pdf"),
  plots[['B02009F6']],
  width = 6,
  height = 6
)
sd_add(obs_df %>% filter(sample_id == "B02009F6") %>% tibble::rownames_to_column("cell") %>%
         dplyr::select(cell, sample_id, annotation, x, y, Ex_k3_score_bin80, Ex_k9_score_bin80),
       "Fig4", "e", description = "Sample B02009F6, all cells. Ex_k3_score_bin80 = stress/apoptosis (left), Ex_k9_score_bin80 = cell-cycle reentry (right): 1 = coloured, 0 = excitatory cell not drawn, NA = non-excitatory cell drawn in grey")



##########################################
###  Fig 6, slices of co-localization ####
##########################################

plot_cell_neighborhood <- function(
    df,
    cell_index,
    col1,
    col2,
    col3,
    sample_id,
    col4 = NULL,          # NEW optional column
    square_size = 400,
    expand_factor = 0.15
) {
  
  target_cell <- df[cell_index, ]
  x0 <- target_cell$x
  y0 <- target_cell$y
  
  xmin <- x0 - square_size
  xmax <- x0 + square_size
  ymin <- y0 - square_size
  ymax <- y0 + square_size
  
  df_square <- df %>%
    dplyr::filter(x >= xmin & x <= xmax,
                  y >= ymin & y <= ymax) %>%
    dplyr::mutate(
      col1_val = factor(.data[[col1]], levels = c(0,1)),
      col2_val = factor(.data[[col2]], levels = c(0,1)),
      col3_val = factor(.data[[col3]], levels = c(0,1)),
      col4_val = if (!is.null(col4))
        factor(.data[[col4]], levels = c(0,1)) else NA
    )
  
  green_cells <- df_square %>% dplyr::filter(col3_val == 1)
  
  # Colors
  col1_color <- "#c7522a"
  col2_color <- "#7B6FD0"
  col3_color <- "#FFC107"
  col4_color <- "#008585"  
  
  p <- ggplot(df_square, aes(x, y)) +
    geom_point_rast(color = "grey80", size = 5, alpha = 0.8,raster.dpi = 300) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", size = 0.1) +
    ggforce::geom_circle(
      data = green_cells,
      aes(x0 = x, y0 = y, r = 200),
      color = col3_color,
      alpha = 0.1,
      size = 0.5,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_point_rast(data = df_square %>% dplyr::filter(col1_val == 1),
                    color = col1_color, size = 6, alpha = 0.7,raster.dpi = 300) +
    geom_point_rast(data = df_square %>% dplyr::filter(col2_val == 1),
                    color = col2_color, size = 6, alpha = 0.7,raster.dpi = 300) +
    geom_point_rast(data = green_cells,
                    color = col3_color, size = 6, alpha = 0.7,raster.dpi = 300)
  
  # ➕ Add col4 only if provided
  if (!is.null(col4)) {
    p <- p +
      geom_point_rast(
        data = df_square %>% dplyr::filter(col4_val == 1),
        color = col4_color,
        size = 6,
        alpha = 0.7,raster.dpi = 300
      )
  }
  
  p +
    coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    labs(title = paste0("Sample ", sample_id),
         x = "X coordinate", y = "Y coordinate") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
}

p=plot_cell_neighborhood(obs_df%>%filter(sample_id=="B01809C2"), cell_index = "65575560697256",col1='Ex_k3_score_bin80', col2 = "AST_k6_score_bin80",
                         col3='Microglia_k8_score_bin80',col4 = 'Ex_k9_score_bin80',sample_id='B01809C2',square_size=250)
ggsave(file.path(SPATIAL_PLOTS_DIR, "zoom_colocalization_B01809C2.pdf"), p, width = 6, height = 6)



plot_full_slice <- function(
    df,
    col1,
    col2,
    col3,
    col4 = NULL,
    sample_id,
    expand_factor = 0.05
) {
  
  df <- df %>%
    mutate(
      col1_val = factor(.data[[col1]], levels = c(0,1)),
      col2_val = factor(.data[[col2]], levels = c(0,1)),
      col3_val = factor(.data[[col3]], levels = c(0,1)),
      col4_val = if (!is.null(col4)) factor(.data[[col4]], levels = c(0,1)) else NA
    )
  
  
  # Colors — MATCH previous function
  col1_color <- "#c7522a"
  col2_color <- "#7B6FD0"
  col3_color <- "#FFC107"
  col4_color <- "#008585"  
  
  cells_list <- list(
    df %>% filter(col1_val == 1) %>% mutate(marker = col1),
    df %>% filter(col2_val == 1) %>% mutate(marker = col2),
    df %>% filter(col3_val == 1) %>% mutate(marker = col3)
  )
  
  color_map <- c(col1_color, col2_color, col3_color)
  names(color_map) <- c(col1, col2, col3)
  
  if (!is.null(col4)) {
    cells_list[[4]] <- df %>% filter(col4_val == 1) %>% mutate(marker = col4)
    color_map[col4] <- col4_color
  }
  
  cells_all <- bind_rows(cells_list)
  
  ggplot(df, aes(x, y)) +
    geom_point_rast(color = "grey80", size = 0.5, alpha = 0.5,raster.dpi = 300) +
    geom_point_rast(data = cells_all, aes(color = marker),
                    size = 0.5, alpha = 0.6,raster.dpi = 300) +
    scale_color_manual(
      name = "Cells",
      values = color_map
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2.5)
    )) +
    coord_equal() +
    labs(
      title = paste0("Sample ", sample_id),
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    expand_limits(
      x = c(min(df$x) - diff(range(df$x)) * expand_factor,
            max(df$x) + diff(range(df$x)) * expand_factor),
      y = c(min(df$y) - diff(range(df$y)) * expand_factor,
            max(df$y) + diff(range(df$y)) * expand_factor)
    )
}

specific_sample="B01809C2"
p=plot_full_slice(
  df = obs_df%>%filter(sample_id==specific_sample),
  col1 = "Ex_k3_score_bin80",
  col2 = "AST_k6_score_bin80",
  col3 = "Microglia_k8_score_bin80",
  col4="Ex_k9_score_bin80",
  sample_id = specific_sample
)
ggsave(file.path(SPATIAL_PLOTS_DIR, "full_slice_colocalization_B01809C2.pdf"), p, width = 8, height = 6)
local({
  d <- obs_df %>% filter(sample_id == specific_sample)
  c0 <- d["65575560697256", c("x", "y")]
  sd_add(d %>% tibble::rownames_to_column("cell") %>%
           mutate(in_zoom = abs(x - c0$x) <= 250 & abs(y - c0$y) <= 250, zoom_centre = cell == "65575560697256") %>%
           dplyr::select(cell, sample_id, annotation, x, y, Ex_k3_score_bin80, AST_k6_score_bin80, Microglia_k8_score_bin80, Ex_k9_score_bin80, in_zoom, zoom_centre),
         "Fig6", "e", description = "Sample B01809C2, all cells (left). 1 = cell in program: Ex_k3 = EXC.p3 stress/apoptosis, AST_k6 = AST.p6, Microglia_k8 = MIC.p8, Ex_k9 = EXC.p9 cell-cycle reentry. Right (zoom): in_zoom cells, a square of +-250 units around zoom_centre (arrow); dashed circles of radius 200 around Microglia_k8 cells")
})



################################################
#### Supp Fig 8- cell type distribution ########
################################################
df_plot <- obs_df %>%
  dplyr::select(sample_id, annotation, levels) %>%
  count(levels, sample_id, annotation)

df_plot$levels <- factor(
  df_plot$levels,
  levels = c("control", "moderate", "advanced", "severe")
)


sample_order <- df_plot %>%
  distinct(sample_id, levels) %>%
  arrange(levels) %>%
  pull(sample_id)

df_plot$sample_id <- factor(df_plot$sample_id, levels = sample_order)

celltype_colors <- c(
  "INH"        = "#c7522a",  
  "Ex"         = "#AD9764",  
  "OPC"        = "#586f44",  
  "OLIGO"      = "#74a892",  
  "Microglia"  = "#008585",  
  "AST"        = "#004343",  
  "END"        = "#B9B0A2",  
  "other"      = "#9A9A9A"   
)


df_plot <- df_plot %>%
  mutate(sample_label = paste(sample_id, levels, sep = " | "))

df_plot$sample_label <- factor(
  df_plot$sample_label,
  levels = unique(paste(sample_order, df_plot$levels[match(sample_order, df_plot$sample_id)], sep = " | "))
)

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "cell_type_distribution_per_sample.pdf"),width = 8,height = 6,{
  print(
    ggplot(df_plot, aes(x = sample_label, y = n, fill = annotation)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = celltype_colors) +
      labs(
        x = "Sample | Diagnosis",
        y = "Number of cells"
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()
      )
    
  )
})
sd_add(df_plot %>% dplyr::select(sample_id, levels, annotation, n), "ED8", "h",
       description = "Number of cells per sample and cell type (levels = diagnosis group)")

##########################################
#### Supp Fig 8- scores distributions ####
##########################################


x <- obs_df$Ex_k3_score
x <- x[!is.nan(x)]
qs <- quantile(x, probs = c(0.75, 0.80, 0.85), na.rm = TRUE)
q_df <- data.frame(
  q = qs,
  label = c("75%", "80%", "85%")
)

withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "Ex_k3_score_density_with_quantiles.pdf"),width = 6,height = 4,{
  print(ggplot(data.frame(x = x), aes(x)) +
          geom_density(
            linewidth = 0.8,
            color = "black",
            #bins=50,
            fill = "#74a892",
            alpha = 0.6
          ) +
          geom_vline( 
            data = q_df,
            aes(xintercept = q),
            linetype = "dashed",
            linewidth = 0.8
          ) +
          geom_text(
            data = q_df,
            aes(x = q, y = 0, label = label),
            angle = 90,
            vjust = -0.4,
            size = 3
          ) +
          theme_classic(base_size = 12) +
          labs(
            x = "Ex_k3 score",
            y = "Density"
          ))
  
})
sd_add(data.frame(Ex_k3_score = x), "ED8", "i", "scores", description = "Stress/apoptosis program (Ex_k3) signature score of every cell with a score (density)")
sd_add(q_df, "ED8", "i", "quantiles", description = "Dashed lines: 75/80/85% quantiles of the score")



# pathway example
pt <- obs_df$Ex_response_to_oxidative_stress_score
pt <- pt[!is.nan(x)]


withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "Ex_oxidative_stress_density_with_quantiles.pdf"),width = 6,height = 4,{
  print(
    ggplot(data.frame(x = pt), aes(x)) +
      geom_density(
        linewidth = 0.8,
        #bins = 50,
        color = "black",
        fill = "#74a892",
        alpha = 0.6
      ) +
      theme_classic(base_size = 12) +
      geom_vline(
        xintercept=0.29, #taken from mean of the second guassian distribution 
        linetype = "dashed",
        linewidth = 0.8
      )
    
  )
})
sd_add(data.frame(Ex_response_to_oxidative_stress_score = pt[!is.na(pt)]), "ED8", "j",
       description = "Response to oxidative stress pathway signature score of every excitatory cell (density; NaN rows, which ggplot drops, removed); dashed line at 0.29")

# TODO: add correlation heatmap

###############################################
### Bin-program proportion boxplots (case vs control) ###
###############################################

#' Plot per-sample proportions for a list of spatial bin programs.
#'
#' For each program column in `programs` that matches `bin{bin_num}`:
#'   - Cells with value == 1 define the "detected" population.
#'   - Proportion = (detected cells per sample) / (total cells of that
#'     annotation in that sample).
#'   - Groups are compared (control vs case) with a Mann-Whitney U test
#'     (wilcox.test), and p-values are FDR-corrected across all programs.
#'   - Returns a patchwork / faceted ggplot.
#'
#' @param obs_df        data.frame with one row per cell. Must contain columns:
#'                        sample_id, diagnosis ("control" / "case"), annotation,
#'                        and all program columns listed in `programs`.
#' @param programs      character vector of bin-column names to include
#'                        (e.g. c("Ex_k3_score_bin80", "AST_k6_score_bin80")).
#' @param bin_num       string; only columns containing "bin{bin_num}" are kept
#'                        (default "80").
#' @param max_cols_per_row  panels per row in the output grid (default 4).
plot_bin_proportion_boxplots <- function(
    obs_df,
    programs,
    bin_num          = "80",
    max_cols_per_row = 4
) {

  bin_pattern <- paste0("bin", bin_num)
  cols_to_plot <- programs[grepl(bin_pattern, programs, fixed = TRUE)]

  if (length(cols_to_plot) == 0) {
    message("No columns matching bin", bin_num, " found in programs.")
    return(invisible(NULL))
  }

  # ---- Phase 1: compute proportions + raw p-values ----
  results <- lapply(cols_to_plot, function(col_i) {
    if (!col_i %in% colnames(obs_df)) {
      warning("Column not found in obs_df: ", col_i)
      return(NULL)
    }

    detected_annotations <- obs_df %>%
      filter(.data[[col_i]] == 1) %>%
      pull(annotation) %>%
      na.omit() %>%
      unique()

    if (length(detected_annotations) == 0) return(NULL)

    annotation_totals <- obs_df %>%
      filter(annotation %in% detected_annotations) %>%
      group_by(sample_id) %>%
      summarise(total_annotation_cells = n(), .groups = "drop")

    proportions_df <- obs_df %>%
      filter(.data[[col_i]] == 1) %>%
      group_by(sample_id, diagnosis) %>%
      summarise(bin_cell_count = n(), .groups = "drop") %>%
      left_join(annotation_totals, by = "sample_id") %>%
      mutate(proportion = bin_cell_count / total_annotation_cells)

    control_vals <- proportions_df %>% filter(diagnosis == "control") %>% pull(proportion)
    case_vals    <- proportions_df %>% filter(diagnosis == "case")    %>% pull(proportion)

    if (length(control_vals) == 0 || length(case_vals) == 0) return(NULL)

    p_val <- wilcox.test(control_vals, case_vals, alternative = "two.sided", exact = FALSE)$p.value

    list(
      col_i          = col_i,
      proportions_df = proportions_df,
      annotations    = detected_annotations,
      p_raw          = p_val
    )
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    message("No valid comparisons found.")
    return(invisible(NULL))
  }

  # ---- Phase 2: FDR correction ----
  p_raw <- sapply(results, `[[`, "p_raw")
  p_adj <- p.adjust(p_raw, method = "BH")

  sig_label <- function(p) {
    dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
  }

  # ---- Phase 3: build combined long data.frame ----
  plot_df <- purrr::imap_dfr(results, function(res, idx) {
    annotation_str <- paste(head(res$annotations, 2), collapse = ", ")
    res$proportions_df %>%
      mutate(
        program   = res$col_i,
        facet_label = paste0(res$col_i, "\n(", annotation_str, ")"),
        adj_p     = p_adj[[idx]],
        sig       = sig_label(p_adj[[idx]])
      )
  })

  # significance bracket positions per facet
  bracket_df <- plot_df %>%
    group_by(facet_label, sig, adj_p) %>%
    summarise(
      y_max   = max(proportion, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_lift  = pmax(y_max * 0.05, 0.005),
      bar_y   = y_max + y_lift,
      text_y  = bar_y + y_lift * 0.3,
      y_upper = text_y + y_lift * 1.5,
      label   = ifelse(
        sig == "ns",
        "n.s.",
        paste0(sig, "\n(FDR p=", formatC(adj_p, digits = 4, format = "f"), ")")
      )
    )

  plot_df <- plot_df %>%
    left_join(bracket_df %>% select(facet_label, y_upper), by = "facet_label")

  # ---- Phase 4: plot ----
  n_cols   <- length(results)
  diagnosis_order <- c("control", "case")
  plot_df  <- plot_df %>%
    mutate(diagnosis = factor(diagnosis, levels = diagnosis_order))

  ggplot(plot_df, aes(x = diagnosis, y = proportion)) +
    geom_boxplot(
      outlier.shape = NA,
      width         = 0.5,
      fill          = "white",
      color         = "black"
    ) +
    geom_jitter(
      width = 0.15,
      size  = 2,
      alpha = 0.7,
      color = "black"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = 1.5, y = text_y, label = label),
      inherit.aes = FALSE,
      size        = 3,
      vjust       = 0
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = 1, xend = 1, y = bar_y, yend = bar_y + y_lift * 0.2),
      inherit.aes = FALSE, linewidth = 0.5
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = 2, xend = 2, y = bar_y, yend = bar_y + y_lift * 0.2),
      inherit.aes = FALSE, linewidth = 0.5
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = 1, xend = 2, y = bar_y + y_lift * 0.2, yend = bar_y + y_lift * 0.2),
      inherit.aes = FALSE, linewidth = 0.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, NA), clip = "off") +
    facet_wrap(~facet_label, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid   = element_blank(),
      strip.text   = element_text(size = 8),
      axis.text.x  = element_text(angle = 20, hjust = 1),
      legend.position = "none",
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black")
    ) +
    labs(
      x = "Diagnosis",
      y = "Proportion within Cell Type"
    )
}
withr::with_pdf(file.path(SPATIAL_PLOTS_DIR, "bin_proportion_boxplots_80.pdf"), width=6, height=12, { 
  g <- plot_bin_proportion_boxplots(obs_df, programs = c("Ex_k3_score_bin80", "AST_k6_score_bin80", "Microglia_k8_score_bin80", "Ex_k9_score_bin80"), bin_num = "80")
  print(g)
  sd_add(g$data %>% dplyr::select(program, sample_id, diagnosis, bin_cell_count, total_annotation_cells, proportion, adj_p, sig), "ED10", "d",
         description = "Per sample: cells in the program (bin_cell_count) / cells of that cell type (total_annotation_cells) = proportion; Wilcoxon test control vs case, BH across the 4 programs (adj_p)")
})

message("Spatial analysis complete. Plots saved to: ", SPATIAL_PLOTS_DIR)

