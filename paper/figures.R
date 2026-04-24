library(rhdf5)
library(dplyr)
library(ComplexHeatmap)
library(slanter)
library(tidyr)
library(ggpubr)
library(ggpointdensity)

source("R/data_500_utils.R")
source("R/figures.funs.R")
source("R/landscape.analysis.funs.R")
source("R/proteomics.R")
source("R/topic_space_functions.R")
source("R/pathways.analysis.R")
source("R/pathway_dynamics.R")
source("R/utils.R")
source("R/visualizations.R")
source("validation/validation.functions.R")


options(ggrastr.default.dpi=900)

embed.width <- embed.height <- 5
embed.width.small <- embed.height.small <- 2.5

white2orange <- colorRampPalette(c("white", "#c7522a"))
blue2orange <- colorRampPalette(c("#008585", "white", "#c7522a"))
cyan2orange <- colorRampPalette(c("#4db6ac", "white", "#ef6c00"))

colorRange <- function(col, n=4) {
  colorRampPalette( colorRampPalette(c("white", col))(n)[(n-1):(n)] )
}

cell.type.to.prefix <- list(
  'astrocytes'='AST',
  'opcs'='OPC',
  'oligo'='OLI',
  'microglia'='MIC',
  'cux2p'='EXC L2-3',
  'cux2m'='EXC L4-6',
  'inhibitory'='INH'
)



AD_TRAITS_NAMES <- c(amyloid_mf_sqrt="Neocortical amyloid",
                     sqrt_amyloid_mf="Neocortical amyloid",
                     sqrt.amyloid_mf="Neocortical amyloid",
                     tangles_mf_sqrt="Neocortical tangles",
                     sqrt_tangles_mf="Neocortical tangles",
                     sqrt.tangles_mf="Neocortical tangles",
                     sqrt.amyloid="amyloid",
                     sqrt.tangles="tangles",
                     cogng_random_slope="Cognitive decline",
                     cogng_demog_slope="Cognitive decline")

figures_dir <- file.path(get_result_dir(), "figures")
dir.create(figures_dir, showWarnings = FALSE)


topics_h5ad_path <- file.path(get_result_dir(), 'topic_space', '500_topics_all.h5ad')

data <- anndata::read_h5ad(topics_h5ad_path)
data.da <- read_h5ad(AD_500_ANDATA_FILE)



#+++++++++++++++++++++++++++
# Figure 1 -----------------
#+++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++
## Figure 1b Structure plot ----
#++++++++++++++++++++++++++++++++

function() {
  # https://colorkit.co/palette/c7522a-e5c185-f0daa5-fbf2c4-b8cdab-74a892-008585-004343/
  colors <- list('inhibitory'="#c7522a", 'cux2p'="#B89658", 'cux2m'="#AD9764", 'opcs'="#586f44",
                 'oligo'="#74a892", 'microglia'= "#008585", 'astrocytes'="#004343")

  ks <- list('inhibitory'="17",
             'cux2p'="14",
             'cux2m'="16",
             'opcs'="7",
             'oligo'="8",
             'microglia'="15",
             'astrocytes'="10")
  fits <- sapply(names(ks), function(ct) {
    L_mat <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[ct], "/L"))

    fit <- list()
    fit$L <- L_mat %>% column_to_rownames('cell') %>% as.matrix()
    fit$F <- matrix(1, ncol = ncol(fit$L))
    class(fit) <- c('poisson_nmf_fit', list)
    
    fit
  }, simplify = FALSE)

  p <- lapply(names(colors), function(ct) {
    # make this graph reproducible
    set.seed(5)
    print(ks[[ct]])
    k.int = as.integer(ks[[ct]])
    cols <- colorRampPalette(colorRampPalette(c("white", colors[[ct]]))(7)[c(2,7)])(k.int)
    # change the order of dark and light colors
    cols <- sample(cols)
    title <- switch(ct, "cux2+"="Exc. L2-3", "cux2-"="Exc. L4-6", opc="OPC", str_to_title(ct))
    ggrastr::rasterise(fastTopics::structure_plot(fits[[ct]], colors=cols,
                   n=5000
    )) + theme_void() + labs(title=title) + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
  }) %>% plot_grid(plotlist=., nrow = 1)

  withr::with_pdf(file.path(figures_dir, "1b_rasterise.pdf"), width = embed.width*2, height = embed.height, {
    print(p)
  })

  rm(fits)
}

## Figure 1c UMAP, Extended Data Figure 2b,e 3a ---------

(function() {
  cell_types_to_topics <- list('astrocytes'=c("p2", "p5", "p6"),
                               'microglia'=c("p9", "p8", "p5"),
                               'opcs'=c("p1", "p2", "p4"),
                               'oligo'=c("p5", "p6", "p8"))

  scales <- list('microglia'=scale_x_continuous(trans = squash_axis(-10, -6, 98)))
  for (ct in names(cell_types_to_topics)) {
    embedding <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/umap.densities"))

    embedding <- embedding %>% rename_with(~gsub("X\\d+.k", "p", .x) %>% gsub(".sim.umap", "", .))

    if (ct %in% names(scales)) {
      scale <- scales[[ct]]
    } else {
      scale <- scale_x_continuous()
    }

    density_plots <- lapply(cell_types_to_topics[[ct]], function (k) plot_umap_density2(embedding, columns = k, title=glue::glue('program {k} density')) + scale)

    withr::with_pdf(file.path(figures_dir, glue::glue("1C_{ct}_umap.pdf")), height = embed.height, width = embed.width * 1.4, {
      print(plot_grid(plot_umap_state(embedding) + scale, plotlist = density_plots, nrow = 2))
    })

    withr::with_pdf(file.path(figures_dir, glue::glue("1C_{ct}_umap_density_only.pdf")), height = embed.height * 1.5, width = embed.width*0.7, {
      print(plot_umap_density2(embedding = embedding, columns = cell_types_to_topics[[ct]], ncol = 1, title = glue::glue("{ct} programs")) + scale)
    })
  }
})()


#+++++++++++++++++++++++++++
## Extended Data Figure 2c,f 3b, 4a-d - F matrix
#+++++++++++++++++++++++++++

(function() {
  shared_neuronal_genes <- c("PRDX5", "MRPL51", "ANAPC4", "HSPH1", "CENPC", "SNAP25", "MT1G", "GRINA", "SYNGR1", "SLC12A5", "SYT9", "SLC25A17")
  cell_type_to_genes <- list(
                            astrocytes=c("GFAP", "SLC1A2", "SERPINA3", "OSMR", "SLC38A2", "HSPH1"),
                            microglia=c("PTPRG", "TOP2A", "IFI6", "MRC1", "CPM", "APOE", "TREM2", "GPNMB", "IL18", "CX3CR1", "SOAT1"),
                            opcs=c("PINK1", "SERPINA3", "OSMR", "VCAN", "PDGFRA", 'FRMD4A'),
                            oligo=c("OSMR", "MBP", "MAL", "SLC38A2", "MOG", "HSPH1", "DNAJB1", "IGF1R", "QDPR"),
                            cux2p=c("RORB", "THEMIS", "FEZF2", "FOS",  shared_neuronal_genes),
                            cux2m=c("RORB", "THEMIS", "FEZF2", "FOS",  shared_neuronal_genes),
                            inhibitory=c("SST", "VIP", "PVALB", "NPY", "KIT", "LAMP5", "FOS", shared_neuronal_genes)
                            )


  for (cell.type in names(cell_type_to_genes)) {
    des <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/de"))
    F_mat <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/F"))
    F_mat <- F_mat %>% dplyr::rename_with(~gsub("^k", "p", .x))
    withr::with_pdf(file.path(figures_dir,  glue::glue("2SuplC_{cell.type}_F.pdf")), {
      plot_F_matrix(F_mat, des, col=cyan2orange(10),
                    column_title=glue::glue("{str_to_sentence(cell.type)} programs gene score"),
                    marked_genes = cell_type_to_genes[[cell.type]],
                    ) %>% print()
    })
  }

  rm(des, F)

})()


#++++++++++++++++++++++++++++++++++++#
## Figure 1d, Extended Data Figure 1d,g, Extended Data Figure 3c - State program score ----
#++++++++++++++++++++++++++++++++++++#

(function (){
  programs_order <- list(
    astrocytes=c("5", "3", "8","2","10", "4","6", "1", "7", "9"),
    opcs=c(3, 7, 5, 6, 4, 2, 1),
    microglia=c(3, 9, 7, 6, 5, 15, 13, 11, 8, 4, 1, 2, 10, 12, 14),
    oligo=c(1, 6, 3, 5, 2, 8, 7, 4),
  )

  for (ct in names(cell.type.to.prefix)) {
    abundance.data <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/abundance.state"))
    abundance.data <- abundance.data %>% mutate(variable=paste0(cell.type.to.prefix[[ct]],".p", variable))

    abundance.data <- abundance.data %>%
      pivot_wider(id_cols = "state", names_from = "variable", values_from = "avg", values_fill = 0) %>%
      arrange(str_extract(state, "\\d+") %>% as.integer()) %>%
      tibble::column_to_rownames("state") %>%
      as.matrix()

    withr::with_pdf(file.path(figures_dir, glue::glue("1d_programs_vs_cluster_{ct}.pdf")), width = embed.width, height = embed.height, {
      Heatmap(abundance.data,
              column_order = paste0(cell.type.to.prefix[[ct]],".p", programs_order[[ct]]),
              row_title = 'States',
              column_title = 'Programs',
              cluster_rows = F,
              cluster_columns = F,
              heatmap_legend_param = list(
                title = "Program mean score",
                legend_direction = "horizontal"
              ),
              col = circlize::colorRamp2(c(0, max(abundance.data)), c("white", "#c7522a"))
      ) %>% draw(heatmap_legend_side = "bottom")
    })
  }

})()


## Figure 1f - Neurons umap --------
(function () {
  MIN_NUMBER_OF_CELLS_LABEL <- 1000

  ct_to_topics <- list("inhibitory"=c("X17.k10.sim.umap", "X17.k6.sim.umap"), "cux2m"=c("X16.k3.sim.umap"), "cux2p"=c("X14.k11.sim.umap"))
  ct_to_group_by <- c("inhibitory"="pruned_subtype", "cux2m"="pruned_subtype", "cux2p"="marker.1")
  for (ct in names(ct_to_topics)) {
    allen <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/allen.annotation")) %>% tibble::column_to_rownames("cell")
    embedding <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/umap.densities")) %>% tibble::column_to_rownames("cell")

    allen <- allen %>% tidyr::separate(col = "pruned_subclass", sep = " ", into = c(NA, NA, "marker.1", NA), extra = "merge", remove = FALSE)

    embedding <- cbind(embedding, allen)
    embedding <- embedding %>%
      dplyr::rename(x=UMAP_1, y=UMAP_2) %>%
      mutate(pruned_subtype=ifelse(pruned_subtype == 'NA', NA, pruned_subtype))


    gb <- ct_to_group_by[[ct]]

    # Plot UMAP embedding of atlas
    withr::with_pdf(file.path(figures_dir, glue::glue("1f_{ct}_{gb}.pdf")), {
      print(ggplot(embedding, aes(x, y, color=!!sym(gb), label=!!sym(gb))) +
              ggrastr::geom_point_rast(size=.01, raster.dpi = 800) +
              geom_text(data = embedding %>%
                          group_by(!!sym(gb)) %>%
                          summarise(across(c(x, y), mean), n = n()) %>%
                          filter(n > MIN_NUMBER_OF_CELLS_LABEL), color="black", size=6) +
              no.labs +
              theme_embedding +
              theme(legend.position = "left") +
              guides(color = guide_legend(override.aes = list(size = 3))))
    })
    topics <- ct_to_topics[[ct]]
    for (topic in topics) { 
      topic_number <- str_extract(topic, "k(\\d+)") %>% gsub("k", "", .)
      title <- glue::glue('{ct} prgoram P{topic_number}')

      withr::with_pdf(file.path(figures_dir, glue::glue("1f_{ct}_program_{topic_number}.pdf")), {
        print(ggplot(embedding, aes_string('x', 'y', color=topic)) +
          ggrastr::geom_point_rast(size=0.01, raster.dpi = 512) +
          scale_color_viridis_c(option='turbo') +
          theme_embedding +
          theme_void() +
          ggtitle(title) +
          labs(color='program density'))
      })
    }
  }
})

## Figure 1e Neurons sub-type -------------------

(function () {
  class_to_panel <- c("inhibitory"='e')
  class_to_plot_dim <- list("inhibitory"=list(w=embed.width * 0.8, h=embed.height))
  for (neuronal_class in c("inhibitory")) {
    if (neuronal_class == 'inhibitory') {
      agg <- 'class'
      df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["inhibitory"]], glue::glue("/abundance.allen.{agg}")))
      df <- df %>% mutate(topic=paste0("INH.p", topic))
    } else {

      df_m <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2m"]], glue::glue("/abundance.allen")))
      df_p <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2p"]], glue::glue("/abundance.allen")))


      df_m <- df_m %>%
        tidyr::separate(col = "label", sep = " ", into = c("type", "layer", "marker.1", "marker.2"), extra = "merge", remove = FALSE) %>%
        group_by(layer, marker.1, topic) %>% summarise(avg = sum(avg * n) / sum(n), n=sum(n)) %>% ungroup

      df_p <- df_p %>%
        tidyr::separate(col = "label", sep = " ", into = c("type", "layer", "marker.1", "marker.2"), extra = "merge", remove = FALSE) %>%
        group_by(layer, marker.1, topic) %>% summarise(avg = sum(avg * n) / sum(n), n=sum(n)) %>% ungroup

      df_p <- df_p %>% mutate(topic=paste0("EXC L2-3.p", topic), cux2='p')
      df_m <- df_m %>% mutate(topic=paste0("EXC L4-6.p", topic), cux2='m')

      df <- rbind(df_m, df_p)
      agg <- "layer_marker.1"
      fractions <- df %>% distinct(layer, marker.1, cux2, n) %>% group_by(layer, marker.1) %>% mutate(frac=n/sum(n)) %>% ungroup %>% arrange(layer, marker.1)

      df <- df %>% dplyr::left_join(fractions, by=c("layer", "marker.1", "cux2", "n")) %>% mutate(avg=avg*frac) %>% mutate(layer_marker.1=paste(layer, marker.1))
    }

    df <- pivot_wider(df, id_cols=all_of(agg), names_from = "topic", values_from = "avg", values_fill = 0) %>%
      column_to_rownames(paste(agg)) %>%
      as.matrix %>%
      t


    class_order <- c(c("CHANDELIER", "PVALB", "LAMP5 LHX6", "LAMP5", "SNCG",
                     "PAX6", "VIP", "SST", "SST CHODL"),
                     c("L2 LAMP5", "L3 LAMP5",
                       "L2 LINC00507", "L2-3 LINC00507",
                       "L3 THEMIS",
                       "L2-3 RORB",  "L3 RORB", "L3-5 RORB", "L5 RORB",
                       "L3-5 FEZF2", "L5 FEZF2", "L5-6 FEZF2", "L6 FEZF2",
                       "L5 THEMIS", "L5-6 THEMIS", "L6 THEMIS"
                       ))


    programs_order <- c(paste0("INH.p", c(17, 6, 4, 12, 3, 15,
                                        13, 5, 9, 2, 16, 1,7,8,10,11,14)),
                        paste0("EXC L2-3.p", c(1,4,9,10,11,14, 13, 5, 6, 12, 2,8, 7, 3)),
                        paste0("EXC L4-6.p", c(10, 1, 12, 15, 16, 6, 13, 2, 14, 11, 3, 4, 5, 7 ,8, 9))
                        )


    withr::with_pdf(file.path(figures_dir, glue::glue("1{class_to_panel[[neuronal_class]]}_allen_label_{agg}.pdf")),
                    height = class_to_plot_dim[[neuronal_class]]$h,
                    width = class_to_plot_dim[[neuronal_class]]$w, {
      Heatmap(df,
              col=test.col,
              cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              cluster_columns = FALSE,
              row_split = rownames(df) %>% str_replace("\\.p\\d+", ""),
              column_order = class_order[class_order %in% colnames(df)],
              row_order = programs_order[programs_order %in% rownames(df)],
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              row_gap = unit(5, "mm"),
              heatmap_legend_param = list(
                title = "Program mean score",
                legend_direction = "horizontal"
              ),
              show_row_names = TRUE, name = "Mean program score") %>% draw(heatmap_legend_side = "bottom")
      })
  }
})()


# Extended Data Figure 1 ----------------------

## Extended Data Figure 1a,b -----------

(function (){
for (ct in names(CELL_TYPE_TO_K)) {
    for (k_diff in c(-1, 1)) {
      message(glue::glue("Plotting L correlation for {ct} diff of {k_diff}"))
      k <- CELL_TYPE_TO_K[[ct]]
      other_k <- as.character(as.integer(k) + k_diff)

      L <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L")) %>% column_to_rownames("cell") %>% rename_with(~ gsub("^k", "p", .x))
      L_other_k <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/fit_", other_k, "/L")) %>% column_to_rownames("cell")
      L_other_k <- L_other_k %>% rename_with(~ gsub("^k", "p", .x))
      c <- correlate_dataframes(L, L_other_k)

      add_significance <- TRUE
      if (add_significance) {
        corr.sig <- c$sig
        corr.sig[corr.sig != ''] <- '*'
        cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)
      } else {
        cell_fun <- NULL
      }

      similarity <- c$corr
      similarity[similarity < 0] <- 0
      order <- slanter::slanted_orders(similarity)

      legend <- Legend(title = "adj.pval", pch = c("*"), type = "points", labels = c("<0.05"), title_gp = gpar(fontface='plain'))
      correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

      withr::with_pdf(file.path(figures_dir, glue::glue("1SuplA_{ct}_{k_diff}.pdf")), {
        ComplexHeatmap::Heatmap(c$corr,
                                heatmap_legend_param = list(title="Programs correlation", direction="vertical"),
                                col = cyan2orange.correlation,
                                clustering_distance_rows = "spearman",
                                clustering_distance_columns = "spearman",
                                show_row_dend = FALSE,
                                show_column_dend = FALSE,
                                row_order = order$rows,
                                column_order = order$cols,
                                column_names_rot = 90,
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                row_title = paste0(k, ' programs'),
                                column_title = paste0(other_k, ' programs'),
                                column_title_side = 'bottom',
                                show_heatmap_legend =FALSE,
                                cell_fun = cell_fun) %>% draw(annotation_legend_list=list(legend, correlation.legend),
                                                              column_title=glue::glue("Different k {ct} cell score correlation"))
      })

    }
  }
})()

## Extended Data Figure 1c -----------

(function (){
  for (ct in c("inhibitory", "astrocytes")) {
    L_all <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L")) %>% column_to_rownames("cell") %>% rename_with(~ gsub("^k", "p", .x))
    L_subset <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/random_subset/L")) %>% column_to_rownames("cell")
    L_subset <- L_subset %>% rename_with(~ gsub("^k", "p", .x))
    c <- correlate_dataframes(L_all[rownames(L_subset),], L_subset)

    add_significance <- TRUE
    if (add_significance) {
      corr.sig <- c$sig
      corr.sig[corr.sig != ''] <- '*'
      cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)
    } else {
      cell_fun <- NULL
    }

    similarity <- c$corr
    similarity[similarity < 0] <- 0
    order <- slanter::slanted_orders(similarity)

    legend <- Legend(title = "adj.pval", pch = c("*"), type = "points", labels = c("<0.05"), title_gp = gpar(fontface='plain'))
    correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

    withr::with_pdf(file.path(figures_dir, glue::glue("1SuplC_{ct}.pdf")), {
      ComplexHeatmap::Heatmap(c$corr,
                              heatmap_legend_param = list(title="Programs correlation", direction="vertical"),
                              col = cyan2orange.correlation,
                              clustering_distance_rows = "spearman",
                              clustering_distance_columns = "spearman",
                              show_row_dend = FALSE,
                              show_column_dend = FALSE,
                              row_order = order$rows,
                              column_order = order$cols,
                              column_names_rot = 90,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              row_title = 'All cells fit',
                              column_title = 'Random subset cells fit',
                              column_title_side = 'bottom',
                              show_heatmap_legend =FALSE,
                              cell_fun = cell_fun) %>% draw(annotation_legend_list=list(legend, correlation.legend), column_title=glue::glue("Cells correlation {ct} all and random subset"))
    })

  }
})()


## Extended Data Figure 4e - trait association to qc, mean count, and number of feature ----------------------
(function() {
  withr::with_pdf(file.path(figures_dir, "4SuplE_trait_association_to_qc.pdf"), {
    QC_NAMES <- c("Total_Genes_Detected" = "Total genes detected",
     "Estimated_Number_of_Cells" = "Estimated number of cells", "Median_UMI_Counts_per_Cell" = "Median UMI counts per cell")
    plot.trait.associations.heatmap(data$uns$trait.analysis.to.sample.qc, 
                                    params=c("sqrt.amyloid_mf",  "sqrt.tangles_mf", "cogng_demog_slope"),
                                    column_labels=AD_TRAITS_NAMES,
                                    row_labels=QC_NAMES,
                                    row.by='covariate', show.only.significant = F)

    plot.trait.associations.heatmap(data$uns$trait.analysis.to.mean.count, 
                                    params=c("sqrt.amyloid_mf",  "sqrt.tangles_mf", "cogng_demog_slope"),
                                    column_labels=AD_TRAITS_NAMES, row_title="Mean count",
                                    row.by='covariate', show.only.significant = F)

    plot.trait.associations.heatmap(data$uns$trait.analysis.to.number.of.features, 
                                    params=c("sqrt.amyloid_mf",  "sqrt.tangles_mf", "cogng_demog_slope"),
                                    column_labels=AD_TRAITS_NAMES, row_title="Number of features",
                                    row.by='covariate', show.only.significant = F)
  })
})()



## Extended Data Figure 1d - program in neuronal subtypes -----------

(function (){
  ct <- "inhibitory"

  L_all <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L")) %>%
    column_to_rownames("cell") %>%
    rename_with(~ gsub("^k", "INH.p", .x))

  L_sst <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/subtype_SST", "/L")) %>%
    column_to_rownames("cell") %>%
    rename_with(~ gsub("^k", "SST.p", .x))

  L_pv <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/subtype_PVALB", "/L")) %>%
    column_to_rownames("cell") %>%
    rename_with(~ gsub("^k", "PV.p", .x))

  c_sst <- correlate_dataframes(L_all[rownames(L_sst),], L_sst)
  c_pv <- correlate_dataframes(L_all[rownames(L_pv),], L_pv)

  c <- list(
    corr=cbind(c_sst$corr, c_pv$corr),
    adj.pval=cbind(c_sst$adj.pval, c_pv$adj.pval),
    sig=cbind(c_sst$sig, c_pv$sig)
  )

  corr.sig <- c$sig
  corr.sig[c$adj.pval < 0.001] <- '*'
  corr.sig[c$adj.pval > 0.001] <- ''
  cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)

  similarity <- c$corr
  similarity[similarity < 0] <- 0
  order <- slanter::slanted_orders(similarity)

  legend <- Legend(title = "adj.pval", pch = c("*"), type = "points", labels = c("<0.001"), title_gp = gpar(fontface='plain'))
  correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

  withr::with_pdf(file.path(figures_dir, "1SuplD_inhibitory_subtype_fit.pdf"), {
    ComplexHeatmap::Heatmap(c$corr,
                            heatmap_legend_param = list(title="Programs correlation", direction="vertical"),
                            col = cyan2orange.correlation,
                            clustering_distance_rows = "spearman",
                            clustering_distance_columns = "spearman",
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            row_order = order$rows,
                            column_order = order$cols,
                            column_names_rot = 90,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            row_title = 'All cells fit',
                            column_title_side = 'bottom',
                            show_heatmap_legend =FALSE,
                            column_split=factor(c(rep('SST fit',ncol(L_sst)), rep('PV fit',ncol(L_pv))),
                                                levels = c('SST fit', 'PV fit')),
                            cell_fun = cell_fun) %>% draw(annotation_legend_list=list(legend, correlation.legend))
  })
})()


## Extended Data Figure 1e,f - scHPF programs -------------
(function() {
  for (ct in c("inhibitory", "astrocytes")) {
    L <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L")) %>% column_to_rownames("cell")
    L <- L %>% rename_with(~ gsub("^k", "p", .x))
    L_schpf <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/schpf/cell_score")) %>% column_to_rownames("cell")
    c <- correlate_dataframes(L[rownames(L_schpf),], L_schpf, method = "pearson")


    de_genes <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/de"))
    de_genes <- de_genes %>% mutate(topic=gsub("^k", "p", topic))
    gene_score <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/schpf/gene_score")) %>% column_to_rownames("gene")

    # Compute overlap for top N genes
    top_N <- 150
    gene_score_top_n <- sapply(colnames(gene_score), function(col) {
        gene_score %>%
            slice_max(order_by = !!sym(col), n = top_N) %>%
            rownames() %>% as_tibble()
    }, simplify = FALSE)

    des_nested <- de_genes %>%
        group_by(topic) %>%
        dplyr::select(gene, topic) %>%
        group_nest()

    gene_overlap_mat <- outer(des_nested %>% pull(data),
                              gene_score_top_n,
                Vectorize(tibble_overlap_coefficient))
    rownames(gene_overlap_mat) <- des_nested %>% pull(topic)
    gene_overlap_mat <- gene_overlap_mat[c$names[[1]], c$names[[2]]] 

    corr.sig <- c$sig
    corr.sig[c$adj.pval < 0.001] <- '*'
    corr.sig[c$adj.pval > 0.001] <- ''
    cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)

    similarity <- c$corr
    similarity[similarity < 0] <- 0
    order <- slanter::slanted_orders(similarity)

    legend <- Legend(title = "FDR", pch = c("*"), type = "points", labels = c("<0.001"), title_gp = gpar(fontface='plain'))
    correlation.legend <- Legend(title = "correlation pearson", col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

    withr::with_pdf(file.path(figures_dir, glue::glue("1SupE_scHPF_{ct}.pdf")), {
      ComplexHeatmap::Heatmap(c$corr,
                              heatmap_legend_param = list(title="Programs correlation", direction="vertical"),
                              col = cyan2orange.correlation,
                              clustering_distance_rows = "spearman",
                              clustering_distance_columns = "spearman",
                              show_row_dend = FALSE,
                              show_column_dend = FALSE,
                              row_order = order$rows,
                              column_order = order$cols,
                              column_names_rot = 90,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              row_title = 'fastTopics programs',
                              column_title = 'scHPF programs',
                              column_title_side = 'bottom',
                              show_heatmap_legend =FALSE,
                              cell_fun = cell_fun,
                              width = ncol(c$corr)*unit(5, "mm"),
                              height = nrow(c$corr)*unit(5, "mm")) %>% draw(annotation_legend_list=list(legend, correlation.legend), column_title=glue::glue("Cells correlation {ct} fastTopics and scHPF"))

      ComplexHeatmap::Heatmap(gene_overlap_mat,
                                col = circlize::colorRamp2(c(0, 1), white2orange(2)),
                                row_title = 'fastTopics programs',
                                row_order = order$rows,
                                column_order = order$cols,
                                column_title = 'scHPF programs',
                                column_title_side = 'bottom',
                                show_column_names = TRUE,
                                show_row_names = TRUE,
                                name = "Gene overlap",
                                show_row_dend = FALSE,
                                show_column_dend = FALSE,
                                width = ncol(gene_overlap_mat)*unit(5, "mm"),
                                height = nrow(gene_overlap_mat)*unit(5, "mm")) %>% draw
    })
  } 
})()

## Extended Data Figure 2a - projid abundance ----------------------
(function() {
  for (ct in names(CELL_TYPE_TO_GROUP)) {
    abundance.projid <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/abundance.projid"))
    withr::with_pdf(file.path(figures_dir, paste0("2SupA_projid_abundance_", ct, ".pdf")), height=embed.height, width=embed.width.small * 0.9, {
      print(abundance.projid %>%
            arrange(as.integer(variable)) %>%
            mutate(topic=paste0(cell.type.to.prefix[[ct]] ,".p", as.character(variable))) %>%
            mutate(topic=factor(topic, levels = unique(topic))) %>%
            ggplot(aes(topic, avg, color="#c7522a", fill="#c7522a")) +
            geom_violin(alpha=.3) +
            scale_color_manual(values="#c7522a") +
            scale_fill_manual(values="#c7522a") +
            scale_y_sqrt(breaks = c(.01, .1, .25, .5, .75), labels = paste0(as.integer(100*c(.01, .1, .25, .5, .75)), "%"), expand = expansion(0)) +
            labs(x=NULL, y=NULL) +
            theme_classic() +
            theme(legend.position = "none") +
            coord_flip()
      )
    })
  }
})()

#+++++++++++++++++++++++++++
# Extended Data Figure 4 -------------
#+++++++++++++++++++++++++++

## Extended Data Figure 4d,e --------------

NEURONAL_SUBTYPE_ORDER <- c("Inh L1 LAMP5 PVRL2", "Inh L1 LAMP5 RAB11FIP1", "Inh L1-6 LAMP5 AARD",
                            "Inh L1-6 LAMP5 NES", "Inh L1-6 LAMP5 CA1", "Inh L5-6 LAMP5 CRABP1",
                            "Inh L3-6 PAX6 LINC01497", "Inh L1 PAX6 MIR101-1", "Inh L1 PAX6 CHRFAM7A",
                            "Inh L1-6 VIP SLC7A6OS", "Inh L2 PAX6 FREM2", "Inh L1-2 VIP HTR3A",
                            "Inh L1-2 VIP WNT4", "Inh L1 LAMP5 BMP2", "Inh L1 LAMP5 NMBR",
                            "Inh L1 PVALB SST ASIC4", "Inh L1 SST P4HA3", "Inh L1 SST DEFB108B",
                            "Inh L1-3 VIP CBLN1", "Inh L1-2 VIP EXPH5", "Inh L1-3 VIP HSPB6",
                            "Inh L1-2 VIP PTGER3", "Inh L1-2 VIP SCML4", "Inh L1-3 VIP FNDC1",
                            "Inh L1-3 VIP CHRNA2", "Inh L2 VIP SLC6A16", "Inh L1 VIP KLHDC8B",
                            "Inh L3-5 VIP IGDCC3", "Inh L1-5 VIP SMOC1", "Inh L3-5 VIP HS3ST3A1",
                            "Inh L5-6 VIP COL4A3", "Inh L1-5 VIP CD27-AS1", "Inh L2-5 VIP BSPRY",
                            "Inh L2-5 VIP SOX11", "Inh L3-6 VIP ZIM2-AS1", "Inh L1-5 VIP PHLDB3",
                            "Inh L1-5 VIP LINC01013", "Inh L3-6 VIP UG0898H09", "Inh L3-5 VIP TAC3",
                            "Inh L1-6 SST NPY", "Inh L3-5 SST CDH3", "Inh L5 SST RPL35AP11",
                            "Inh L5-6 SST ISX", "Inh L5-6 PVALB SST CRHR2", "Inh L3-5 SST OR5AH1P",
                            "Inh L3-5 SST GGTLC3", "Inh L2-3 SST NMU", "Inh L5-6 SST PAWR",
                            "Inh L5-6 SST PIK3CD", "Inh L1-2 SST CCNJL", "Inh L1-3 SST FAM20A",
                            "Inh L1-2 SST PRRT4", "Inh L5-6 SST BEAN1", "Inh L5-6 SST DNAJC14",
                            "Inh L5-6 SST C4orf26", "Inh L5-6 SST FBN2", "Inh L5-6 SST KLHL1",
                            "Inh L6 SST TH", "Inh L5-6 PVALB KCNIP2", "Inh L5-6 PVALB ZFPM2-AS1",
                            "Inh L3-5 PVALB ISG20", "Inh L5 PVALB LRIG3", "Inh L2-5 PVALB HHIPL1",
                            "Inh L2-5 PVALB RPH3AL", "Inh L3 PVALB SAMD13", "Inh L1-2 PVALB CDK20",
                            "Inh L2 PVALB FRZB", "Inh L1-2 SST CLIC6", "Inh L5-6 PVALB GAPDHP60",
                            "Inh L5-6 PVALB FAM150B", "Inh L5-6 PVALB MEPE", "Inh L1-6 PVALB COL15A1",
                            "Exc L2 LAMP5 KCNG3", "Exc L2 LINC00507 ATP7B", "Exc L2 LINC00507 GLRA3",
                            "Exc L2-3 RORB RTKN2", "Exc L2-3 LINC00507 DSG3", "Exc L3 LAMP5 CARM1P1",
                            "Exc L2-3 RORB CCDC68", "Exc L2-3 RORB PTPN3",
                            "Exc L3 THEMIS ENPEP", "Exc L3-5 RORB TNNT2", "Exc L3-5 RORB LAMA4",
                            "Exc L3 RORB OTOGL", "Exc L5 THEMIS VILL", "Exc L3-5 RORB LINC01202",
                            "Exc L5 THEMIS SLC22A18", "Exc L5-6 THEMIS TNFAIP6", "Exc L3-5 RORB LNX2",
                            "Exc L3-5 RORB RPRM", "Exc L5 RORB MED8", "Exc L5 THEMIS FGF10",
                            "Exc L6 THEMIS LINC00343", "Exc L6 THEMIS SLN", "Exc L6 THEMIS SNTG2",
                            "Exc L5-6 THEMIS SMYD1", "Exc L5-6 FEZF2 OR1L8", "Exc L5 THEMIS LINC01116",
                            "Exc L5 THEMIS RGPD6", "Exc L5-6 FEZF2 C9orf135-AS1", "Exc L5-6 FEZF2 FILIP1L",
                            "Exc L5-6 FEZF2 SH2D1B", "Exc L6 FEZF2 PDYN", "Exc L6 FEZF2 FFAR4",
                            "Exc L6 FEZF2 PROKR2", "Exc L5-6 FEZF2 CFTR", "Exc L6 FEZF2 KLK7",
                            "Exc L6 FEZF2 POGK", "Exc L5 FEZF2 CSN1S1", "Exc L3-5 FEZF2 ASGR2",
                            "Exc L3-5 FEZF2 LINC01107", "Exc L5 FEZF2 PKD2L1", "Exc L5 FEZF2 NREP-AS1",
                            "Exc L5 FEZF2 RNF144A-AS1", "Exc L5-6 FEZF2 IFNG-AS1", "Exc L5-6 FEZF2 LPO")


MARKER_1_COLORS <- c("SST"="#e5c185", "PVALB"="#b8cdab", "VIP"="#008585", "PAX6"="#004343", "LAMP5"="#c7522a",
                    "LINC00507"="#e5c185", "RORB"="#004343", "THEMIS"="#008585", "FEZF2"="#c7522a")

LAYERS_COLORS <- c("L1"="#1B9E77", "L1-2"="#8DD3C7", "L1-3"="#66A61E", "L1-5"="#0072B2", "L1-6"="#80B1D3",
                   "L2"="#E6AB02",  "L2-3"="#FFFFB3", "L2-5"="#FDB462",
                   "L3"="#E7298A",   "L3-5"="#BEBADA", "L3-6"="#7570B3",
                   "L5"="#D95F02",   "L5-6"="#FB8072",
                   "L6"="#666666")


(function() {
  test.col <- colorRampPalette(c("white", "#4db6ac", "#ef6c00"))(21)

  df_m <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2m"]], "/abundance.allen"))
  df_p <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2p"]], "/abundance.allen"))

  df_p <- df_p %>% mutate(topic=paste0("EXC L2-3 p", topic), cux2='p')
  df_m <- df_m %>% mutate(topic=paste0("EXC L4-6 p", topic), cux2='m')

  df_exc <- rbind(df_m, df_p)
  fractions <- df_exc %>% distinct(label, cux2, n) %>% group_by(label) %>% mutate(frac=n/sum(n)) %>% ungroup %>% arrange(label)

  df_exc <- df_exc %>% dplyr::left_join(fractions, by=c("label", "cux2", "n")) %>% mutate(avg=avg*frac)

  df_inhibitory <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["inhibitory"]], "/abundance.allen"))
  df_inhibitory <- df_inhibitory %>% mutate(topic=paste0("INH p", topic))

  neuronal_class_to_df <- list("inhibitory"=df_inhibitory, "excitatory"=df_exc)

  for (neuronal_class in c("inhibitory", "excitatory")) {
    df <- neuronal_class_to_df[[neuronal_class]]

    neuronal.subtypes <- df %>%
      dplyr::select(label) %>%
      unique %>%
      tidyr::separate(col = "label", sep = " ", into = c("type", "layer", "marker.1", "marker.2"), extra = "merge", remove = FALSE) %>%
      mutate(marker.1 = factor(marker.1, rev(c("LAMP5", "PAX6", "SST", "VIP", "PVALB","RORB", "LINC00507","THEMIS","FEZF2"))),
             across(c(layer, marker.2), ~factor(., sort(unique(.))))) %>%
      remove_rownames %>%
      column_to_rownames("label")

    df <- df  %>%  pivot_wider(id_cols="label", names_from = "topic", values_from = "avg", values_fill = 0) %>%
      column_to_rownames("label") %>%
      as.matrix %>%
      t

    inhibitory_programs_order <- paste0("INH p", c(14, 11, 10, 8, 7, 1, 17, 6, 2, 9, 16, 5, 13, 15, 3, 4, 12))
    cux2p_programs_order <- paste0("EXC L2-3 p", c(14, 11, 10, 9, 4, 1, 3, 7, 6, 8, 12, 2, 5, 13))
    cux2m_programs_order <- paste0("EXC L4-6 p", c(), c(6,16, 2,13, 14,11, 12, 15, 1,10,
                                                        3,4,5,7,8,9))

    programs_order <- c(inhibitory_programs_order,
                        cux2m_programs_order,
                        cux2p_programs_order)

    marker_1_colors <- MARKER_1_COLORS[unfactor(unique(neuronal.subtypes$marker.1))]

    withr::with_pdf(file.path(figures_dir, glue::glue("4SupD_E_allen_label_{neuronal_class}.pdf")),
                    height=embed.height*1.5, width = embed.width*3, {
          Heatmap(df,
                  col=test.col,
                  height = unit(nrow(df)*.35,"cm"),
                  # width = unit(ncol(df)*.35,"cm"),
                  border=T,
                  cluster_row_slices = FALSE,
                  show_column_dend = FALSE,
                  show_row_dend = FALSE,
                  column_order = NEURONAL_SUBTYPE_ORDER[NEURONAL_SUBTYPE_ORDER %in% colnames(df)],
                  cluster_rows = F,
                  row_order = programs_order[programs_order %in% rownames(df)],
                  bottom_annotation = columnAnnotation(
                    `Cortical layer` = neuronal.subtypes[colnames(df), "layer"],
                    `Marker 1` = neuronal.subtypes[colnames(df), "marker.1"],
                    `Cortical layer ` = anno_text(paste(neuronal.subtypes[colnames(df), "layer"], ""), just = "left", location=0),
                    `Marker 1 ` = anno_text(paste(" ", neuronal.subtypes[colnames(df), "marker.1"], "  ")),
                    `Marker 2 ` = anno_text(neuronal.subtypes[colnames(df), "marker.2"]),
                    col = list(
                      `Marker 1`=marker_1_colors,
                      `Cortical layer`=LAYERS_COLORS)
                  ),
                  left_annotation=rowAnnotation(
                    row_names = anno_text(
                      rownames(df), just = "left", location = 0
                    )
                  ),
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  name="Mean program score") %>% draw()
    })
  }

  for (ct in c("inhibitory", "cux2p", "cux2m")) {
    for (agg in c("marker.1", "class")) {
      df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], glue::glue("/abundance.allen.{agg}")))
      df <- df %>% filter(n > 1000)
      df <- pivot_wider(df, id_cols=agg, names_from = "topic", values_from = "avg", values_fill = 0) %>%
        column_to_rownames(agg) %>%
        as.matrix

      order <- slanter::slanted_orders(df)

      withr::with_pdf(file.path(figures_dir, glue::glue("allen_{agg}_{ct}.pdf")), {
        Heatmap(df,
                col=test.col,
                width = unit(ncol(df)*.35,"cm"),
                border=T,
                cluster_rows = FALSE,
                cluster_row_slices = FALSE,
                column_order = order$cols,
                row_order = order$rows,
                show_column_dend = FALSE,
                show_row_dend = FALSE,
                # column_order = programs.order[[ct]],
                show_row_names = TRUE, name = "Mean program score") %>% draw()
      })
    }
  }
})()


## Figure 1g - trait association --------------

(function() {
  programs_trait_analsysis <- list(
    "main" = data$uns$trait.analysis,
    "batch_mixed_model" = data$uns$trait.analysis_batch_mixed_model
  )

  trait_analysis_to_file <- list(main="1g_trait_assocation_test.pdf",
                                 batch_mixed_model="4SuplF_trait_assocation_batch_mixed_model.pdf")

  for (program_trait_analysis_name in names(programs_trait_analsysis)) {
    program_trait_analysis <- programs_trait_analsysis[[program_trait_analysis_name]]

    df <- program_trait_analysis %>%
      py_to_r %>%
      mutate(trait=AD_TRAITS_NAMES[trait])

    params <- c("Neocortical amyloid",
                "Neocortical tangles",
                "Cognitive decline")

    df$sig <- cut(df$adj.pval, c(-.1, .001, .01, .05, Inf), c("***", "**", "*", ""))

    neuornal_topics <- unlist(lapply(neuronal_topics_to_annotation, names), use.names = FALSE)

    df <- df %>% group_by_at(c("covariate")) %>%
      filter(if_any("adj.pval", ~sum(. < .05) > 0) | covariate %in% neuornal_topics) %>%
      ungroup() %>%
      data.frame(check.names=FALSE)

    df <- df %>% mutate(covariate=topic_names_mapping(covariate, use_static_mapping=FALSE))

    vals <- tidyr::pivot_wider(df, id_cols = covariate, names_from = trait, values_from = 'tstat', values_fill = NA_real_) %>%
      tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))
    sig <- tidyr::pivot_wider(df, id_cols =covariate, names_from = trait, values_from = "sig", values_fill = "") %>%
      tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))

    vals <- vals %>% as.matrix()
    sig <- sig %>% as.matrix()

    v <- max(abs(vals), na.rm = T)
    col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])

    annotation <- rowAnnotation(`cell type`=rownames(vals) %>% str_extract("\\w+"),
                                    col = list(`cell type`=c(
                                      "AST" = "#004343",
                                      "EXC" = "#f0daa5",
                                      "INH" = "#c7522a",
                                      "MIC" = "#008585",
                                      "OLI" = "#74a892",
                                      "OPC" = "#586f44"
                                    )))

    hm <- Heatmap(vals,
                  name = "t",
                  col = circlize::colorRamp2(col.vals, window2fox(length(col.vals))),
                  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                  cluster_columns = FALSE,
                  right_annotation = annotation,
                  clustering_distance_columns = "euclidean",
                  column_names_rot = 45,
                  show_column_dend = FALSE,
                  show_row_dend = FALSE)

    legend <- Legend(title = "FDR", pch = c("***","**","*"),
                    type = "points", labels = c("<0.001","<0.01", "<0.05"))

    withr::with_pdf(file.path(figures_dir, trait_analysis_to_file[[program_trait_analysis_name]]), width=embed.width * 0.7, height = embed.height * 2, {
        draw(hm, annotation_legend_list = list(legend), merge_legend=T, column_title=program_trait_analysis_name)
    })
  }
 })()


#+++++++++++++++++++++++++++
# Figure 2 -----------------
#+++++++++++++++++++++++++++

#+++++++++++++++++++++++++++
## Figure 2a neuronal topic correlation --------------
#+++++++++++++++++++++++++++

(function() {
  topics <- unlist(unname(neuronal_topics_to_annotation))
  df <- data$layers[['sqrt.prev']][,names(topics)] %>% as.data.frame() %>% dplyr::rename(!!!setNames(names(topics), topics))

  c <- correlate_dataframes(df ,df)

  add_significance <- TRUE
  if (add_significance) {
    corr.sig <- c$sig
    corr.sig[corr.sig == '***'] <- '**'
    cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y, gp = grid::gpar(fontsize=6, fontface="bold"))
  } else {
    cell_fun <- NULL
  }


  legend <- Legend(title = "adj.pval", pch = c("*", "**"), type = "points", labels = c("<0.05", "<0.01"))
  correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation)

  withr::with_pdf(file.path(figures_dir, glue::glue("2a_neuronal_correlation.pdf")), {
    ComplexHeatmap::Heatmap(c$corr,
                            column_title_gp = gpar(fontsize=10, fontface="bold"),
                            heatmap_legend_param = list(title="Programs correlation across individual", direction="vertical"),
                            col = cyan2orange.correlation,
                            clustering_distance_rows = "spearman",
                            clustering_distance_columns = "spearman",
                            row_labels = topic_names_mapping(colnames(df), use_static_mapping=FALSE),
                            column_labels = topic_names_mapping(colnames(df), use_static_mapping=FALSE),
                            cluster_rows=FALSE,
                            cluster_columns=FALSE,
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            column_names_rot = 90,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            width=unit(length(topics) * 7, "mm"),
                            height=unit(length(topics) * 7, "mm"),
                            cell_fun = cell_fun, show_heatmap_legend =FALSE) %>%
                            draw(annotation_legend_list=list(legend, correlation.legend),
                            column_title=glue::glue("Programs correlation across individual"))
  })

})

#+++++++++++++++++++++++++++
## Figure 2a neuronal topic jaccard --------------
#+++++++++++++++++++++++++++

(function () {
  topics <- unlist(lapply(neuronal_topics_to_annotation, names), use.names = FALSE)
  des <- read.dataset.for.topics(topics, dataset_path = "/de")

  des_nested <- des %>%
    group_by(topic_name) %>%
    dplyr::select(gene, topic_name) %>%
    group_nest() %>%
    arrange(match(topic_name, topics))

  mat <- outer(des_nested %>% pull(data),
               des_nested %>% pull(data),
               Vectorize(tibble_jaccard))

  rownames(mat) <- des_nested %>% pull(topic_name)
  colnames(mat) <- des_nested %>% pull(topic_name)

  withr::with_pdf(file.path(figures_dir,
                            glue::glue("2a_neuronal_jaccard.pdf")), {
    Heatmap(mat,
            cluster_rows = FALSE, cluster_columns = FALSE,
            name="jaccard",
            row_labels = topic_names_mapping(rownames(mat), use_static_mapping=FALSE),
            column_labels = topic_names_mapping(colnames(mat), use_static_mapping=FALSE),
            # cell_fun = function(j, i, x, y, width, height, fill) {
            #   grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))},
            col=circlize::colorRamp2(c(0, 1), c("white", "#ef6c00")),
            show_column_dend = FALSE, show_row_dend = FALSE,) %>%
                                draw(column_title=glue::glue("Programs signature jaccard index"))
  })

})()


## Figure 2b - phate probability --------

(function() {
  df_branch <- data$uns$trajectories$palantir$branch.probs %>% py_to_r %>% dplyr::select(prAD, ABA)
  pallet <- colorRampPalette(c(
    "#56B4E9", "#A3C8DA", "lightgray", "#db836b", "#c7522a", "#ad4623"
  ), bias=0.4)

  withr::with_pdf(file.path(figures_dir, "2b_phate.pdf"), width=embed.width, height=embed.height, {
    plot.landscape(df_branch %>% mutate(diff=prAD-ABA) %>% dplyr::select(diff),
                   smoothened = F, size = 3, cols = pallet, 
                   embedding='X_phate', enforce.same.color.scale = FALSE, legend.position = 'bottom'
    ) %>% print()
  })


})()


## Figure 2c - trait dynamics --------

AD.traits        <- c(sqrt.amyloid_mf="neocortical amyloid", sqrt.tangles_mf="neocortical tangles", cogng_demog_slope="cognitive decline rate")
AD.traits.colors <- c(sqrt.amyloid_mf="midnightblue", sqrt.tangles_mf="olivedrab4", cogng_demog_slope="firebrick3")

(function() {
  withr::with_pdf(file.path(figures_dir, glue::glue("2c_trait_dynamics.pdf")), width=embed.width, height=embed.height*1.75, {
    plot_grid(plotlist = lapply(names(AD.traits), function(trait)
      plot.dynamics.wrapper(dynamics = data$uns$trajectories$palantir$dynamics,
                            features=c(trait),
                            cols=AD.traits.colors,
                            overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME, ncol=2, strip.position="top", scales="free",
                            label = F, legend.position = c(2.5,0), include.points=F) +
        labs(x="Pseudotime", y=AD.traits[[trait]], title=NULL)),
      ncol=1) %>% print()
  })
})()


## Extended Data Figure 5c - trait dynamics with points --------

(function() {
  withr::with_pdf(file.path(figures_dir, glue::glue("5SuplC_trait_dynamics_dots.pdf")), width=embed.width, height=embed.height*1.75, {
    plot_grid(plotlist = lapply(names(AD.traits), function(trait)
      plot.dynamics.wrapper(dynamics = data$uns$trajectories$palantir$dynamics,
                            features=c(trait),
                            cols=AD.traits.colors,
                            overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                            ncol=2, strip.position="top", scales="free",
                            label = F, legend.position = c(2.5,0), include.points=T) +
        labs(x="Pseudotime", y=AD.traits[[trait]], title=NULL)),
      ncol=1) %>% print()
  })
})()

## Figure 2d - Glia dynamics --------
(function() {
  dynamics_options <- list('paper'= list(panel='2d', dynamics=data$uns$trajectories$palantir$dynamics,
                                      overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                      xlab='Pseudotime'),
                           'green et al'=list(panel='5SuplD', dynamics=data$uns$trajectories$da_fits_sqrt,
                                           overlap.pseudotime=GREEN_ET_AL_OVERLLAPING_PSEUDOTIME,
                                           xlab='Green et al pseudotime'))

  groups <- list("Astrocytes"=c("Ast.10.k6"="AST.p6", "Ast.10.k5"="AST.p5"),
                 "Oligodendrocytes"=c("Oligo.8.k8"="OLI.p8", "Oligo.8.k6"="OLI.p6"),
                 "Microglia"=c("Micro.15.k8"="MIC.p8", "Micro.15.k9"="MIC.p9"))

  for (option in dynamics_options) {
    for (group_name in names(groups)) {
      withr::with_pdf(file.path(figures_dir, glue::glue("{option$panel}_{group_name}.pdf")), width=embed.width, height=embed.height.small, {
        print(plot.dynamics.wrapper(option$dynamics,
                                    features = names(groups[[group_name]]),
                                    labels=groups[[group_name]],
                                    cols=c("#ef6c00", "#4db6ac"),
                                    overlap.pseudotime=option$overlap.pseudotime,
                                    ncol=2, strip.position="top", scales="free_x",
                                    ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=TRUE, include.points=FALSE, show.label.legend=T) +
                theme(
                  legend.position="none",
                  axis.line = element_line(),
                ) +
                labs(x=option$xlab, y="sqrt(Program abundance)", title=""))
      })
    }
  }
})()


## Figure 2e Neuronal dynamics--------------

(function () {
  dynamics_options <- list('paper'= list(panel='2e', dynamics=data$uns$trajectories$palantir$dynamics,
                                    overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                    xlab='Pseudotime'),
                         'green et al'=list(panel='5SuplE', dynamics=data$uns$trajectories$da_fits_sqrt,
                                         overlap.pseudotime=GREEN_ET_AL_OVERLLAPING_PSEUDOTIME,
                                         xlab='Green et al pseudotime'))

  group_names <- c("decreasing_early", "decreasing_late", "increasing_early", "increasing_late")
  for (option in dynamics_options) {
    for (group_name in group_names) {
      features <- names(neuronal_topics_to_annotation[[group_name]])
      labels <- topic_names_mapping(features, use_static_mapping=FALSE)
      names(labels) <- features
      withr::with_pdf(file.path(figures_dir, glue::glue("{option$panel}_{group_name}_dynamics.pdf")), {
        print(plot.dynamics.wrapper(option$dynamics,
                                    features = features,
                                    labels=labels,
                                    cols=neuronal_topic_to_color(features),
                                    overlap.pseudotime=option$overlap.pseudotime,
                                    ncol=1, strip.position="left", scales="free_x", ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=TRUE, include.points=FALSE, trajectories='prAD', show.label.legend=T) +
          theme(strip.text = element_blank(),
                legend.position="none",
                axis.line = element_line(),
                axis.text.y = element_text()) +
          labs(x=option$xlab , y="sqrt(Program abundance)", title=""))

        print(plot.dynamics.wrapper(option$dynamics,
                                    features = features,
                                    labels=labels,
                                    cols=neuronal_topic_to_color(features),
                                    overlap.pseudotime=option$overlap.pseudotime,
                                    ncol=1, strip.position="left", scales="free_x", ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=FALSE, include.points=FALSE, trajectories='prAD', show.label.legend=T) +
          theme(strip.text = element_blank(),
                legend.position="bottom",
                axis.line = element_line(),
                axis.text.y = element_text()) +
          labs(x=option$xlab , y="sqrt(Program abundance)", title=""))

        })
    }
  }
})()


#+++++++++++++++++++++++++++
# Extended Data Figure 5  -----------------
#+++++++++++++++++++++++++++


# Sup Figure 5a branch probability correlation -----------------
(function () {
  data.da <- read_h5ad(AD_500_ANDATA_FILE)

  c <- correlate_dataframes(data$uns$trajectories$palantir$branch.probs %>% py_to_r(),
                            data.da$uns$trajectories$palantir$branch.probs %>% py_to_r())


  add_significance <- TRUE
  if (add_significance) {
    corr.sig <- c$sig
    corr.sig[corr.sig == '***'] <- '**'
    cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y, gp = grid::gpar(fontsize=6, fontface="bold"))
  } else {
    cell_fun <- NULL
  }

  legend <- Legend(title = "adj.pval", pch = c("*", "**"), type = "points", labels = c("<0.05", "<0.01"))
  correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation)

  withr::with_pdf(file.path(figures_dir, glue::glue("5SuplA_branch_probability_cor.pdf")),
                  height = embed.height.small,
                  width = embed.width.small,
                  {
    ComplexHeatmap::Heatmap(c$corr,
                            col = cyan2orange.correlation,
                            cluster_rows=FALSE,
                            cluster_columns=FALSE,
                            show_row_dend = FALSE,
                            row_order = c('prAD', 'ABA'),
                            column_order = c('prAD', 'ABA'),
                            show_column_dend = FALSE,
                            column_names_rot = 90,
                            show_row_names = TRUE,
                            column_title = 'Green et al',
                            column_title_side = 'bottom',
                            show_column_names = TRUE,
                            cell_fun = cell_fun,
                            show_heatmap_legend =FALSE) %>%
      draw(annotation_legend_list=list(legend, correlation.legend),
           column_title=glue::glue("Branch probability correlation"))
  })
})()

# Sup Figure 5b pseudotime comparison -----------------

(function () {
  data.da <- read_h5ad(AD_500_ANDATA_FILE)

  cor(data.da$uns$trajectories$palantir$pseudotime,
      data$uns$trajectories$palantir$pseudotime,
      method = 'spearman', use = 'pairwise.complete.obs')

  df <- data.frame(green_et_al=data.da$uns$trajectories$palantir$pseudotime,
                   green_et_al=data.da$uns$trajectories$palantir$branch.probs %>% py_to_r,
                   our=data$uns$trajectories$palantir$pseudotime,
                   our=data$uns$trajectories$palantir$branch.probs %>% py_to_r)
  
  pallet <- colorRampPalette(c(
    "#56B4E9", "#A3C8DA", "lightgray", "#db836b", "#ad4623"
  ), bias=0.4)
  
  withr::with_pdf(file.path(figures_dir, glue::glue("5SuplB_pseudotime_comparison.pdf")),
                  height = embed.height,
                  width = embed.width, {
      (ggplot(df, aes(x=green_et_al, y=our, color=our.prAD - our.ABA)) +
         geom_point() +
         scale_color_gradientn(colors=pallet(20)) +
         theme_minimal() +
         labs(x="Green et al pseudotime", y='Pseudotime', color='prAD branch prob.', title = "Pseudotime comparision")
      )  %>% print
  })

})()



#+++++++++++++++++++++++++++
# Figure 3 -----------------
#+++++++++++++++++++++++++++


## Figure 3.b pathway bubble plot --------------

(function () {
  inhibitory_pts <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["inhibitory"]], "/pathways"), read.attributes = TRUE)


  pathways_to_show <- list(k10=c("Metallothioneins bind metals", "actin binding", "CaM pathway", "microtubule", #?
                                 "synaptic vesicle", "Synaptic vesicle cycle", "Translation",
                                 "oxidative phosphorylation", "regulation of cell morphogenesis",
                                 "response to oxidative stress", "neuron apoptotic process", "regulation of transporter activity",
                                 "rRNA binding", "neuron spine"),
                           k8=c("M Phase", "Centrosome maturation", "DNA Repair", "chromatin remodeling", "sister chromatid segregation", "RNA degradation",
                                "RNA splicing"))

  group_names <- list(k10=c('1'="response to metal ions",
                              '2'="actin binding",
                              '3'="ca-dependent events,synaptic vesicle",
                              '4'="apoptotic process",
                              '5'="translation",
                              '6'="oxidative phosphorylation",
                              '7'="cell morphogensis",
                              '8'="oxidative stress",
                              '9'="ion channel",
                              '10'="ribosome proteins"),
                      k8=c('1'='centrosome',
                           '2'="DNA repair,chromatin remodeling",
                           '3'='M Phase',
                           '4'='sister chromatid segregation',
                           '5'='RNA degradation',
                           '6'='response to UV',
                           '7'='RNA splicing'))

  withr::with_pdf(file.path(figures_dir, "3b_inhibitory_k8_pathway_network_rast_edges_only.pdf"), {
    print(plot_pathways_overview_graph(inhibitory_pts, k='k8', max_distance = 0.7,
                                       rasterise_edges = TRUE,
                                       group_names=group_names$k8,
                                       seed=2, named_nodes=pathways_to_show[["k8"]]))
  })

  withr::with_pdf(file.path(figures_dir, "3b_inhibitory_k10_pathway_network_rast_edges_only.pdf"), width = embed.width * 2,
                  height = embed.height * 2, {
      print(plot_pathways_overview_graph(inhibitory_pts, k='k10', max_distance = 0.3,
                                         rasterise_edges = TRUE,
                                         group_names=group_names$k10,
                                         seed=13, named_nodes=pathways_to_show[["k10"]]))
  })

})()



## Extended Data Figure 6c - antisene per program  --------------

(function () {
  for (ct in c("inhibitory", "cux2p", "cux2m")) {
    des <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/de"), read.attributes = TRUE)
    des <- des %>% mutate(is_as=grepl("-AS\\d+", gene))

    universe <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/F"), read.attributes = TRUE) %>% pull(gene)

    number_of_universe_genes <- length(universe)
    number_of_universe_as <- sum(grepl('-AS\\d+$', universe))
    as_stats <- des %>% group_by(topic) %>% summarise(number_of_de=n(), number_of_as_de=sum(is_as))

    as_stats <-as_stats %>%  mutate(p.value= 1 - stats::phyper(number_of_as_de, number_of_universe_as,
                                                               number_of_universe_genes - number_of_universe_as,
                                                               number_of_de)) %>%
      mutate(p.adjust=p.adjust(p.value, method = "BH"))

    wider_as_stats <- as_stats %>%
      mutate(number_of_not_as=number_of_de - number_of_as_de) %>%
      pivot_longer(cols=c(number_of_not_as, number_of_as_de), names_to = 'de_type', values_to = 'n') %>%
      mutate(topic=factor(str_replace(topic, "k", "p"), levels=paste0('p', 1:length(topic))))

    ct_name <- case_when(ct == 'inhibitory' ~ "Inhibitory",
                         ct == 'cux2p' ~ "Exc L2-3",
                         ct == 'cux2m' ~ "Exc L4-6")

    withr::with_pdf(file.path(figures_dir, glue::glue("6SupC_anti_sense_{ct}.pdf")), width = embed.width * 1.3, height = embed.height, {
      print(
        ggplot(wider_as_stats, aes(x = topic, y = n, fill=de_type)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_manual(values = c("number_of_not_as" = "#005F60", "number_of_as_de" = "lightblue"),
                            labels = c("Non-antisense", "Antisense"),
                            name="Gene type") +
          annotate("text", x = wider_as_stats$topic, y = wider_as_stats$number_of_de + 30,
                   label = ifelse(wider_as_stats$p.adjust < 0.001, "***","") , size = 6) +
          labs(title=glue::glue("Number of de anti-sense genes - {ct_name}"),
               x="Programs",
               y="Number of genes") +
          theme_classic()
      )
    })
  }
})


## Figure 2f, Extended data figure 6g,i --------------


(function(){
  pathways_to_name_and_topic_group <- list(
    "GO:0005761" = c("mitochondrial ribosome", "decreasing_early"),
    "hsa04721" = c("Synaptic vesicle cycle", "increasing_early"),
    "GO:0051402" = c("neuron apoptotic process", "increasing_early"),
    "GO:0006979" = c("response to oxidative stress", "increasing_early"),
    "R-HSA-68886" = c("M Phase", "increasing_late"),
    "GO:0071103" = c("DNA conformation change", "increasing_late"),
    "R-HSA-73894" = c("DNA Repair", "increasing_late"),
    "GO:0008380" = c("RNA splicing", "increasing_late"),
    "GO:0051028" = c("mRNA transport", "decreasing_late"),
    "GO:0099572" = c("postsynaptic specialization", "decreasing_late"), # excitatory only
    "R-HSA-3371453" = c("Regulation of HSF1-mediated heat shock response", "increasing_late"),
    "R-HSA-5205647" = c("Mitophagy", "decreasing_early")
  )

  for (pathway_id in names(pathways_to_name_and_topic_group)) {
    group_name <- pathways_to_name_and_topic_group[[pathway_id]][[2]]
    pathway_name <- pathways_to_name_and_topic_group[[pathway_id]][[1]]

    message(pathway_name)

    group <- neuronal_topics_to_annotation[[group_name]]
    de <- read.dataset.for.topics(names(group), dataset_path = "/de")
    pathway.genes <- get.pathway.genes(pathway_id)
    de <- de %>% filter(gene %in% pathway.genes)
    topic.to.genes <- sapply(de %>% pull(cell.type) %>% unique, function(ct) {
      de %>% filter(cell.type == ct) %>% pull(gene)
    }, simplify = FALSE)

    all_genes <- de$gene %>% unique

    # Mark significant cell types with *
    pt <- read.dataset.for.topics(names(group), dataset_path = "/pathways")
    significant.cell.type <- pt %>% filter(ID == pathway_id) %>% pull(cell.type)
    cell.types <- names(topic.to.genes)
    cell.type.labels <- case_when(cell.types %in% significant.cell.type ~ paste0(cell.types, " (*)"), TRUE ~ cell.types)
    cell.type.labels <- cell.type.labels %>% gsub("cux2m", "Exc L4-6", .) %>% gsub("cux2p", "Exc L2-3", .)
    names(cell.type.labels) <- cell.types

    dynamics <- pathways.dynamics(data=data, genes=all_genes,
                                  pseudobulks=pseudobulks, scale=T, signature=T,
                                  cell.types=cell.types)

    withr::with_pdf(file.path(figures_dir, glue::glue("2f_{group_name}_{pathway_name}.pdf")), width=embed.width, height=embed.height * 1.5, {
      print(plot.dynamics.wrapper(dynamics,
                                  features = cell.types,
                                  cols=neuronal_topic_to_color(cell.types),
                                  labels=cell.type.labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=1, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=T, include.points=FALSE, trajectories='prAD',
                                  show.label.legend=T, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.position="none",
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="Mean pathway expresssion", title=glue::glue("{pathway_id} - {pathway_name} dynamics")))


      print(plot.dynamics.wrapper(dynamics,
                                  features = cell.types,
                                  cols=neuronal_topic_to_color(cell.types),
                                  labels=cell.type.labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=2, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=T, include.points=FALSE,
                                  show.label.legend=T, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.position="none",
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="Mean pathway expresssion", title=glue::glue("{pathway_id} - {pathway_name} dynamics")))

      for (ct in names(topic.to.genes)) {
        print(fit.plot.genes.dynamics.heatmap(data=data,
                                              pseudobulks=pseudobulks,
                                              ct=ct,
                                              genes=all_genes,
                                              title=glue::glue("{pathway_name} - {cell.type.labels[[ct]]}")))
      }
    })
  }

})()

# validation across snuc datasets --------------

VALIDATION_PROJECTS <- c("mit", "cuimc2", "psychAD", "seaad")
PROGRAM_GROUP_TO_PATHWAYS <- list(
    `decreasing early` = c("mitochondrial ribosome", "antinsense"),
    `increasing early` = c("postsynaptic specialization", "mRNA transport"),
    `cell cycle` = c("M Phase", "DNA conformation change", "DNA Repair", "Regulation of HSF1-mediated heat shock response"),
    `stress/apoptosis` = c("Synaptic vesicle cycle", "neuron apoptotic process", "response to oxidative stress", "Metallothioneins bind metals")
)


## Figure 3c - Coherence summary heatmap across all datasets
(function(){
  # Create mapping from kx to program names based on cell type
  # Based on neuronal_topics_to_annotation structure
  program_name_mapping <- function(gene_set, cell_type, project = NULL) {
    # Map based on neuronal_topics_to_annotation
    # For inhibitory: k1=decreasing_early, k7=decreasing_late, k10=increasing_early, k8=increasing_late
    if (cell_type == "inhibitory") {
      mapping <- c("k1" = "decreasing-early", "k7" = "decreasing-late", 
                   "k8" = "increasing-2", "k10" = "increasing-1", 
                   "decreasing_early" = "decreasing-early", "decreasing_late" = "decreasing-late", 
                   "increasing_late" = "increasing-2", "increasing_early" = "increasing-1")
    } else if (cell_type == "excitatory") {
     mapping <- c("decreasing_early" = "decreasing-early", "decreasing_late" = "decreasing-late", 
                   "increasing_late" = "increasing-2", "increasing_early" = "increasing-1")
    } else {
      return(gene_set)  # Return original if unknown
    }
    
    if (gene_set %in% names(mapping)) {
      return(mapping[gene_set])
    }
    return(gene_set)  # Return original if not in mapping
  }
  
  # Collect coherence results for all datasets
  all_coherence <- list()
  
  for (project in VALIDATION_PROJECTS) {
    for (cell.type in c("excitatory", "inhibitory")) {
      coherence_path <- paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project, "/coherence_results")
      
      tryCatch({
        coherence_results <- h5read(TOPICS_DATA_H5_PATH, coherence_path)
        
        # Fix gene_set if it's stored as bytes (b'k1' -> 'k1')
        if (is.character(coherence_results$gene_set) || is.factor(coherence_results$gene_set)) {
          # Remove 'b' prefix and quotes if present (from bytes encoding)
          coherence_results$gene_set <- gsub("^b['\"]|['\"]$", "", as.character(coherence_results$gene_set))
          coherence_results$gene_set <- gsub("^b", "", coherence_results$gene_set)
        }
        # Convert to character if needed
        coherence_results$gene_set <- as.character(coherence_results$gene_set)
        
        coherence_results$project <- project
        coherence_results$cell_type <- cell.type
        
        # Map gene_set (kx) to program name
        coherence_results$program_name <- sapply(coherence_results$gene_set, function(gs) {
          program_name_mapping(gs, cell.type, project)
        })
        
        # Create column name: program_name + cell_type
        coherence_results$program_celltype <- paste0(coherence_results$program_name, "_", cell.type)
        
        all_coherence[[paste0(project, "_", cell.type)]] <- coherence_results
      }, error = function(e) {
        message(glue::glue("Could not load coherence for {project}_{cell.type}: {e$message}"))
      })
    }
  }
  
  if (length(all_coherence) == 0) {
    message("No coherence results found")
    return()
  }
  
  # Combine all results
  coherence_df <- bind_rows(all_coherence)
  
  # Get all unique projects and program_celltype combinations
  projects <- sort(unique(coherence_df$project))
  
  # Define the order of program_celltypes based on the covariate mapping in figures.R (lines 1866-1886)
  # Order: inhibitory decreasing-early, excitatory decreasing-early, inhibitory decreasing-late,
  #        excitatory decreasing-late, inhibitory increasing-1, excitatory increasing-1,
  #        inhibitory increasing-2, excitatory increasing-2
  program_celltype_order <- c(
    "decreasing-early_inhibitory",
    "decreasing-early_excitatory",
    "decreasing-late_inhibitory",
    "decreasing-late_excitatory",
    "increasing-1_inhibitory",
    "increasing-1_excitatory",
    "increasing-2_inhibitory",
    "increasing-2_excitatory"
  )
  
  # Filter to only those present in the data
  program_celltypes <- program_celltype_order[program_celltype_order %in% unique(coherence_df$program_celltype)]
  
  # Create matrices for effect size and p-values
  # Rows = program_celltype (signatures), Columns = project (datasets)
  effect_size_mat <- matrix(NA, nrow = length(program_celltypes), ncol = length(projects),
                            dimnames = list(program_celltypes, projects))
  p_adj_mat <- matrix(NA, nrow = length(program_celltypes), ncol = length(projects),
                      dimnames = list(program_celltypes, projects))
  
  # Fill matrices
  # Rows = program_celltypes, Columns = projects
  for (i in seq_along(program_celltypes)) {
    for (j in seq_along(projects)) {
      subset_data <- coherence_df %>% 
        filter(program_celltype == program_celltypes[i], project == projects[j])
      
      if (nrow(subset_data) > 0) {
        # Use cohens_d as effect size (or cliffs_delta if cohens_d is NA)
        cohens_d <- if (!is.na(subset_data$cohens_d[1]) && !is.nan(subset_data$cohens_d[1])) {
          subset_data$cohens_d[1]
        } else if (!is.na(subset_data$cliffs_delta[1]) && !is.nan(subset_data$cliffs_delta[1])) {
          subset_data$cliffs_delta[1]
        } else {
          NA
        }
        
        # Use log of Cohen's d (preserve sign)
        # For positive: log(d + 1), for negative: -log(|d| + 1)
        if (!is.na(cohens_d)) {
          if (cohens_d > 0) {
            log_effect_size <- log(cohens_d + 1)
          } else if (cohens_d < 0) {
            log_effect_size <- -log(abs(cohens_d) + 1)
          } else {
            log_effect_size <- 0
          }
        } else {
          log_effect_size <- NA
        }
        
        effect_size_mat[i, j] <- log_effect_size
        p_adj_mat[i, j] <- subset_data$p_adj[1]
      }
    }
  }
  
  # Create asterisk labels based on p-value significance
  # * p < 0.05, ** p < 0.01, *** p < 0.001
  asterisk_mat <- matrix("", nrow = nrow(effect_size_mat), ncol = ncol(effect_size_mat),
                         dimnames = dimnames(effect_size_mat))

  asterisk_mat[!is.na(p_adj_mat) & p_adj_mat < 0.01] <- "**"
  asterisk_mat[!is.na(p_adj_mat) & p_adj_mat >= 0.01 & p_adj_mat < 0.05] <- "*"
  
  # Create cell labels with asterisks
  cell_labels <- matrix("", nrow = nrow(effect_size_mat), ncol = ncol(effect_size_mat),
                        dimnames = dimnames(effect_size_mat))
  for (i in seq_len(nrow(effect_size_mat))) {
    for (j in seq_len(ncol(effect_size_mat))) {
      if (!is.na(effect_size_mat[i, j])) {
        cell_labels[i, j] <- sprintf("%s", asterisk_mat[i, j])
      }
    }
  }
  
  # Create color mapping for log effect size
  # Use a diverging color scale: negative (blue) to positive (red)
  effect_range <- range(effect_size_mat, na.rm = TRUE)
  max_abs <- max(abs(effect_range))
  
  # Color function: blue (negative) -> white (zero) -> red (positive)
  col_fun <- circlize::colorRamp2(
    breaks = seq(-max_abs, max_abs, length.out = 100),
    colors = colorRampPalette(c("#008585", "white", "#c7522a"))(100)
  )
  
  # Create p-value annotation for legend
  p_value_legend <- Legend(
    title = "p-value",
    pch = c("*", "**"),
    type = "points",
    labels = c("p < 0.05", "p < 0.01"),
    title_gp = gpar(fontface = "plain")
  )
  
  # Create heatmap
  plot_file_path <- file.path(figures_dir, "3C_validation_coherence_summary_heatmap.pdf")
  withr::with_pdf(plot_file_path, width = 12, height = 8, {
    ht <- Heatmap(
      effect_size_mat,
      name = "log(Cohen's d)",
      col = col_fun,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (!is.na(effect_size_mat[i, j])) {
          grid.text(cell_labels[i, j], x, y, gp = gpar(fontsize = 8))
        }
      },
      column_title = "Coherence Summary Across Validation Datasets",
      row_title = "Program",
      column_title_side = "top",
      column_names_rot = 45,
      column_names_side = "bottom",
      row_names_side = "left",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      na_col = "gray90",
      heatmap_legend_param = list(
        title = "log(Cohen's d)",
        title_position = "topcenter",
        legend_direction = "horizontal",
        legend_width = unit(4, "cm"),
        at = pretty(effect_range, n = 5)
      )
    )
    
    draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right",
         merge_legend = TRUE)
    
    # Add p-value legend
    draw(p_value_legend, x = unit(0.98, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
  })
  
  message(glue::glue("Coherence summary heatmap saved to {plot_file_path}"))
})()


## Extended Data Figure 7a,b - Pathway correlation heatmap --------
(function(){

  for (project in VALIDATION_PROJECTS) {
    for (cell.type in c("inhibitory", "excitatory")) {
      c <- h5_read_correlation(cell.type=cell.type, group_name=paste0(project, "/pathway_correlation"))

      dimnames(c$cor) <- list(gsub("^Pathway_", "", c$names[[1]]), gsub("^Pathway_", "", c$names[[2]]))
      dimnames(c$cor) <- list(gsub("_", " ", dimnames(c$cor)[[1]]), gsub("_", " ", dimnames(c$cor)[[2]]))

      # fix HSF1 mediated to HSF1-mediated
      dimnames(c$cor) <- list(gsub("HSF1 mediated", "HSF1-mediated", dimnames(c$cor)[[1]]), gsub("HSF1 mediated", "HSF1-mediated", dimnames(c$cor)[[2]]))

      pathway_names <- rownames(c$cor)

      # Order pathways according to program_group_to_pathways
      ordered_pathways <- unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names=FALSE)
      available_pathways <- pathway_names[pathway_names %in% ordered_pathways]

      c$cor <- c$cor[available_pathways, available_pathways]
      
      # Build group annotation (factor for group)
      group_assignment <- setNames(rep(names(PROGRAM_GROUP_TO_PATHWAYS), lengths(PROGRAM_GROUP_TO_PATHWAYS)),
                                   unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names=FALSE))
      annotation_vec <- group_assignment[available_pathways]

      pathway_group_anno <- ComplexHeatmap::rowAnnotation(
        "Group" = factor(annotation_vec, levels = names(PROGRAM_GROUP_TO_PATHWAYS)),
        col = list(Group = c(
          `decreasing early` = "#b8cdab",
          `increasing early` = "#004343",
          `cell cycle` = "#e5c185",
          `stress/apoptosis` = "#c7522a"
        )),
        annotation_name_side = "top"
      )

      withr::with_pdf(file.path(figures_dir, glue::glue("7SupA_validation_{project}_pathway_correlation_{cell.type}.pdf")), {
        Heatmap(c$cor,
                column_title = glue::glue("{cell.type} - Pathway correlation"),
                row_title = glue::glue("{project} - Pathway correlation"),
                col = cyan2orange.correlation,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_order = unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names = FALSE),
                column_order = unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names = FALSE),
                name='correlation',
                left_annotation = pathway_group_anno) %>% draw
      })
    }
  }
})()



# Figure 3c - Validation Signature trait asscocication heatmap --------

(function() {
  score_type <- "signature"
  for (project in c("mit", "cuimc2", "seaad", "psychAD")) {
    inhibitory_trait_association <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["inhibitory"], "/", project, "/trait_association_", score_type))
    excitatory_trait_association <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["excitatory"], "/", project, "/trait_association_", score_type))

    df <- rbind(inhibitory_trait_association %>% mutate(cell.type="inhibitory"), excitatory_trait_association %>% mutate(cell.type="excitatory"))

    if (project == "seaad") {
      trait_to_pretty_name <- c("amyloid"="Cortical amyloid",
                "tangles"="Cortical tangles",
                "Dementia"="Dementia")
      df <- df %>% mutate(trait=trait_to_pretty_name[trait])
      params <- unname(trait_to_pretty_name)
    } else if (project == "psychAD") {
      trait_to_pretty_name <- c("AD"="AD",
                "Braak_4_5_6"="Braak 4-6",
                "HIGH_CERAD"="CERAD 3-4")
      df <- df %>% filter(trait %in% names(trait_to_pretty_name)) %>% mutate(trait=trait_to_pretty_name[trait])
      params <- unname(trait_to_pretty_name)
    } else {
      df <- df %>% mutate(trait=AD_TRAITS_NAMES[trait])
      params <- c("Neocortical amyloid",
                "Neocortical tangles",
                "Cognitive decline")
    }
    
    df$sig <- cut(df$adj.pval, c(-.1, .001, .01, .05, Inf), c("***", "**", "*", ""))

    
    
    df <- df %>% mutate(covariate=glue::glue("{cell.type}_{covariate}") %>% gsub("SignatureUnique_", "", .))
    covariate_mapping <- c(
      "inhibitory_k1" = "Inhibitory decreasing early",
      "excitatory_decreasing_early" = "Excitatory decreasing early",
      "inhibitory_k7" = "Inhibitory decreasing late",
      "excitatory_decreasing_late" = "Excitatory decreasing late",
      "inhibitory_k8" = "Inhibitory increasing-1",
      "excitatory_increasing_early" = "Excitatory increasing-1",
      "inhibitory_k10" = "Inhibitory increasiing-2",
      "excitatory_increasing_late" = "Excitatory increasing-2"
    )

    covariate_levels <- unname(covariate_mapping)

    df <- df %>%
      mutate(
        covariate = ifelse(
          covariate %in% names(covariate_mapping),
          covariate_mapping[covariate],
          covariate
        ),
        covariate = factor(covariate, levels = covariate_levels)
      )

    vals <- tidyr::pivot_wider(df, id_cols = covariate, names_from = trait, values_from = 'tstat', values_fill = NA_real_) %>%
      tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))
    sig <- tidyr::pivot_wider(df, id_cols =covariate, names_from = trait, values_from = "sig", values_fill = "") %>%
      tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))

    vals <- vals %>% as.matrix()
    sig <- sig %>% as.matrix()

    v <- max(abs(vals), na.rm = T)
    col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])

    hm <- Heatmap(vals,
                name = "t",
                col = circlize::colorRamp2(col.vals, window2fox(length(col.vals))),
                cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_order = covariate_levels,
                clustering_distance_columns = "euclidean",
                column_names_rot = 45,
                show_column_dend = FALSE,
                show_row_dend = FALSE)

    legend <- Legend(title = "FDR", pch = c("***","**","*"),
                    type = "points", labels = c("<0.001","<0.01", "<0.05"))

    withr::with_pdf(file.path(figures_dir, glue::glue("3c_validation_{project}_trait_association_{score_type}.pdf")), width=embed.width, height = embed.height * 2, {
        draw(hm, annotation_legend_list = list(legend), merge_legend=T, column_title=glue::glue("{project} - Trait association - {score_type}"))
    })
  }
})()


# Extended Data Figure 7c - Validation pathway trait asscocication heatmap --------

(function() {
  score_type <- "pathway"
  for (project in c("mit", "cuimc2", "seaad", "psychAD")) {
    print(glue::glue("Processing {project} - {score_type}"))
    inhibitory_trait_association <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["inhibitory"], "/", project, "/trait_association_", score_type))
    excitatory_trait_association <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["excitatory"], "/", project, "/trait_association_", score_type))

    df <- rbind(inhibitory_trait_association %>% mutate(cell.type="inhibitory"), excitatory_trait_association %>% mutate(cell.type="excitatory"))

    if (project == "seaad") {
      trait_to_pretty_name <- c("amyloid"="Cortical amyloid",
                "tangles"="Cortical tangles",
                "Dementia"="Dementia")
      df <- df %>% mutate(trait=trait_to_pretty_name[trait])
      params <- unname(trait_to_pretty_name)
    } else if (project == "psychAD") {
      trait_to_pretty_name <- c("AD"="AD",
                "Braak_4_5_6"="Braak 4-6",
                "HIGH_CERAD"="CERAD 3-4")
      df <- df %>% filter(trait %in% names(trait_to_pretty_name)) %>% mutate(trait=trait_to_pretty_name[trait])
      params <- unname(trait_to_pretty_name)
    } else {
      df <- df %>% mutate(trait=AD_TRAITS_NAMES[trait])
      params <- c("Neocortical amyloid",
                "Neocortical tangles",
                "Cognitive decline")
    }
    
    df <- df %>% mutate(covariate=covariate %>%
     gsub("Pathway_", "", .) %>% gsub("_", " ", .) %>% gsub("HSF1 mediated", "HSF1-mediated", .)) %>%
     filter(covariate != 'Cellular response to heat stress')

    df$sig <- cut(df$adj.pval, c(-.1, .001, .01, .05, Inf), c("***", "**", "*", ""))

    group_assignment <- setNames(rep(names(PROGRAM_GROUP_TO_PATHWAYS), lengths(PROGRAM_GROUP_TO_PATHWAYS)),
                                unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names=FALSE))
    
    annotation_vec <- group_assignment[df$covariate %>% unique]
    pathway_group_anno <- ComplexHeatmap::rowAnnotation(
      "Group" = factor(annotation_vec, levels = names(PROGRAM_GROUP_TO_PATHWAYS)),
      col = list(Group = c(
        `decreasing early` = "#b8cdab",
        `increasing early` = "#004343",
        `cell cycle` = "#e5c185",
        `stress/apoptosis` = "#c7522a"
      )),
      annotation_name_side = "top"
    )

    covariate_levels <- unlist(PROGRAM_GROUP_TO_PATHWAYS, use.names = FALSE)

    for (neuronal_class in c("inhibitory", "excitatory")) {
      df_neuronal_class <- df %>% filter(cell.type == neuronal_class)

      vals <- tidyr::pivot_wider(df_neuronal_class, id_cols = covariate, names_from = trait, values_from = 'tstat', values_fill = NA_real_) %>%
        tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))
      sig <- tidyr::pivot_wider(df_neuronal_class, id_cols =covariate, names_from = trait, values_from = "sig", values_fill = "") %>%
        tibble::column_to_rownames('covariate') %>% dplyr::select(all_of(params))

       vals <- vals %>% as.matrix()
      sig <- sig %>% as.matrix()

      v <- max(abs(vals), na.rm = T)
      col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])

      hm <- Heatmap(vals,
                  name = "t",
                  col = circlize::colorRamp2(col.vals, window2fox(length(col.vals))),
                  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  row_order = covariate_levels,
                  clustering_distance_columns = "euclidean",
                  column_names_rot = 45,
                  show_column_dend = FALSE,
                  show_row_dend = FALSE,
                  left_annotation = pathway_group_anno)

      legend <- Legend(title = "FDR", pch = c("***","**","*"),
                      type = "points", labels = c("<0.001","<0.01", "<0.05"))

      withr::with_pdf(file.path(figures_dir, glue::glue("7SuplC_validation_{project}_trait_association_{score_type}_{neuronal_class}.pdf")), width=embed.width, height = embed.height * 2, {
          draw(hm, annotation_legend_list = list(legend), merge_legend=T, column_title=glue::glue("{project} - Trait association - {score_type} - {neuronal_class}"))
      })
    }
  }
})()

plot_by_subtype <- function(
  df_long_by_subtype,
  cols_to_plot,
  comparisons = NULL,
  sample_id_column = "projid",
  subtype_colors = NULL,
  project = NULL
) {
  project <- if (is.null(project)) get("project", envir = parent.frame()) else project
  # Filter out panels where all values are NA or the variance is 0
  df_long_by_subtype <- df_long_by_subtype %>%
    filter(Name %in% cols_to_plot) %>%
    group_by(Name) %>%
    filter(var(Score, na.rm = TRUE) > 0) %>%
    ungroup()

  if (!is.null(df_long_by_subtype$subtype)) {
    subtype_order <- sort(unique(df_long_by_subtype$subtype))
    df_long_by_subtype$subtype <- factor(df_long_by_subtype$subtype, levels = subtype_order)
  }

  if (!is.null(comparisons)) {
    df_long_by_subtype <- df_long_by_subtype %>%
      filter(subtype %in% unlist(comparisons, use.names = FALSE))
  }

  # Order pathways alphabetically
  df_long_by_subtype$Name <- factor(df_long_by_subtype$Name, levels = cols_to_plot)

  p <- ggplot(df_long_by_subtype, aes(x = subtype, y = Score, fill = subtype)) +
    # geom_boxplot(position = position_dodge(preserve = "single"), outlier.shape = NA, width = 0.65) +
    geom_violin(alpha = 1) +
    stat_summary(
      fun = mean, geom = "crossbar", width = 0.6, fatten = 2, color = "black", size = 0.6,
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~Name, scales = "free", ncol = 1) +
    theme_classic() +
    labs(
      title = paste(
        glue::glue("Score by Subtype (n={n_distinct(df_long_by_subtype[[sample_id_column]])})"),
        glue::glue("{project}")
      ),
      x = "Subtype",
      y = "Mean signature Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(size = 10, face = "bold")
    )
  
  # Add option for color for each subtype
  if (!is.null(subtype_colors)) {
    p <- p + scale_fill_manual(values = subtype_colors)
  }

  # Add t.test comparison if needed
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(
      method = "t.test",
      comparisons = comparisons,
      label = "p.signif",
      hide.ns = FALSE
    )
  }
  p
}


## Figure 5c, Extended Data Figure 9b - validation Vulnerability analysis --------------
(function(){
  project_to_donor_metadata_path <- list(
    "mit" = file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv"),
    "cuimc2" = file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv"),
    "psychAD" = file.path(get_result_dir(), "psychAD", "inhibitory", "donor_metadata.csv"),
    "seaad" = file.path(get_result_dir(), "seaad", "inhibitory", "donor_metadata.csv")
  )
  
  for (project in VALIDATION_PROJECTS) {
    donor.meta.data <- read.csv(project_to_donor_metadata_path[[project]], check.names = FALSE)
    
    df_by_subtype <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP["inhibitory"], "/", project, "/mean_scores_by_subtype"))

    # Per project config
    if (project == "psychAD") {
      df_by_subtype$subtype <- gsub("IN_", "", df_by_subtype$subclass)
      ad_only_donor.meta.data <- donor.meta.data %>% filter(AD_status == "Yes")
      sample_id_column <- "donor_id"
    } else if (project == "seaad") {
      df_by_subtype$subtype <- str_to_upper(df_by_subtype$Subclass)
      ad_only_donor.meta.data <- donor.meta.data %>% filter(`Overall AD neuropathological Change` == 'High')
      sample_id_column <- "Donor ID"
    } else if (project %in% c("cuimc2", "mit")) {
      ad_only_donor.meta.data <- donor.meta.data %>% filter(cogdx %in% c(4, 5))
      sample_id_column <- "projid"
      df_by_subtype$subtype <- gsub("Inh ", "", df_by_subtype$subtype)
    } else {
      stop(glue::glue("Project {project} not supported"))
    }
  
    df_by_subtype <- df_by_subtype %>% filter(.data[[sample_id_column]] %in% ad_only_donor.meta.data[[sample_id_column]])
    
    df_by_subtype$SignatureUnique_k8_Normalized <- (df_by_subtype$SignatureUnique_k8 - min(df_by_subtype$SignatureUnique_k8)) / (max(df_by_subtype$SignatureUnique_k8) - min(df_by_subtype$SignatureUnique_k8))
    df_by_subtype$SignatureUnique_k10_Normalized <- (df_by_subtype$SignatureUnique_k10 - min(df_by_subtype$SignatureUnique_k10)) / (max(df_by_subtype$SignatureUnique_k10) - min(df_by_subtype$SignatureUnique_k10))
    df_by_subtype$k8_by_k10 <- (df_by_subtype$SignatureUnique_k8_Normalized + 1) / (df_by_subtype$SignatureUnique_k10_Normalized + 1)

    
    # Rename the columns before pivoting to get the new names in df_long_by_subtype
    col_rename_map <- c(
      "SignatureUnique_k8" = "Cell-cycle reentry",
      "SignatureUnique_k10" = "Stress/apoptosis",
      "k8_by_k10" = "Cell-cycle reentry by Stress/apoptosis"
    )

    df_by_subtype <- df_by_subtype %>%
      dplyr::rename(!!!setNames(names(col_rename_map), col_rename_map))
    cols_to_plot <- unname(col_rename_map)

    df_long_by_subtype <- df_by_subtype %>%
        pivot_longer(
          cols = all_of(cols_to_plot),
          names_to = "Name",
          values_to = "Score",
        ) %>%
        filter(!is.na(subtype))

    comparisons <- list(c("SST", "PVALB"))
    subtype_colors <- c("SST" = "#ef6c00", "PVALB" = "#4db6ac")

    withr::with_pdf(file.path(figures_dir, glue::glue("5c_validation_{project}_vulnerability_by_subtype_AD.pdf")), height = 12 + 2*length(cols_to_plot), {
        print(plot_by_subtype(df_long_by_subtype, cols_to_plot=cols_to_plot, comparisons=comparisons, sample_id_column=sample_id_column, subtype_colors=subtype_colors))
    })
  }
})()


## Figure 3e, Extended data Figures 7d - proteomics validation -------------


(function () {
  protein_groups <- list(
      "3e" = list(
        "GO:0005761" = c("mitochondrial ribosome", "decreasing_early"),
        "GO:0051028" = c("mRNA transport", "decreasing_late"),
        "GO:0099572" = c("postsynaptic specialization", "decreasing_late"),
        "GO:0051402" = list("neuron apoptotic process", "increasing_early", exclude_genes=c("APOE", "GRIN1")),
        "R-HSA-68886" = c("M Phase", "increasing_late")
        
      ),
      "7SupD"=list("R-HSA-73894" = c("DNA Repair", "increasing_late"), "GO:0006979" = c("response to oxidative stress", "increasing_early"),
                   "R-HSA-3371453" = c("Regulation of HSF1-mediated heat shock response", "increasing_late", exclude_genes=c("ATM")))
  )

  for (panel in names(protein_groups)) {
    pathways_to_name_and_topic_group <- protein_groups[[panel]]

    genes_group <- sapply(names(pathways_to_name_and_topic_group), function(pathway_id) {
      group_name <- pathways_to_name_and_topic_group[[pathway_id]][[2]]
      pathway_name <- pathways_to_name_and_topic_group[[pathway_id]][[1]]
      
      if ('exclude_genes' %in% names(pathways_to_name_and_topic_group[[pathway_id]])) {
        
        exclude_genes <- pathways_to_name_and_topic_group[[pathway_id]][['exclude_genes']] 
        print(exclude_genes)
      } else {
        exclude_genes <- c()
      }
       
      message(pathway_name)

      group <- neuronal_topics_to_annotation[[group_name]]
      de <- read.dataset.for.topics(names(group), dataset_path = "/de")
      pathway.genes <- get.pathway.genes(pathway_id)

      de <- de %>% filter(gene %in% pathway.genes & !(gene %in% exclude_genes))
      topic.to.genes <- sapply(de %>% pull(cell.type) %>% unique, function(ct) {
        de %>% filter(cell.type == ct) %>% pull(gene)
      }, simplify = FALSE)

      de$gene %>% unique
    }, simplify = FALSE)

    pathway_names <- lapply(pathways_to_name_and_topic_group, function(g) g[[1]])


    names(genes_group) <- unname(unlist(pathway_names[names(genes_group)]))

    genes <- unlist(genes_group, use.names = FALSE)

    traits_columns <- c( "sqrt_amyloid_mf", "sqrt_tangles_mf", "cogng_random_slope")
    colData(proteom)$sqrt_amyloid_mf <- sqrt(colData(proteom)$amyloid_mf)
    colData(proteom)$sqrt_tangles_mf <- sqrt(colData(proteom)$tangles_mf)


    all_prot_traits <- associate.proteomic.traits(proteom, genes=genes, traits = traits_columns)
    all_prot_traits <- all_prot_traits %>% left_join(stack(genes_group) %>% dplyr::rename(gene = values, group = ind)) %>% arrange(match(group, names(genes_group)))

    row_labels <- all_prot_traits %>% pull(gene, UniProt)
    names(row_labels) <- make.names(names(row_labels))

    hm <- plot.trait.assocation.multiple.groups(all_prot_traits,
                                                row.by = 'UniProt', params = traits_columns,
                                                show.only.significant=FALSE,
                                                row_gap=unit(5, "mm"),
                                                cluster_row_slices=FALSE,
                                                row_labels=all_prot_traits$gene,
                                                column_labels = AD_TRAITS_NAMES,
                                                row_group_mapping=stack(genes_group) %>% pull(ind, values), plot=FALSE
                                                )
    print(panel)
    withr::with_pdf(file.path(figures_dir, glue::glue("{panel}_proteomics_all.pdf")), height = unit(18, 'mm'), width = unit(4, "mm"), {
      draw(hm$hm, annotation_legend_list=hm$legend, merge_legend=T)
    })
  }
})()



## Figure 3d - ELISA to programs assoication --------------


(function() {
  synapse_trait_col_to_name <- list(synap_3cort_vamp='VAMP',
                                    synap_3cort_syntaxin='SYNTAXIN1',
                                    synap_3cort_synaptophys='SYNAPTOPHYSIN',
                                    synap_3cort_stagmin='SYNAPTOTAGMIN',
                                    synap_3cort_snap25='SNAP25',
                                    synap_3cort_sept5='SEPTIN5',
                                    synap_3cort_complex2='COMPLEXIN II',
                                    synap_3cort_complex1='COMPLEXIN I',
                                    synap_3cort_capture4='SNAP25 - VAMP',
                                    synap_3cort_capture3='SNAP25 - SYNTAXIN1',
                                    synap_3cort_capture2='SYNTAXIN1 - SNAP25',
                                    synap_3cort_capture1='SYNTAXIN1 - VAMP',
                                    zcapture_syn_3cort='Total interaction')
  
  topics_names <- unlist(lapply(neuronal_topics_to_annotation, names), use.names = FALSE)
  
  
  withr::with_pdf(file.path(figures_dir, "3d.pdf"), heidht = embed.height * 1.2, width = embed.width *1.2, {
    plot.trait.associations.heatmap(data$uns$syn.trait.analysis %>%
                                      py_to_r %>%
                                      filter(grepl("Inh|Exc", covariate)) %>% 
                                      filter(covariate %in% topics_names) %>% 
                                      dplyr::arrange(match(covariate, topics_names)) %>% 
                                      mutate(covariate = topic_names_mapping(covariate, use_static_mapping=FALSE))
                                    ,
                                    params = synapse_trait_col_to_name %>% enframe %>% dplyr::select(2:1) %>% deframe(),
                                    row.by='covariate',
                                    cluster_rows=FALSE,
                                    show.only.significant = F,
                                    column_gap = unit(2, "mm"),
                                    column_split=factor(c(rep('Protien level',6), rep('Complexin',2), rep('Interaction',5)),
                                                        levels = c('Complexin', 'Protien level', 'Interaction')),
                                    column_title_gp=grid::gpar(fontsize = 12),
                                    column_names_gp=grid::gpar(fontsize = 8),
                                    row_names_gp=grid::gpar(fontsize = 8),
                                    # row_gap=unit(7, "mm"),
                                    cols=cyan2orange)
  })
})()


## Extended Data Figure 4a - Cell cycle checkpoints ----------

(function() {

  cell_cycle_checkpoints_genes_to_topic <- list("ATR"=neuronal_topics_to_annotation$decreasing_late,
                                                "BRCA1"=neuronal_topics_to_annotation$decreasing_early,
                                                "RAD9A"=neuronal_topics_to_annotation$decreasing_late)

  cell.types <- c("inhibitory", 'cux2m', 'cux2p')


  for (gene in names(cell_cycle_checkpoints_genes_to_topic)) {
    print(gene)
    des <- read.dataset.for.topics(names(cell_cycle_checkpoints_genes_to_topic[[gene]]), dataset_path = "/de")

    significant.cell.type <- des %>% filter(.data[['gene']] == .env[['gene']]) %>% pull(cell.type)
    cell.type.labels <- case_when(cell.types %in% significant.cell.type ~ paste0(cell.types, " (*)"), TRUE ~ cell.types)
    cell.type.labels <- cell.type.labels %>% gsub("cux2m", "Exc L4-6", .) %>% gsub("cux2p", "Exc L2-3", .)
    names(cell.type.labels) <- cell.types

    dynamics <- pathways.dynamics(data=data, genes=c(gene),
                                  pseudobulks=pseudobulks, scale=T, signature=T,
                                  cell.types=cell.types)


    withr::with_pdf(file.path(figures_dir, glue::glue("4a_{gene}.pdf")), width=embed.width, height=embed.height, {
      print(plot.dynamics.wrapper(dynamics,
                                  features = cell.types,
                                  cols=neuronal_topic_to_color(cell.types),
                                  labels=cell.type.labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=1, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=T, include.points=FALSE, trajectories='prAD',
                                  show.label.legend=T, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.position="none",
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="expresssion", title=glue::glue("{gene} dynamics")))

        print(plot.dynamics.wrapper(dynamics,
                                    features = cell.types,
                                    cols=neuronal_topic_to_color(cell.types),
                                    labels=cell.type.labels,
                                    overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                    ncol=1, strip.position="left", scales="free_x",
                                    ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=T, include.points=FALSE,
                                    show.label.legend=T, show_n=FALSE) +
                theme(strip.text = element_blank(),
                      legend.position="none",
                      axis.line = element_line(),
                      axis.text.y = element_text()) +
                labs(x="pseudotime" , y="expresssion", title=glue::glue("{gene} dynamics")))
    })
  }
})()


## Extended Data Figures 6 -  Genes dynamics across cell types ----------

(function() {
  gene_groups <- list(
    "mitophagy"=list(genes=c("PINK1","TOMM5","TOMM7"), group=neuronal_topics_to_annotation$decreasing_early),
    "mitochondrial ribosome"=list(genes=c("MRPL51","MRPL46","MRPL54"), group=neuronal_topics_to_annotation$decreasing_early),
    "neuron apoptotic process"=list(genes=c("BNIP3", "HIPK2"), group=neuronal_topics_to_annotation$increasing_early),
    "postsynaptic specialization"=list(genes=c("DLGAP4", "GRIN2D", "SHANK1"), group=neuronal_topics_to_annotation$decreasing_late),
    "mRNA transport"=list(genes=c("NUP133","NSUN2","DDX19A"), group=neuronal_topics_to_annotation$decreasing_late),
    "response to oxidative stress"=list(genes=c("PRDX1", "PRDX2", "GPX3"), group=neuronal_topics_to_annotation$increasing_early),
    "Metallothioneins bind metals"=list(genes=c("MT3", "MT1E",  "MT1G", "MT1F"), group=neuronal_topics_to_annotation$increasing_early),
    "neurotransmitter transporters"=list(genes=c("SLC1A3", "SLC17A7"), group=neuronal_topics_to_annotation$increasing_early),
    "vesicle release cycle"=list(genes=c("VAMP2", "SNAP25", "CPLX1", "CPLX2"), group=neuronal_topics_to_annotation$increasing_early),
    "heat shock response"=list(genes=c("HSPH1", "HSPA12A", "HSPA4L"), group=neuronal_topics_to_annotation$increasing_late),
    "DNA damage repair"=list(genes=c("ERCC5", "ATM", "POLK", "REV1"), group=neuronal_topics_to_annotation$increasing_late),
    "M phase"=list(genes=c("CENPC", "ANAPC16", "CEP290"), group=neuronal_topics_to_annotation$increasing_late),
    "chromatin remodelling"=list(genes=c("SMARCA1", "SMARCAD1", "CHD1"), group=neuronal_topics_to_annotation$increasing_late),
    "RNA-splicing"=list(genes=c("RNPC3", "SRPK1", "SRSF1"), group=neuronal_topics_to_annotation$increasing_late)
  )

  cell.types <- c("inhibitory", 'cux2m', 'cux2p')


  for (group_name in names(gene_groups)) {
    print(group_name)
    des <- read.dataset.for.topics(names(gene_groups[[group_name]]$group), dataset_path = "/de")
    genes <- gene_groups[[group_name]]$genes

    dynamics <- pathways.dynamics(data=data, genes=genes,
                                  pseudobulks=pseudobulks, scale=T, signature=F,
                                  cell.types=cell.types)

    dynamics_not_scaled <- pathways.dynamics(data=data, genes=genes,
                                  pseudobulks=pseudobulks, scale=F, signature=F,
                                  cell.types=cell.types)

    features_names <- dynamics$fitted.vals$feature %>% unique

    labels <- sapply(features_names, function(name) {
      ct <- str_split(name, "\\.")[[1]][[1]]
      gene_name <- str_split(name, "\\.")[[1]][[2]]

      if (nrow(des %>% filter(gene == gene_name & cell.type == ct)) > 0) {
        significant <- '*'
      } else {
        significant <- ''
      }

      ct_name <- case_when(ct == 'inhibitory' ~ "Inh",
                           ct == 'cux2p' ~ "Exc L2-3",
                           ct == 'cux2m' ~ "Exc L4-6")


      label <- glue::glue("{gene_name} {ct_name} {significant}")

      label
    })

    withr::with_pdf(file.path(figures_dir, glue::glue("5Sup_gene_dynamics_{group_name}.pdf")), width=embed.width, height=embed.height, {
      print(plot.dynamics.wrapper(dynamics,
                                  features = features_names,
                                  labels=labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=1, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=T, include.points=FALSE, trajectories='prAD',
                                  show.label.legend=T, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.position="none",
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="scaled expresssion", title=glue::glue("{group_name} dynamics")))


      plot.genes.dynamics.heatmap(dynamics = dynamics_not_scaled, genes=features_names,
                                  row_labels=labels,
                                  cluster_rows = FALSE,
                                  row_order = order(labels),
                                  title=glue::glue("{group_name} dynamics")) %>% draw

    })
  }
})()

# Figure 4 -----------------

## Figure 4a - cell cycle and stress/apoptosis programs UMAP --------------

(function () {
  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct_number in 1:length(cts)) {
    ct <- cts[[ct_number]]
    embedding <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/umap.densities")) %>% tibble::column_to_rownames("cell")
    embedding <- embedding %>%
      dplyr::rename(x=UMAP_1, y=UMAP_2) %>%
      rename_at(vars(ends_with("sim.umap")), ~ gsub("X\\d+\\.(k\\d+).sim.umap",  "\\1", .)) # keep only the k and topic number


    apoptosis_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_early)[[ct_number]])
    cell_cycle_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_late)[[ct_number]])
    
    labels <- c("cell-cycle reentry program", "stress/apoptosis program")
    names(labels) <- c(cell_cycle_topic, apoptosis_topic)
    
    withr::with_pdf(file.path(figures_dir, glue::glue("4b_{ct}_umap_cell_cycle_vs_apoptosis.pdf")), height = embed.height.small*1.2, width = embed.width, {
      print(plot_umap_density2(embedding, columns = c(cell_cycle_topic, apoptosis_topic), ncol=2, labels = labels, title=ct))
    })
  }
})()




#+++++++++++++++++++++++++++
## Figure 4c, Exetended Data Figure 8a,b - programs correlation --------------
#+++++++++++++++++++++++++++

(function() {
  cell.types <- c("inhibitory", "cux2p", "cux2m")
  cell.types.names <- c("Inh", "Exc upper layers", "Exc lower layers")

  for (ct_index in 1:length(cell.types)) {
    ct <- cell.types[[ct_index]]
    names_to_topics <- lapply(neuronal_topics_to_annotation, function(x) str_extract(names(x)[ct_index], "k.*")) %>% unlist
    L <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L")) %>% column_to_rownames("cell")

    L <- L %>% dplyr::select(names_to_topics)

    c <- correlate_dataframes(L, L)

    add_significance <- TRUE
    if (add_significance) {
      corr.sig <- c$sig
      corr.sig[corr.sig == '****'] <- '**'
      corr.sig[corr.sig == '***'] <- '**'
      cell_fun <- function(j, i, x, y, w, h, fill) {
        grid.text(sprintf("%.2f\n%s", c$corr[i, j], corr.sig[i,j]), x, y, gp = gpar(fontsize = 8))
      }

    } else {
      cell_fun <- NULL
    }


    legend <- Legend(title = "adj.pval", pch = c("*", "**"), type = "points", labels = c("<0.05", "<0.01"))
    correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation)

    withr::with_pdf(file.path(figures_dir, glue::glue("4c_{ct}.pdf")), height = embed.height, width = embed.width, {
      ComplexHeatmap::Heatmap(c$corr,
                              column_title_gp = gpar(fontsize=10, fontface="bold"),
                              col = cyan2orange.correlation,
                              clustering_distance_rows = "spearman",
                              clustering_distance_columns = "spearman",
                              show_row_dend = FALSE,
                              show_column_dend = FALSE,
                              show_heatmap_legend=FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              column_names_rot = 90,
                              show_row_names = TRUE,
                              width = ncol(c$corr)*unit(5*5, "mm"),
                              height = ncol(c$corr)*unit(5*5, "mm"),
                              show_column_names = TRUE,
                              cell_fun = cell_fun
      ) %>% draw(annotation_legend_list=list(legend, correlation.legend),
                 column_title=glue::glue("{cell.types.names[[ct_index]]} - programs correlation across cells"))
    })
  }


})()


## Extended Data Figure 8c,d,e - program density and thresholds -------------
(function (){
  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct in cts) {
    quantiles_apoptosis <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/apoptosis_quantiles"))
    quantiles_cell_cycle <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_quantiles"))
    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))

    withr::with_pdf(file.path(figures_dir, glue::glue("8SuplD_{ct}_program_density.pdf")), height = embed.height*0.75, width = embed.width*1.25, {
      print(plot_topic_quantiles(df, "apoptosis_topic_score", quantiles_apoptosis, title = glue::glue('{ct}\napoptosis program histogram')))
      print(plot_topic_quantiles(df, "cell_cycle_topic_score", quantiles_cell_cycle, title = glue::glue('{ct}\ncell cycle program histogram')))
    })
  }
})()


#+++++++++++++++++++++++++++
## Figure 4g - validation increasing program scores correlation across datasets 
#+++++++++++++++++++++++++++
(function (){
  results <- list()

  for (project in VALIDATION_PROJECTS) {
    for (cell.type in c("inhibitory", "excitatory")) {
      # Read correlation matrix from HDF5
      corr_path <- paste0(project, "/signature_unique_correlation")
      
      corr_obj <- tryCatch({
        h5_read_correlation(cell.type = cell.type, group_name = corr_path)
      }, error = function(e) {
        message(glue::glue("Could not read signature_unique_correlation for {project} {cell.type}: {e$message}"))
        NULL
      })
      
      # Extract correlation matrix and program names
      cor_mat <- corr_obj$cor
      prog_names <- corr_obj$names[[1]]  # Program names from correlation matrix
      
      # Clean program names (remove SignatureUnique_ prefix if present, handle bytes)
      prog_names <- as.character(prog_names)
      prog_names <- gsub("^SignatureUnique_", "", prog_names)
      prog_names <- gsub("^k10$", "increasing-1", prog_names)
      prog_names <- gsub("^k8$", "increasing-2", prog_names)
      prog_names <- gsub("^increasing_early$", "increasing-1", prog_names)
      prog_names <- gsub("^increasing_late$", "increasing-2", prog_names)
      
      dimnames(cor_mat) <- list(prog_names, prog_names)
      dimnames(corr_obj$adj.pval) <- list(prog_names, prog_names)

      results[[paste(project, cell.type, sep = "_")]] <- data.frame(
        project = project,
        cell_type = cell.type,
        correlation = cor_mat["increasing-1", "increasing-2"],
        p_adj = corr_obj$adj.pval["increasing-1", "increasing-2"],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) {
    message("No correlations computed for increasing programs.")
    return()
  }

  plot_df <- bind_rows(results)
  
  # Add asterisks based on p_adj significance
  # * p < 0.05
  plot_df$asterisk <- ""
  plot_df$asterisk[!is.na(plot_df$p_adj) & plot_df$p_adj < 0.05] <- "*"
  
  # Create label with correlation and asterisk
  plot_df$label <- sprintf("%.2f%s", plot_df$correlation, plot_df$asterisk)

  plot_file_path <- file.path(figures_dir, "4g_increasing_program_correlation_barplot.pdf")
  withr::with_pdf(plot_file_path, width = 8, height = 5, {
    p <- ggplot(plot_df, aes(x = project, y = correlation)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6, fill = "#008585") +
      geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
      geom_text(aes(label = label),
                position = position_dodge(width = 0.7), 
                vjust = ifelse(plot_df$correlation >= 0, -0.3, 1.3), 
                size = 3.5) +
      facet_wrap(~cell_type, scales = "free_x") +
      labs(
        title = "Correlation between increasing programs (per validation dataset)",
        x = "Validation dataset",
        y = "spearman correlation"
      ) +
      theme_classic(base_size = 11) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    print(p)
  })

  message(glue::glue("Increasing program correlation barplot saved to {plot_file_path}"))
})()


## Figure 4f - spatial increasing pathways correlation 

(function () {
  anndata <- import("anndata", convert = FALSE)
  spatial_base_dir <- file.path(get_plots_dir(), "spatial")
  spatial_anndata_path <- file.path(spatial_base_dir, "for_paper/final_obj/spatial_obj_for_paper.h5ad")
  spatial_anndata <- anndata$read_h5ad(spatial_anndata_path,
                             backed = "r")
  
  obs <- spatial_anndata$obs %>% py_to_r
  obsm <- spatial_anndata$obsm %>% py_to_r

  df <- obs %>% filter(annotation=="Ex") %>% dplyr::select(c(
    'Ex_Metallothioneins_bind_metals_score', 'Ex_cell_cycle_score', 'Ex_DNA_Repair_score', 'Ex_MT_genes_score','Ex_response_to_oxidative_stress_score'
  ))
  
  rm(anndata, obs, obsm)
  
  c <- correlate_dataframes(df, method = "spearman")

  add_significance <- FALSE
  if (add_significance) {
    corr.sig <- c$sig
    corr.sig[corr.sig == '****'] <- '**'
    corr.sig[corr.sig == '***'] <- '**'
    cell_fun <- function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%.2f\n%s", c$corr[i, j], corr.sig[i,j]), x, y, gp = gpar(fontsize = 8))
    }

  } else {
    cell_fun <- function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%.2f\n", c$corr[i, j]), x, y, gp = gpar(fontsize = 8))
    }
  }
  
  col <- circlize::colorRamp2(
    c(-1, 0, 0.15, 0.3, 1),
    c("#4db6ac", "white", "#ffd79c", "#ef8c2b", "#ef6c00")
  )
  # legend <- Legend(title = "adj.pval", pch = c("*", "**"), type = "points", labels = c("<0.05", "<0.01"))
  correlation.legend <- Legend(title = "correlation", col_fun = col)

  withr::with_pdf(file.path(figures_dir, glue::glue("4f_spatial_increasing_pathways_correlation.pdf")), height = embed.height *1.5, width = embed.width*1.5, {
    ComplexHeatmap::Heatmap(c$corr,
                            column_title_gp = gpar(fontsize=10, fontface="bold"),
                            col = col,
                            clustering_distance_rows = "spearman",
                            clustering_distance_columns = "spearman",
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            show_heatmap_legend=FALSE,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            column_names_rot = 90,
                            show_row_names = TRUE,
                            width = ncol(c$corr)*unit(5*5, "mm"),
                            height = ncol(c$corr)*unit(5*5, "mm"),
                            show_column_names = TRUE,
                            cell_fun = cell_fun
    ) %>% draw(annotation_legend_list=correlation.legend,
                column_title=glue::glue("Spatial increasing pathways correlation"))
  })

  breaks <- seq(-1, 1, length.out = 101)
  pheatmap::pheatmap(
    df,
    breaks = breaks,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )
  
  rm(anndata)
})

## Figure 4d percent cell total -----

(function () {
  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct in cts) {
    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))

    withr::with_pdf(file.path(figures_dir, glue::glue("4d_{ct}_precent_expressing.pdf")), {
      print(df %>% mutate(apoptosis_cell=apoptosis_95 == 1,
                          cell_cycle_cell=cell_cycle_95 == 1,
                          both=apoptosis_cell & cell_cycle_cell) %>% ggvenn::ggvenn(
                            data = .,
                            columns = c("apoptosis_cell", "cell_cycle_cell"),
                            auto_scale=TRUE,
                            fill_color = c("#ef6c00", "#4db6ac"), show_outside='always'
                          ) + ggtitle(glue::glue("Percent expressing cells - {ct}"), glue::glue("quantiles - {quantile_threshold}")))
    })
  }
})()

# Figure 4h, Extended Data Figure 8j,k -  numbers of cell cycle cells start vs prAD

(function (){
  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct in cts) {
    prop_cells_df <- get_AD_programs_proportion(ct)
    prop_cells_df <- prop_cells_df %>% filter(!is.na(group) & group != 'NA')
    quantiles <- c("95", "90", "80")
    
    withr::with_pdf(file.path(figures_dir, glue::glue("4h_number_of_cell_start_vs_prAD_{ct}.pdf")), width = embed.width, height = embed.height, {
      for (q in quantiles) {
        print(ggbarplot(prop_cells_df, x="group",
                        y=glue::glue("prop_cell_cycle_{q}"), fill="group",  width = 0.9,
                        add = c("mean_se")) + 
                stat_compare_means(comparisons=list(c('prAD', 'start')), label = "p.signif") +
                theme_minimal() +
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
                labs(title=glue::glue("cell cycle cells propotion - {q}%"), x='group', y='% of neurons') +
                scale_fill_manual(values=c("start"="#74a892",
                                           "prAD"="#c7522a")))
      }

    })
  }
})()


## Figure 4i - DNA copies data --------------

(function () {
  # taken from https://doi.org/10.1523/JNEUROSCI.0379-07.2007
  data_string <- "
Group,2n_Mean,2n_SEM,2n-4n_Mean,2n-4n_SEM,4n_Mean,4n_SEM,n_samples
Control I,298.3,13.4,36.8,13.2,1.2,1.0,6
Control II,250.2,14.3,34.1,14.1,1.1,0.8,7
AD Early,186.8,4.8,50.9,4.8,5.8,1.7,6
AD Advanced,123.1,12.2,51.0,14.1,2.1,1.2,7"

  data_df <- read.csv(textConnection(data_string), header = TRUE, check.names = FALSE)

  data_df <- data_df %>%
    pivot_longer(
      cols = starts_with("2n") | starts_with("4n"),  # Columns to pivot
      names_to = c("n", ".value"),  # Split into two columns
      names_sep = "_"  # Separator between "2n" and "Mean/SEM"
    )

  data_df <- data_df %>% mutate(Group=factor(Group, levels=c("Control I", "Control II", "AD Early", "AD Advanced")))

  withr::with_pdf(file.path(figures_dir, "4i_hyperploidity.pdf"), width = embed.width.small * 1.5, height = embed.height.small, {
 `   print(data_df %>% group_by(Group) %>% mutate(prop = Mean /sum(Mean)) %>%
      ungroup %>%
      ggplot(aes(x=n, y=prop, fill=Group)) + geom_bar(stat = 'identity', position = 'dodge') +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(plot.title = element_text(hjust = 0.5, size=12), legend.key.size = unit(0.2, "cm")) +
      labs(title="hyperploidity propotion", x='hyperploidity', y='% of neurons') +
      scale_fill_manual(values=c("Control I"="#74a892",
                                 "Control II"="#008585",
                                 "AD Early"="#c7522a",
                                 "AD Advanced"="#f0daa5")))`
  })

})


## Figure 5a - vulnerability vs cell cycle and apoptosis -------
(function (){
  SUBTYPE_MIN_N_CELLS <- 1000
  
  neuronal_classes <- c("inhibitory", "excitatory") 
  res <- sapply(neuronal_classes, simplify = FALSE, function(neuronal_class) {
    if (neuronal_class == 'inhibitory') {
      df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["inhibitory"]], "/cell_cycle_apoptosis"))
    } else {
      df_cux2p <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2p"]], "/cell_cycle_apoptosis"))
      df_cux2m <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["cux2m"]], "/cell_cycle_apoptosis"))
      
      df <- rbind(df_cux2p, df_cux2m)
      rm(df_cux2m, df_cux2p)
    }
    
    removed_subtypes <- df %>%
      dplyr::count(pruned_subtype) %>%
      filter(n < SUBTYPE_MIN_N_CELLS) %>%
      pull(pruned_subtype)
    
    print(glue::glue("Removing subtypes: {paste(removed_subtypes)}"))


    summary_data <- lapply(c(`cell cycle`="cell_cycle_95", apoptosis="apoptosis_95"), function(column) {
      percent_cells_df <- df %>%
        filter(pruned_subtype != 'NA' & !(pruned_subtype %in% removed_subtypes)) %>%
        group_by(projid, group, pruned_subtype) %>% summarise(above=mean(!!sym(column)))
      
      summary_data <- percent_cells_df %>%
        filter(!is.na(group) & group != 'NA') %>%
        mutate(group=factor(group, levels=c('start', 'prAD'))) %>%
        group_by(group, pruned_subtype) %>%
        summarise(
          mean = mean(above),
          sem = sd(above) / sqrt(n())
        )
    })
    
    subtype_order <- summary_data$`cell cycle` %>%
      filter(group == 'prAD') %>%
      arrange(mean) %>% pull(pruned_subtype)
    
    traits <- associate.allen.traits(df = df, data=data, exclude_list = removed_subtypes)
    vuln_mat <- cbind(apoptosis=summary_data$apoptosis %>% pivot_wider(names_from = 'group', id_cols = 'pruned_subtype', values_from = 'mean') %>% column_to_rownames('pruned_subtype'),
                      cell_cycle=summary_data$`cell cycle` %>% pivot_wider(names_from = 'group', id_cols = 'pruned_subtype', values_from = 'mean') %>% column_to_rownames('pruned_subtype'))
    
    vuln_mat <- vuln_mat[subtype_order,
                         c('cell_cycle.start', 'apoptosis.start', 'cell_cycle.prAD', 'apoptosis.prAD'),
                         drop=FALSE]
    
    list(vuln_mat=vuln_mat, traits=traits)
  })
  
  max_percent_cells <- max(res$inhibitory$vuln_mat, res$excitatory$vuln_mat)
  max_tstat <- max(abs(res$inhibitory$traits$tstat), abs(res$excitatory$traits$tstat))
  
  for (neuronal_class in neuronal_classes) {
    vuln_mat <- res[[neuronal_class]]$vuln_mat
    
    vuln_mat <- vuln_mat[!grepl(".start", colnames(vuln_mat))]
    traits <- res[[neuronal_class]]$traits
    
    withr::with_pdf(file.path(figures_dir, glue::glue("5a_{neuronal_class}_subtypes_cell_cycle_apoptosis.pdf")), {
      annotation <- HeatmapAnnotation(`neuronal state`=sapply(str_split(colnames(vuln_mat), "\\."), `[`, 1),
                                      annotation_name_side='left',
                                      col=list(`neuronal state`=c(apoptosis="#ef6c00", `cell_cycle`="#4db6ac")))
      
      tratit_hm_and_legend <- plot.trait.associations.heatmap(traits,
                                                              max_color_value = max_tstat,
                                                              row_order=rownames(vuln_mat),
                                                              cluster_rows=FALSE,
                                                              column_title ='T.A',
                                                              params = c('sqrt.amyloid_mf', 'sqrt.tangles_mf'),
                                                              width=3*unit(4, "mm"),
                                                              row.by = 'covariate', show.only.significant = FALSE,
                                                              plot=FALSE)
      
      draw(Heatmap(vuln_mat %>% as.matrix(), 
                   cluster_rows = FALSE, name="Cells propotion",
                   cluster_columns = FALSE,
                   width=4*unit(7, "mm"), column_title='Cells propotion',
                   height = length(subtype_order)*unit(7, "mm"),
                   col = circlize::colorRamp2(c(0, max_percent_cells), c("white", "red")),
                   bottom_annotation = annotation) +tratit_hm_and_legend$hm
      )
    }) 
  }
})()


## Figure 5b - programs dynamics in SST and PV ------
(function (){
  ct <- 'inhibitory'
  df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))
  dynamics_features <- df %>%
    filter(pruned_subtype != 'NA' ) %>%
    group_by(projid, pruned_subtype) %>% summarise(cell_cycle_topic_score=mean(cell_cycle_topic_score),
                                                   cell_cycle_topic_n=mean(`cell_cycle_95`),
                                                   apoptosis_topic_score=mean(apoptosis_topic_score),
                                                   apoptosis_topic_n=mean(`apoptosis_95`),
                                                   decreasing_late_topic_score=mean(decreasing_late_topic_score),
                                                   decreasing_early_topic_score=mean(decreasing_early_topic_score))
  
  withr::with_pdf(file.path(figures_dir, "5b_sst_pv_topic_dynamics.pdf"), {
    for (c in names(tidyselect::eval_select(quote(c(-projid, -pruned_subtype)), dynamics_features))) {
      features <- dynamics_features %>%
        dplyr::select(projid, pruned_subtype, all_of(c)) %>%
        mutate(pruned_subtype=str_to_upper(pruned_subtype)) %>%
        tidyr::pivot_wider(names_from=pruned_subtype, values_from = all_of(c), values_fill=NA) %>%
        filter(!is.na(projid)) %>%
        column_to_rownames("projid") %>%
        filter(rownames(.) %in% rownames(data)) %>%
        arrange(match(rownames(.), rownames(data))) %>%
        as.matrix()
      
      if (str_detect(c, "score")) {
        features <- features %>% sqrt
      }
      
      dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                               features = features,
                               trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                               trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                               evaluate.fit = T,
                               bootstrap = F)
      
      print(plot.dynamics.wrapper(dynamics,
                                  features =  c("SST", "PVALB"),
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  cols=c("SST"="#ef6c00", "PVALB"="#4db6ac"),
                                  ncol=2, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=T, include.points=FALSE, trajectories='prAD',
                                  show.label.legend=F, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.text = element_text(margin = margin(0,0,0,0)),
                    legend.position="none", legend.box="vertical", legend.margin=margin(0, 0, 0, 0),
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y=glue::glue("{c}"), title=glue::glue("{ct} {c} dynamics")))
    }
  })
})()

# Figure 4j -  senescence like genes signature -----
(function (){
  cts <- c("inhibitory", "cux2m"
    #, "cux2p"
  )
  for (ct in cts) {
    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))
    senescence_like_genes_df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/senescence_like_genes.signatures"))
    stopifnot(all(df$cell == senescence_like_genes_df$cell))

    df <- df %>% cbind(senescence_like_genes_df %>% dplyr::select(-cell, -projid))

    df <- df %>% mutate(cell_fate=dplyr::case_when(`cell_cycle_95` & `apoptosis_95` ~ 'both',
                                                   `cell_cycle_95` & !`apoptosis_95` ~ 'cell-cycle reentry',
                                                   !`cell_cycle_95` & `apoptosis_95` ~ 'stress/apoptosis',
                                                   TRUE ~ "neither")) %>%
      filter(cell_fate!='both') %>%
      mutate(cell_group = case_when(group == 'prAD' ~ paste(group, cell_fate, sep=" "),
                                    group == 'start' &  cell_fate == 'neither' ~ 'start neither',
                                    TRUE ~ 'other')) %>%
      filter(cell_group != 'other') %>% 
      mutate(cell_group=factor(cell_group, levels=c("start neither", "prAD neither", "prAD cell-cycle reentry", "prAD stress/apoptosis")))
    

    senescence_like_signatures <- list(
      `Checkpoints`=c("ATR", "RAD9A", "BRCA1"),
      `DNA damage repair`=c("ATM"),
      `Translision`=c("REV1", "POLK"),
      `Cell cycle`=c("CENPC", "ANAPC16", "CEP290", "STAG2"),
      `PI3K/AKT pathway`=c("PIK3CA", "AKT3"),
       NFKB1=c("NFKB1", "TRAF6"),
       `anti apoptotic`=c("BANP"),
       PUMA=c("BBC3")
    )
    column_group <- unlist(mapply(function(group, genes) rep(group, length(genes)), names(senescence_like_signatures), senescence_like_signatures), use.names = FALSE) 
    senescence_like_genes <- unlist(senescence_like_signatures, use.names = FALSE)
    withr::with_pdf(file.path(figures_dir, glue::glue("4j_SenescenceCellCycleDotplot_by_group_{ct}.pdf")), width = embed.width*1.5, height = embed.height, {
      df <- df %>% filter(cell_group != 'prAD neither')
      print(dot_plot(df , senescence_like_genes, group.vector = df[["cell_group"]], alpha.threshold = 0, scale = TRUE, mean_threshold = -0.1, 
        column.group=column_group, order=FALSE) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title=glue::glue("Senescence related genes in {ct}"), x="", y=""))
    })
  }
})()

## Figure 5d - synaptic signatures in cell cycle or apoptosis -----
(function (){
  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct in cts) {
    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))
    
    df <- df %>% mutate(cell_fate=dplyr::case_when(`cell_cycle_95` & `apoptosis_95` ~ 'both',
                                                   `cell_cycle_95` & !`apoptosis_95` ~ 'cell-cycle reentry',
                                                   !`cell_cycle_95` & `apoptosis_95` ~ 'stress/apoptosis',
                                                   TRUE ~ "neither")) %>%
      filter(cell_fate!='both') %>%
      mutate(cell_group = case_when(group == 'prAD' ~ paste(group, cell_fate, sep=" "),
                                    group == 'start' &  cell_fate == 'neither' ~ 'start neither',
                                    TRUE ~ 'other')) %>%
      filter(cell_group != 'other') %>% 
      mutate(cell_group=factor(cell_group, levels=c("start neither", "prAD neither", "prAD cell-cycle reentry", "prAD stress/apoptosis")))
    
    cell_group_counts <- df %>%
      group_by(cell_group) %>%
      summarise(n=n())
    
    features <- c("Synaptic vescile signature"="presynaptic_genes",
                  "Postsynaptic component signature"="postsynaptic_genes",
                  "Decreasing early program"="decreasing_early_topic_score",
                  "Decreasing late program"="decreasing_late_topic_score")
    
    for (feature_name in names(features)) {
      feature <- features[[feature_name]]  
      
      legend_df <- data.frame(signif = c("*", "**", "***"), x = 1, y = 1)
      
      comparisons_list <- list(c("start neither", "prAD cell-cycle reentry"),
                               c("start neither", "prAD stress/apoptosis"),
                               c("prAD cell-cycle reentry", "prAD stress/apoptosis"))
      p_values <- sapply(comparisons_list, function(comp) {
        test_result <- wilcox.test(formula=as.formula(glue::glue("{feature} ~ cell_group")),
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
      
      withr::with_pdf(file.path(figures_dir, glue::glue("5d_siganture_{ct}_{feature}.pdf")), width = embed.width, height = embed.height, {
        print(ggviolin(df, x="cell_group", fill="cell_group", y=feature, legend.title='cell group') + 
                stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") + 
                geom_signif(
                  comparisons = comparisons_list,
                  annotations = annotation_df$label,
                  y_position = c(1, 1.2, 1.4),
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
                  aes(x = cell_group, y = max(df[[feature]]) + 1.5, label = paste0("n=", n)),
                  inherit.aes = FALSE
                ) +
                scale_fill_manual(values=c('start neither'='#9dc56c',
                                           'prAD neither'='red',
                                           'prAD cell-cycle reentry'='#4db6ac',
                                           'prAD stress/apoptosis'="#ef6c00")) + 
                scale_x_discrete(guide = guide_axis(angle = 45)) +
                theme(legend.position = "bottom", legend.box = "vertical") +
                labs(y=feature_name) +
                ggtitle(glue::glue("Signature score in prAD vs start neither cells  - {ct}")))
      })
      
    }
  }
})()


## Extended Data Figure 9c - neurotransmitters expression violins
(function (){
  ct <- 'inhibitory'
  neurotransmitters <- c("SST", "VIP")
  df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))

  df <- df %>% mutate(cell_fate=dplyr::case_when(`cell_cycle_95` & `apoptosis_95` ~ 'both',
                                                `cell_cycle_95` & !`apoptosis_95` ~ 'cell-cycle reentry',
                                                !`cell_cycle_95` & `apoptosis_95` ~ 'stress/apoptosis',
                                                TRUE ~ "neither")) %>%
              filter(cell_fate!='both') %>%
              mutate(cell_group = case_when(group == 'prAD' ~ paste(group, cell_fate, sep=" "),
                                            group == 'start' &  cell_fate == 'neither' ~ 'start neither',
                                            TRUE ~ 'other')) %>%
              filter(cell_group != 'other') %>%
              mutate(cell_group=factor(cell_group, levels=c("start neither", "prAD neither", "prAD cell-cycle reentry", "prAD stress/apoptosis")))

  for (nt in neurotransmitters) {
    subtype_df <- df %>% filter((pruned_subtype == str_to_title(nt)) | (nt == 'NPY' & pruned_subclass == 'Inh L1-6 SST NPY'))
    subtype_counts <- subtype_df %>%
        group_by(cell_group) %>%
        summarise(n=n())

    legend_df <- data.frame(signif = c("*", "**", "***"), x = 1, y = 1)

    comparisons_list <- list(c("start neither", "prAD cell-cycle reentry"),
                             c("start neither", "prAD stress/apoptosis"),
                             c("prAD cell-cycle reentry", "prAD stress/apoptosis"))
    p_values <- sapply(comparisons_list, function(comp) {
      test_result <- wilcox.test(formula=as.formula(glue::glue("{nt} ~ cell_group")),
                                 data = subtype_df %>% filter(cell_group %in% comp))
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

    MAX_EXPRESSION <- 4
    stopifnot(nrow(subtype_df %>% filter(.data[[nt]] > MAX_EXPRESSION)) < 30)

    withr::with_pdf(file.path(figures_dir, glue::glue("9SupC_neurotransmitters_levels_in_group_{nt}.pdf")), {
      filtered_subtype_df <- subtype_df %>% filter(.data[[nt]] < MAX_EXPRESSION)
      print(ggviolin(filtered_subtype_df, x="cell_group", fill="cell_group", y=nt, legend.title='cell group') +
        stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") +
        geom_signif(
          comparisons = comparisons_list,
          annotations = annotation_df$label,
          y_position = c(4.5, 5, 5.5),
          tip_length = 0.02
        ) +
        geom_point(data = legend_df, aes(x = x, y = y, shape = signif), color = "black", size = 4) +
        scale_shape_manual(
          name = "Significance",
          values = c("Significance" = 8),  # Use star shape for legend
          labels = c(glue::glue("* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001"))
        ) +
        geom_text(
          data = subtype_counts,
          aes(x = cell_group, y = max(filtered_subtype_df[[nt]]) + 1, label = paste0("n=", n)),
          inherit.aes = FALSE
        ) +

        scale_fill_manual(values=c('prAD cell-cycle reentry'='#4db6ac',
                                   'prAD stress/apoptosis'="#ef6c00",
                                   'start neither'='#9dc56c',
                                   'prAD neither'='red')) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(legend.position = "bottom", legend.box = "vertical") +
        labs(y=glue::glue("{nt} expression")) +
        ggtitle(glue::glue("{nt} expression in prAD vs start neither cells")))
    })
  }
})()



## Extended Data Figure 7d - pre synaptic proteins --------------

(function () {
  pathway_ids <- c("GO:0098693", #  regulation of synaptic vesicle cycle
                   "GO:0099504", # Synaptic vesicle cycle
                   "hsa04721" # Synaptic vesicle cycle
                   )


  de <- read.dataset.for.topics(names(neuronal_topics_to_annotation$increasing_early), dataset_path = "/de")
  pathway.genes <- unlist(lapply(pathway_ids, get.pathway.genes)) %>% unique
  de <- de %>% filter(gene %in% pathway.genes)
  topic.to.genes <- sapply(de %>% pull(cell.type) %>% unique, function(ct) {
    de %>% filter(cell.type == ct) %>% pull(gene)
  }, simplify = FALSE)

  all_genes <- de$gene %>% unique


  traits_columns <- c( "sqrt_amyloid_mf", "sqrt_tangles_mf", "cogng_random_slope")
  colData(proteom)$sqrt_amyloid_mf <- sqrt(colData(proteom)$amyloid_mf)
  colData(proteom)$sqrt_tangles_mf <- sqrt(colData(proteom)$tangles_mf)

  prot_traits <- associate.proteomic.traits(proteom, genes=all_genes, traits = traits_columns)


  withr::with_pdf(file.path(figures_dir, "7SuplD_presynaptic_proteins.pdf"), height = embed.height * 1.5, width = embed.width, {
    plot.trait.assocation.multiple.groups(prot_traits,
                                          row.by = 'UniProt', params = traits_columns,
                                          show.only.significant=FALSE,
                                          row_gap=unit(5, "mm"),
                                          cluster_row_slices=FALSE,
                                          column_labels=AD_TRAITS_NAMES,
                                          row_labels=prot_traits$gene,
                                          column_title="presynaptic proteins")

  })
})()


# Figure 5 Program Communities -----------------

## Extended Data Figure 9d - correlation and dynamics similarity  --------------

(function (){
  library(dendextend)
  communities   <- data$uns$communities

  var_names <- data$var %>% filter(!is.na(community)) %>% rownames()
  membership    <- split(data$var_names, data$var$community)
  dynamics      <- communities$similarities$dynamics %>% `dimnames<-`(list(communities$names, communities$dynamics.colnames))
  corr          <- communities$similarities$correlation %>% `dimnames<-`(list(communities$names, communities$names))
  dyn.adjacency <- communities$similarities$dynamics.adjacency %>% `dimnames<-`(list(communities$names, communities$names))


  # States to annotate
  mark.states <- c("Ast.10.k2", "Ast.10.k6", "Ast.10.k5",
                   "Micro.15.k11", "Micro.15.k5", "Oligo.8.k8",
                   "Inh.17.k10", "Inh.17.k1", "Inh.17.k7", "Inh.17.k8")

  # Specify dynamics column-order: zero pseudotime in the middle
  dynamics <- dynamics[,c(rev(colnames(dynamics)[grepl("ABA", colnames(dynamics))]),
                          colnames(dynamics)[grepl("prAD", colnames(dynamics))])]

  # Hierarchically split community's membership, creating dynamics heatmap within each community/sub-community
  hm.comms <- data$var %>% split(., as.character(.$community)) %>%
    lapply(., function(inner) {

      if(any(!is.na(inner$sub.community)))
        lst <- split(rownames(inner), as.character(inner$sub.community))
      else
        lst <- list(rownames(inner))

      lapply(names(lst), function(sub.community) {
        states <- lst[[sub.community]]
        mtx  <- dyn.adjacency[states,states, drop=FALSE] + corr[states,states, drop=FALSE]
        if (nrow(mtx) > 1) {
          dend <- dendsort::dendsort(hclust(dist(mtx))) %>% as.dendrogram() %>% dendextend::set("labels_to_character")
        } else {
          dummy.mtx <- rbind(cbind(mtx, mtx), cbind(mtx, mtx))
          dimnames(dummy.mtx) <- list(c(rownames(mtx), "dummy"), c(rownames(mtx), "dummy"))
          dend <- dendsort::dendsort(hclust(dist(dummy.mtx))) %>% as.dendrogram() %>% dendextend::set("labels_to_character") %>% prune(c("dummy"))
        }

        prepare(Heatmap(dynamics[states,, drop=FALSE],
                        cluster_rows = dend,
                        cluster_row_slices = F,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = T,
                        show_row_dend = T,
                        name=sub.community,
                        row_dend_width = unit(10,"pt"),
                        col = circlize::colorRamp2(seq(-4,4,length.out=21),
                                                   colorRampPalette(c("darkorchid4","white","#E65100"))(21)),
                        column_title = NULL,
                        row_names_side = "left",
                        column_gap = unit(1, "pt"),
                        ))

      }) %>% unlist(., recursive=FALSE)
    })  %>% unlist(., recursive=FALSE)

  # Retrieve program order to be used in composite adjacency matrix - Fig 5c
  program.order <- do.call(rbind, lapply(hm.comms, function(hm) data.frame(comm = hm@name, program=rownames(hm@matrix)[row_order(hm)])))

  write.csv(program.order, file.path(figures_dir, "community_programs_order.csv"))

  # Plot dynamics adjacency and state-state correlations in lower/upper triangles
  dyn.adjacency <- dyn.adjacency[program.order$program, program.order$program]
  corr <- corr[program.order$program, program.order$program]

  marks <- sapply(mark.states, function(s) which(rownames(dyn.adjacency) == s))
  matrices <- list(dyn.adjacency, corr)
  colors   <- list(circlize::colorRamp2(seq(0,1,length.out=21), colorRampPalette(c("white","salmon","red","firebrick4"))(21)),
                   circlize::colorRamp2(seq(-1,1,length.out=21), blue2orange(21)))

  
  program.order$comm <- factor(program.order$comm, levels=program.order$comm %>% unique)

  legends <- list(
    Legend(col_fun = colors[[1]], title = "Dynamics similarty"),
    Legend(col_fun = colors[[2]], title = "Spearman correlation")

  )

  withr::with_pdf(file.path(figures_dir, "9SupD.pdf"), width=embed.width*2, height=embed.height*2, {
    Heatmap(
      matrices[[i]],
      row_split = program.order$comm,
      column_split = program.order$comm,
      cluster_rows = F,
      cluster_columns = F,
      col = colors[[i]],
      rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, w, h, col) {
        if(i==j) grid.rect(x,y,w,h, gp = gpar(fill="black", col=NA))
        else if(i > j) grid.rect(x,y,w,h, gp = gpar(fill=colors[[1]](matrices[[1]][i, j]), col=NA))
        else if(i<j) grid.rect(x,y,w,h, gp = gpar(fill=colors[[2]](matrices[[2]][i, j]), col=NA))

      },
      column_title = NULL,
      show_heatmap_legend = F,
      show_row_names = F,
      show_column_names = F,
      row_gap = unit(5, "pt"),
      column_gap = unit(5, "pt"),
      right_annotation = rowAnnotation(states = anno_mark(marks, topic_names_mapping(names(marks), use_static_mapping = FALSE)))
    ) %>% draw(annotation_legend_list=legends)
  })
})()


## Extended Data Figure 9e, Extended Data Figure 10a - Community analysis Validation  --------------

(function(){
   convert_program_name <- function(name) {
    program_name_to_group <- inner_to_root(neuronal_topics_to_annotation)
    if (name %in% names(program_name_to_group) & grepl("^Exc\\.", name) ) {
      name <- paste0("excitatory_", program_name_to_group[[name]])
      return(name)
    }

    name <- gsub("Exc.", "", name)
    name <- gsub("Inh", "inhibitory", name)
    name <- gsub("Micro", "microglia", name)
    name <- gsub("Ast", "astrocytes", name)
    name <- gsub("Oligo", "oligo", name)
    name <- gsub("OPC", "opcs", name)
    name <- tolower(name)
    name <- gsub("^([^.]+)\\.\\d+\\.(k\\d+)$", "\\1_\\2", name)

    
    name
  }

  programs_community <- data$var %>% filter(!is.na(community))
  programs_community <- programs_community %>%
    mutate(program_name=rownames(.)) %>%
    rowwise %>% 
    mutate(short_program_name=convert_program_name(program_name)) %>% 
    ungroup %>%
    distinct(short_program_name, .keep_all = TRUE)
  
  programs_order <- read.csv(file.path(figures_dir, "community_programs_order.csv")) %>%
   mutate(community=gsub("\\.\\d+$", "", comm)) %>%
   rowwise %>%
   mutate(short_program_name=convert_program_name(program)) %>%
   ungroup %>% 
   distinct(short_program_name, .keep_all = TRUE)

  for (project in VALIDATION_PROJECTS) {
    print(project)
    c <- h5_read_correlation(group_name=paste0("/validation/", project, "/programs_correlation"))
    
    not_in_community_index <- c$names[[1]] %in% (data$var %>% filter(is.na(community)) %>% rownames() %>% sapply(convert_program_name, USE.NAMES = FALSE))
    
    c_subset <- filter_correlation_by_index(c, !not_in_community_index)

    NonRPLRPS_programs <- grepl("RPLRPS", c_subset$names[[1]])
    if (any(NonRPLRPS_programs)) {
      c_subset <- filter_correlation_by_index(c_subset, !NonRPLRPS_programs)
    }

    community_colors <- c("C1" = "#008585", "C2" = "#7c242c", "C3" = "#194a7a")
    subcommunity_colors <- c("C2.1"= "#c7522a", "C2.2"="#ef6c00", "C2.3"="#7c242c", "C1.1"="#74a892", "C1.2"="#008585", "C3.1" = "#194a7a", "C3.2"="#7593af")
    mark_programs <- c("astrocytes_k2"= "AST.p2", "astrocytes_k6"= "AST.p6", "astrocytes_k5"= "AST.p5",
                   "microglia_k11"= "MIC.p11", "microglia_k5"= "MIC.p5", "oligo_k8"= "OLI.p8",
                   "inhibitory_k10"= "INH.p10", "inhibitory_k1"= "INH.p1", "inhibitory_k7"= "INH.p7", "inhibitory_k8"= "INH.p8")

    mark_programs <- mark_programs[names(mark_programs) %in% c_subset$names[[1]]]
    current_programs <- programs_order %>% filter(short_program_name %in% c_subset$names[[1]])

    row_anno <- rowAnnotation(
      Community = programs_community$community[match(c_subset$names[[1]], programs_community$short_program_name)],
      Subcommunity = programs_community$sub.community[match(c_subset$names[[1]], programs_community$short_program_name)],
      program = anno_mark(at = which(c_subset$names[[1]] %in% names(mark_programs)), 
      labels = mark_programs[c_subset$names[[1]][c_subset$names[[1]] %in% names(mark_programs)]]),
      col = list(Community = community_colors, Subcommunity = subcommunity_colors)
    )
    
    dimnames(c_subset$cor) <- list(c_subset$names[[1]], c_subset$names[[2]])
    
    withr::with_pdf(file.path(figures_dir, glue::glue("9SupE_{project}_community_validation_with_sub.pdf")), width=embed.width * 2, height=embed.height * 2, {
      Heatmap(
        c_subset$cor, 
        col = circlize::colorRamp2(seq(-1,1,length.out=21), blue2orange(21)),
        row_order = match(current_programs$short_program_name, c_subset$names[[1]]), 
        column_order = match(current_programs$short_program_name, c_subset$names[[2]]),
        left_annotation = row_anno,
        width = unit(8, "cm"),
        height = unit(8, "cm"),
        name="correlation",
        cell_fun = NULL,
        show_column_names = FALSE,
        show_row_names = FALSE,
      ) %>% draw(column_title=glue::glue("{project} - Community correlation across individuals"))
    })
  }
})()



## Figure 6a, Extended Data Figure 6a  - Communities heatmap --------------

(function (){
  communities <- data$var %>% filter(!is.na(community)) %>% pull(community) %>% unique
  all.sub.communities <- data$var %>% filter(!is.na(sub.community)) %>% pull(sub.community) %>% unique %>% sort


  in_community_var <- data$var %>% filter(!is.na(community))


  withr::with_pdf(file.path(figures_dir, "6a_topics_heatmap_by_communities.pdf"), width=embed.width, height=embed.height * 1.5, {
    plot.genes.dynamics.heatmap(data$uns$trajectories$palantir$dynamics,
                                genes = rownames(in_community_var),
                                title="@{x}", row_names_mapping_func = pryr::partial(topic_names_mapping, use_static_mapping=FALSE),
                                both_trajectories = TRUE, show_lfc = F, name = 'Scaled program score',
                                row_split = unfactor(in_community_var$sub.community),
                                clustering_distance_rows='spearman',
                                use_slanter_ordering=TRUE,
                                in_community_var=in_community_var,
                                row_gap = unit(2, "mm")
                                ) %>% draw(column_title='Pseudotime',
                                           column_title_side='bottom',
                                           row_title='Sub-communities',
                                           row_title_side='left')
  })
})()


## Figure 6b - sub communities assoication --------------

(function (){
  all.sub.communities <- data$var %>% filter(!is.na(sub.community)) %>% pull(sub.community) %>% unique %>% sort
  withr::with_pdf(file.path(figures_dir, "6b_communities_trait_association.pdf"), width = embed.width.small, height = embed.height, {
    plot.trait.associations.heatmap(data$uns$communities$trait.association,
                                    row.by = 'covariate',
                                    height=10*unit(8, "mm"),
                                    width=unit(30, "mm"),
                                    column_labels=AD_TRAITS_NAMES,
                                    show.only.significant = F)
  })
})()


## Figure 6c,d and 7b,c - cross cell types pathways  --------------

(function() {
  pathways <- list(
    "Increasing" = list(
      c("GO:0090077", "foam cell differentiation", "microglia", "k11"),
      c("GO:0051402", "neuron apoptotic process", "inhibitory", "k10"),
      # c("GO:0051402", "neuron apoptotic process", "cux2p", "k11"),
      # c("GO:0051402", "neuron apoptotic process", "cux2m", "k3"),

      c("GO:0006979", "response to oxidative stress", "inhibitory", "k10"),
      # c("GO:0006979", "response to oxidative stress", "cux2p", "k11"),
      # c("GO:0006979", "response to oxidative stress", "cux2m", "k3"),

      c("GO:0032365", "intracellular lipid transport", "microglia", "k11"),
      c("GO:0050663", "cytokine secretion", "microglia", "k8"),
      c("GO:0006911", "phagocytosis, engulfment", "microglia", "k8")
    ),
    "Cellular response to heat stress" = list(
      #  c("R-HSA-3371556", "Cellular response to heat stress", "inhibitory", "k10"),
      #  c("R-HSA-3371556", "Cellular response to heat stress", "cux2p", "k11"),
      #  c("R-HSA-3371556", "Cellular response to heat stress", "cux2m", "k3"),

       c("R-HSA-3371556", "Cellular response to heat stress", "inhibitory", "k8"),
       c("R-HSA-3371556", "Cellular response to heat stress", "cux2p", "k4"),
       c("R-HSA-3371556", "Cellular response to heat stress", "cux2m", "k9"),

       c("R-HSA-3371556", "Cellular response to heat stress", "oligo", "k8"),
       c("R-HSA-3371556", "Cellular response to heat stress", "astrocytes", "k6"),
       c("R-HSA-3371556", "Cellular response to heat stress", "microglia", "k13")
     ),

     "Metallothioneins bind metals" = list(
       c("R-HSA-5661231", "Metallothioneins bind metals", "inhibitory", "k10"),
       c("R-HSA-5661231", "Metallothioneins bind metals", "cux2p", "k11"),
       c("R-HSA-5661231", "Metallothioneins bind metals", "cux2m", "k3"),

       # not in communities
       c("R-HSA-5661231", "Metallothioneins bind metals", "astrocytes", "k7"),
       c("R-HSA-5661231", "Metallothioneins bind metals", "microglia", "k6"),
       c("R-HSA-5661231", "Metallothioneins bind metals", "oligo", "k2"),
       c("R-HSA-5661231", "Metallothioneins bind metals", "opcs", "k3")
     ),

     # synaptic
     "Synaptic pathways" = list(
       c("GO:0006836", "neurotransmitter transport", "astrocytes", "k5"),
       c("GO:0098883", "synapse pruning", "microglia", "k12"),

       c("hsa04721", "Synaptic vesicle cycle", "inhibitory", "k10"),
       # c("hsa04721", "Synaptic vesicle cycle", "cux2p", "k11"),
       # c("hsa04721", "Synaptic vesicle cycle", "cux2m", "k3"),
       # c("GO:0099572", "postsynaptic specialization", "inhibitory", "k7"),
       # c("GO:0099572", "postsynaptic specialization", "cux2p", "k10"),
       c("GO:0099572", "postsynaptic specialization", "cux2m", "k8")
       # c("GO:0007215", "glutamate receptor signaling pathway", "opcs", "k7")

     ),

    "OligoMylin" = list(
      c("GO:0043209", "myelin sheath", "oligo", "k6"),
      c("GO:0031641", "regulation of myelination", "oligo", "k6"),
      c("R-HSA-191273", "Cholesterol biosynthesis", "oligo", "k8"),
      c("GO:0010878", "cholesterol storage", "microglia", "k11")
    )
  )

  for (group in names(pathways)) {
    message(group)

    current_pathways <- pathways[[group]]
    topics <- lapply(current_pathways, function(c) paste0(c[[3]], ".", c[[4]])) %>% unlist %>% unique
    des <- read.dataset.for.topics(topics, dataset_path = "/de")

    pathways_scores <- sapply(1:length(current_pathways), function(i) {
      pathway_id <- current_pathways[[i]][[1]]
      cell.type <- current_pathways[[i]][[3]]
      k <- current_pathways[[i]][[4]]

      pathway.genes <- get.pathway.genes(pathway_id)
      de <- des %>% filter(gene %in% pathway.genes & .env[['cell.type']] == .data[['cell.type']] & topic == k)

      if (cell.type == 'cux2p') {
        bulk.cell.type <- 'cux2+'
      } else if (cell.type == 'cux2m') {
        bulk.cell.type <- 'cux2-'
      } else if (cell.type == 'oligo') {
        bulk.cell.type <- 'oligodendrocytes'
      } else {
        bulk.cell.type <- cell.type
      }

      features <- pseudobulks[[bulk.cell.type]] %>% dplyr::select(de %>% pull(gene)) %>%
        filter(rownames(.) %in% rownames(data)) %>%
        arrange(match(rownames(.), rownames(data))) %>% scale


      features <- data.frame(pathways.score=features %>% rowMeans()) %>% dplyr::rename(!!glue::glue("{current_pathways[[i]][[2]]} - {cell.type}"):=pathways.score)
    }) %>% as.data.frame(check.names=FALSE)

    pretty_colnames <- pathways_scores %>% colnames() %>% gsub("cux2m", "Exc l4-6", .) %>% gsub("cux2p", "Exc L2-3", .)
    colnames(pathways_scores) <- pretty_colnames

    dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                             features = pathways_scores,
                             trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                             trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                             evaluate.fit = T,
                             bootstrap = F)


    withr::with_pdf(file.path(figures_dir, glue::glue("6c_{group}.pdf")), width=embed.width, height=embed.height, {
      print(plot.dynamics.wrapper(dynamics,
                                  features = pathways_scores %>% colnames(),
                                  cols=neuronal_topic_to_color(c("inhibitory", "cux2p", "cux2m")),
                                  labels=labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=1, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=F, include.points=FALSE, trajectories='prAD',
                                  show.label.legend=T, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.text = element_text(margin = margin(0,0,0,0), size = 12),
                    legend.position="bottom", legend.box="vertical", legend.margin=margin(0, 0, 0, 0),
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="Mean pathway expresssion", title=glue::glue("{group} dynamics")))


      print(plot.dynamics.wrapper(dynamics,
                                  features = pathways_scores %>% colnames(),
                                  cols=neuronal_topic_to_color(c("inhibitory", "cux2p", "cux2m")),
                                  labels=labels,
                                  overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                  ncol=2, strip.position="left", scales="free_x",
                                  ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=FALSE, include.points=FALSE,
                                  show.label.legend=TRUE, show_n=FALSE) +
              theme(strip.text = element_blank(),
                    legend.text = element_text(margin = margin(0,0,0,0), size = 12),
                    legend.position="bottom", legend.box="vertical", legend.margin=margin(0, 0, 0, 0),
                    axis.line = element_line(),
                    axis.text.y = element_text()) +
              labs(x="pseudotime" , y="Mean pathway expresssion", title=glue::glue("{group} dynamics")))

      plot.genes.dynamics.heatmap(dynamics,
                                  genes = (pathways_scores %>% colnames()),
                                  scale = FALSE, show_lfc = FALSE, cluster_rows=T, both_trajectories=F) %>% draw
    })
  }
})()

# Figure 7 prAD VS ABA --------

## Figure 7a - neuronal dynamics in the ABA

(function () {
  group_names <- c("decreasing_early", "decreasing_late", "increasing_early", "increasing_late")
  withr::with_pdf(file.path(figures_dir, glue::glue("7a_dynamics_prad_aba.pdf")), width = embed.width, height = embed.height.small, {
    for (group_name in group_names) {
        print(plot.dynamics.wrapper(data$uns$trajectories$palantir$dynamics,
                                    features = names(neuronal_topics_to_annotation[[group_name]]),
                                    labels=neuronal_topics_to_annotation[[group_name]],
                                    cols=neuronal_topic_to_color(names(neuronal_topics_to_annotation[[group_name]])),
                                    overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                    ncol=2, strip.position="left", scales="free_x", ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=FALSE, include.points=FALSE , show.label.legend=T) +
                theme(legend.position="bottom",
                      axis.line = element_line(),
                      axis.text.y = element_text()) +
                labs(x="pseudotime" , y="sqrt(Program abundance)", title=glue::glue("{group_name}")))

      plot.genes.dynamics.heatmap(data$uns$trajectories$palantir$dynamics,
                                  genes = names(neuronal_topics_to_annotation[[group_name]]),
                                  title="@{x}", row_names_mapping_func = pryr::partial(topic_names_mapping, use_static_mapping=FALSE),
                                  both_trajectories = TRUE, show_lfc = F, name = 'Scaled',
                                  clustering_distance_rows='spearman',
                                  use_slanter_ordering=TRUE
      ) %>% draw(column_title='pseudotime',
                 column_title_side='bottom')
    }
  })
})()

## Extended Data Figure 10d - neuronal dynamics in the ABA in green et al trajectories

(function () {
  group_names <- c("decreasing_early", "decreasing_late", "increasing_early", "increasing_late")
  withr::with_pdf(file.path(figures_dir, glue::glue("10SuplD_green_dynamics_prad_aba.pdf")), width = embed.width * 1.2, height = embed.height.small * 1.2, {
    for (group_name in group_names) {
      features <- names(neuronal_topics_to_annotation[[group_name]])
      labels <- topic_names_mapping(features, use_static_mapping = FALSE)
      names(labels) <- features

      print(plot.dynamics.wrapper(data$uns$trajectories$da_fits_sqrt,
                                  features = features,
                                  labels=labels,
                                  cols=neuronal_topic_to_color(features),
                                  overlap.pseudotime=GREEN_ET_AL_OVERLLAPING_PSEUDOTIME,
                                  ncol=2,
                                  scales="free_x", ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  legend.position = "none",
                                  #label=FALSE,
                                  include.points=FALSE , show.label.legend=T) +
              labs(x="Green et al pseudotime" , y="sqrt(Program abundance)", title=glue::glue("{group_name}")))

      print(plot.dynamics.wrapper(data$uns$trajectories$da_fits_sqrt,
                                  features = features,
                                  labels=labels,
                                  cols=neuronal_topic_to_color(features),
                                  overlap.pseudotime=GREEN_ET_AL_OVERLLAPING_PSEUDOTIME,
                                  ncol=2,
                                  scales="free_x", ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                  label=FALSE,
                                  include.points=FALSE , show.label.legend=T) +
              theme(legend.position = "bottom") +
              labs(x="Green et al pseudotime" , y="sqrt(Program abundance)", title=glue::glue("{group_name}")))

      plot.genes.dynamics.heatmap(data$uns$trajectories$da_fits_sqrt,
                                  genes = names(neuronal_topics_to_annotation[[group_name]]),
                                  title="@{x}", row_names_mapping_func = pryr::partial(topic_names_mapping, use_static_mapping=FALSE),
                                  both_trajectories = TRUE, show_lfc = F, name = 'Scaled',
                                  clustering_distance_rows='spearman',
                                  use_slanter_ordering=TRUE
      ) %>% draw(column_title='Green et al pseudotime',
                 column_title_side='bottom')
    }
  })
})()


#+++++++++++++++++++++++++++++++++
## Figure 6d ELISA dynamics -----
#+++++++++++++++++++++++++++++++++

(function() {
  synapse_trait_col_to_name <- list(synap_3cort_vamp='VAMP',
                                    synap_3cort_syntaxin='SYNTAXIN1',
                                    synap_3cort_synaptophys='SYNAPTOPHYSIN',
                                    synap_3cort_stagmin='SYNAPTOTAGMIN',
                                    synap_3cort_snap25='SNAP25',
                                    synap_3cort_sept5='SEPTIN5',
                                    synap_3cort_complex2='COMPLEXIN II',
                                    synap_3cort_complex1='COMPLEXIN I',
                                    synap_3cort_capture4='SNAP25 - VAMP',
                                    synap_3cort_capture3='SNAP25 - SYNTAXIN1',
                                    synap_3cort_capture2='SYNTAXIN1 - SNAP25',
                                    synap_3cort_capture1='SYNTAXIN1 - VAMP',
                                    zcapture_syn_3cort='Total interaction')

  features_groups <- list("Interaction"=c("synap_3cort_capture1", "synap_3cort_capture2", "synap_3cort_capture3", "synap_3cort_capture4", "zcapture_syn_3cort"),
                          "Proteins"=c("synap_3cort_vamp", "synap_3cort_snap25",
                                       "synap_3cort_sept5", "synap_3cort_synaptophys",
                                        "synap_3cort_stagmin"),
                          "Complexin"=c("synap_3cort_complex2", "synap_3cort_complex1"))
  trajectories_list <- list(c('prAD'), c('prAD', 'ABA'))
  withr::with_pdf(file.path(figures_dir, "6d_elisa_dynamics.pdf"), height = embed.height, width = embed.width * 1.6, {
    for (group in names(features_groups)) {
      for (trajectories in trajectories_list) {
        print(plot.dynamics.wrapper(data$uns$trajectories$palantir$syn_dynamics,
                                    features = features_groups[[group]],
                                    labels=synapse_trait_col_to_name,
                                    overlap.pseudotime=TOPIC_OVERLLAPING_PSEUDOTIME,
                                    ncol=2, strip.position="left", scales="free_x",
                                    ribbon_alpha=0.1, line_size=0.8, label_size=6,
                                    label=T, include.points=FALSE, trajectories = trajectories,
                                    show.label.legend=T, show_n=FALSE) +
                theme(strip.text = element_blank(),
                      legend.position="none",
                      axis.line = element_line(),
                      axis.text.y = element_text()) +
                labs(x="pseudotime" , y="scaled ELISA", title=glue::glue("ELISA {group} dynamics")))
      }
    }
  })
})()
