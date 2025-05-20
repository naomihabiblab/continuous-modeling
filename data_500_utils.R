library(stringr)
library(dplyr)
library(anndata)
library(tibble)
library(rhdf5)
library(reshape2)

source("utils.R")
source("pathways.analysis.R")

H5_500_PATH <- file.path(get_result_dir(), "../../data/500_29_8_23.h5")
AD_500_ANDATA_FILE_V1 <- "../data/500_19_12_2022.h5ad"
AD_500_ANDATA_FILE_V1.1 <- "../data/500_26_08_2023.h5ad"
AD_500_ANDATA_FILE <- AD_500_ANDATA_FILE_V1.1

cell_type_to_group <- function (cell.type) {
  group <- switch(cell.type,
                  oligodendrocytes="glia/oligodendroglia",
                  opcs="glia/oligodendroglia",
                  microglia="glia/microglia",
                  astrocytes="glia/astrocytes",
                  inhibitory="neuronal/inhibitory",
                  endothelial="vascular.niche",
                  smc="vascular.niche",
                  fibroblast="vascular.niche",
                  pericytes="vascular.niche",
                  pericytes_fibroblast_smc="vascular.niche",
                  vascular.niche="vascular.niche",
                  oligodendroglia="glia/oligodendroglia",
                  excitatory="neuronal/excitatory",
                  `excitatory_cux-`="neuronal/excitatory",
                  `excitatory_cux+`="neuronal/excitatory",
                  stop("Not handled cell_type"))
  filter_function <-  switch(cell.type,
                             endothelial=function (df) {df %>% filter(str_detect(state, '^End.|Venule|Arteriole'))},
                             fibroblast=function (df) {df %>% filter(str_detect(state, '^Fib.'))},
                             smc=function (df) {df %>% filter(str_detect(state, '^SMC.'))},
                             pericyte=function (df) {df %>% filter(str_detect(state, '^Peri.'))},
                             pericytes_fibroblast_smc=function (df) {df %>% filter(str_detect(state, '^Peri.|^SMC.|^Fib.'))},
                             `excitatory_cux+`=function (df) {df %>% filter(str_detect(state, "Exc.10|Exc.11|Exc.12|Exc.13|Exc.14|Exc.15|Exc.16||Exc.8|Exc.9"))},
                             `excitatory_cux-`=function (df) {df %>% filter(str_detect(state, "Exc.1|Exc.2|Exc.3|Exc.4|Exc.5|Exc.6|Exc.7"))},
                             opcs=function (df) {df %>% filter(!grepl("Oli.", state))},
                             oligodendrocytes=function (df) {df %>% filter(grepl("Oli.", state))},
                             function (df) df)

  list(group=group, filter_function=filter_function)
}

read_expressing_cells_stats <- function (cell.type, data_path=H5_500_PATH) {
  group_and_filter <- cell_type_to_group(cell.type)
  gene.exp <- h5read(data_path, paste0(group_and_filter$group, '/gene.exp')) %>% group_and_filter$filter_function()

  number_of_cell <- gene.exp %>%
    dplyr::group_by(state) %>%
    dplyr::slice(1) %>%
    dplyr::pull(n) %>%
    sum()

  gene.exp %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(n.exp=sum(n.exp)) %>%
    dplyr::mutate(p.exp = n.exp / number_of_cell)
}

read_cell_dataset <- function (cell.type, dataset_name, data_path=H5_500_PATH) {
  group_and_filter <- cell_type_to_group(cell.type)
  gene.exp <- h5read(data_path, paste0(group_and_filter$group, str_glue('/{dataset_name}'))) %>% group_and_filter$filter_function()
}

read_states_de <- function (cell.type, ...) {
  read_cell_dataset(cell.type, dataset_name = 'de', ...)
}

read_states_pathways <- function (cell.type, ...) {
  read_cell_dataset(cell.type, dataset_name = 'pa', ...)
}

# Fixing up an object

fix_similarity_names <- function (object) {
  for(s in c("sim.umap", "sim.phate")) {
    if (s %in% names(Misc(object))) {colnames(object@misc[[s]]) = rownames(object@misc[[s]]) = colnames(object) }
  }
  object
}

enrich_meta_data <- function (object) {
  for (col_name in c('tangles', 'amyloid')) {
    new_colname <- paste0('sqrt.', col_name)
    if (!new_colname %in% colnames(object@meta.data)) {
      object <- AddMetaData(object, sqrt(FetchData(object, col_name)), paste0('sqrt.', col_name))
    }
  }

  if (!'cogng_demog_slope' %in% colnames(object@meta.data)) {
    data.da <- read_h5ad(AD_500_ANDATA_FILE)
    cogng_demog_slope <- dplyr::left_join(object@meta.data %>% dplyr::select(projid) %>%  mutate(projid=as.character(projid)),
                   data.da$obsm$meta.data %>% dplyr::select(cogng_demog_slope) %>% rownames_to_column(var = "projid"), by='projid') %>% pull(cogng_demog_slope)

    object <- AddMetaData(object = object, cogng_demog_slope, 'cogng_demog_slope')
  }

  if (!is.factor(object@meta.data$sex)) {
    object@meta.data$sex <- factor(object@meta.data$sex, levels = c("Female", "Male"))
  }

  object
}


TOPICS_DATA_H5_PATH <- file.path(get_result_dir(), "topics_data.h5")


save_abundance_stat <- function (obj, group, topics_columns=topics_columns_list[[k]], h5_path=TOPICS_DATA_H5_PATH) {
  message(paste("writing summary output to", group, "at path", h5_path))

  h5createGroup(h5_path, group = group)

  df <- obj@meta.data %>%
    dplyr::select(state, topics_columns) %>%
    reshape2::melt(id.vars="state") %>%
    group_by(state, variable) %>%
    summarise(avg=mean(value), `pct.abv.05`=mean(value>0.05), n=n())

  h5write(df, h5_path, paste0(group, "/abundance.state"))

  df <- obj@meta.data %>%
    dplyr::select(projid, topics_columns) %>%
    reshape2::melt(id.vars="projid") %>%
    group_by(projid, variable) %>%
    summarise(avg=mean(value), `pct.abv.05`=mean(value>0.05), n=n())

  h5write(df, h5_path, paste0(group, "/abundance.projid"))

  df <- obj@meta.data %>%
    dplyr::select(projid, state, topics_columns) %>%
    reshape2::melt(id.vars=c("projid", "state")) %>%
    group_by(projid, state, variable) %>%
    summarise(avg=mean(value), `pct.abv.05`=mean(value>0.05), n=n())

  h5write(df, h5_path,paste0(group, "/abundance.projid.state"))
}

read_object_metadata <- function(cell.type, field, objects_dir=SEURAT_OBJECT_DATA_PATH) {
  cell.type.name_mapping <- list('cux2m'='cux2-', 'cux2p'='cux2+')

  if (cell.type %in% names(cell.type.name_mapping)) {
    cell.type <- cell.type.name_mapping[[cell.type]]
  }

  meta.data <- h5read(file.path(objects_dir, paste0(cell.type, '.h5Seurat')), paste0('meta.data/', field))
  cell.names <- h5read(file.path(objects_dir, paste0(cell.type, '.h5Seurat')), 'cell.names')
  names(meta.data) <- cell.names

  meta.data
}

save_topic_per_neuronal_label <- function (cell.types=c('cux2m', 'cux2p', 'inhibitory') , h5_path=TOPICS_DATA_H5_PATH, verbose=FALSE) {
  for (ct in cell.types) {
    h5_group <- CELL_TYPE_TO_GROUP[[ct]]

    if (verbose) message(glue::glue("Writing allen label abundance to {h5_group} for cell type {ct}"))

    L <- h5read(h5_path, paste0(h5_group, "/L"))
    cell.type.name <- list('cux2m'='cux-',
                           'cux2p'='cux+',
                           'inhibitory'='inh'
                          )[[ct]]
    if (verbose) message("Reading labels")
    labels <- readRDS(file.path(get_result_dir(), glue::glue('../../data/preds_{cell.type.name}_subclass.rds')))
    class <- readRDS(file.path(get_result_dir(), glue::glue('../../data/preds_{cell.type.name}.rds')))

    df <- L %>%
      column_to_rownames("cell")

    df$label <- labels[rownames(df), 'pruned.labels']
    df$class <- str_to_upper(class[rownames(df), 'pruned.labels'])
    df$marker.1 <- sapply(str_split(df$label, " "), "[", 3)

    projid <- read_object_metadata(ct, field = "projid")
    df$projid <- projid[rownames(df)]

    if (verbose) message("Computing stats")

    group_terms <- c("label", "marker.1", "class")
    res <- sapply(group_terms, function(group.by) {
      df %>%
        filter(!is.na(get(group.by))) %>%
        group_by(get(group.by), projid) %>%
        summarise(n=n(), across(starts_with('k'), list(mean), .names = "{.col}")) %>%
        rename(!!group.by := `get(group.by)`) %>%
        reshape2::melt(id.vars=c(group.by, "n", 'projid'), variable.name = "topic", value.name = "avg")
    }, simplify = F)

    h5write(res$label, h5_path, paste0(h5_group, "/abundance.allen.projid"))
    h5write(res$marker.1, h5_path, paste0(h5_group, "/abundance.allen.projid.marker.1"))
    h5write(res$class , h5_path, paste0(h5_group, "/abundance.allen.projid.class"))

    df.labels <- df %>%
      filter(!is.na(label)) %>%
      group_by(label) %>%
      summarise(n=n(), across(starts_with('k'), list(mean), .names = "{.col}")) %>%
      melt(id.vars=c("label", "n"), variable.name = "topic", value.name = "avg")

    df.marker.1 <- df %>%
      filter(!is.na(marker.1)) %>%
      group_by(marker.1) %>%
      summarise(n=n(), across(starts_with('k'), list(mean), .names = "{.col}")) %>%
      melt(id.vars=c("marker.1", "n"), variable.name = "topic", value.name = "avg")

    df.class <- df %>%
      filter(!is.na(class)) %>%
      group_by(class) %>%
      summarise(n=n(), across(starts_with('k'), list(mean), .names = "{.col}")) %>%
      melt(id.vars=c("class", "n"), variable.name = "topic", value.name = "avg")

    if (verbose) message("Writing to h5")
    h5write(df.labels, h5_path, paste0(h5_group, "/abundance.allen"))
    h5write(df.marker.1, h5_path, paste0(h5_group, "/abundance.allen.marker.1"))
    h5write(df.class, h5_path, paste0(h5_group, "/abundance.allen.class"))

  }
}

CELL_TYPE_TO_GROUP <- list(
  "inhibitory"="neuronal/inhibitory",
  "oligo"="glia/oligodendrocytes",
  "microglia"="glia/microglia",
  "opcs"="glia/opcs",
  "cux2p"="neuronal/excitatory/cux2plus",
  "cux2m"="neuronal/excitatory/cux2minus",
  "astrocytes"="glia/astrocytes"
)

save_umap_densities <- function(obj, columns, h5_path=TOPICS_DATA_H5_PATH, cell.type=obj@project.name) {
  obj <- FeatureDensity(obj, features=columns, return_object = TRUE)
  umap_densities <- obj@meta.data %>%
    dplyr::select(starts_with(columns) & ends_with("sim.umap"), state)

  umap_densities <- cbind(umap_densities, obj@reductions$umap@cell.embeddings) %>%
    rownames_to_column('cell')
  h5delete(h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/umap.densities"))
  h5write(umap_densities, h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/umap.densities"))
}


save_allen_annotation <- function(obj, h5_path=TOPICS_DATA_H5_PATH, cell.type=obj@project.name) {
  classes_level <- obj@meta.data$pruned_subclass %>% levels %>% as.data.frame()

  df <- obj@meta.data %>%
    dplyr::select(pruned_subtype, pruned_subclass) %>%
    mutate(pruned_subclass=as.character(pruned_subclass)) %>%
    rownames_to_column("cell")
  h5write(df, h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/allen.annotation"))
  h5write(classes_level, h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/allen.annotation.level"))
}


read.dataset.for.topics <- function(topics, dataset_path) {
    df <- lapply(topics, function(topic_name) {
      k <- str_extract(topic_name, "k(.*)")
      ct <- case_when(str_detect(topic_name, "Inh")  ~ "inhibitory",
                str_detect(topic_name, "Exc.Cux2P") ~ "cux2p",
                str_detect(topic_name, "Exc.Cux2M") ~ "cux2m",
                TRUE ~ gsub("\\.k.*", "", topic_name))

      df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], dataset_path))
      df <- df %>% filter(topic == k) %>% mutate(topic_name=topic_name)

      df
  }) %>% do.call(rbind, .)

  df
}

neuronal_cell_type_to_index <- function (cell.type) {
  switch(cell.type, "inh"=1, "inhibitory"=1, "cux2p"=2, "cux2+"=2, "cux2m"=3, "cux2-"=3, stop(glue::glue("unkown cell type - {cell.type }")))
}

get_AD_programs_proportion <- function(ct) {
  df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))

  prop_cells_df <- df %>%
        mutate(group=factor(group, levels=c('start', 'prAD'))) %>%
        group_by(projid, group) %>%
        summarise(prop_cell_cycle_95=mean(cell_cycle_95),
                  prop_cell_cycle_90=mean(cell_cycle_90),
                  prop_cell_cycle_80=mean(cell_cycle_80),
                  prop_apoptosis_95=mean(apoptosis_95),
                  prop_apoptosis_90=mean(apoptosis_90),
                  prop_apoptosis_80=mean(apoptosis_80)
      )

  prop_cells_df
}


compute_signatures_score <- function(obj, cell.type=obj@project.name) {
  get_topic_de_for_pathway <- function(de, topic, pathway_id) {
    pathway_genes <- get.pathway.genes(pathway_id)
    de %>% filter(.data[['topic_name']] == .env[['topic']] & gene %in% pathway_genes) %>% pull(gene)
  }
  cell_type_index <- neuronal_cell_type_to_index(cell.type)
  postsynaptic_topic <- names(neuronal_topics_to_annotation$decreasing_late)[[cell_type_index]]
  presynaptic_topic <- names(neuronal_topics_to_annotation$increasing_early)[[cell_type_index]]
  des <- read.dataset.for.topics(c(postsynaptic_topic, presynaptic_topic), dataset_path = "/de")

  signatures <- list("presynaptic_genes"=c(get_topic_de_for_pathway(des, presynaptic_topic, "GO:0098693"), # regulation of synaptic vesicle cycle
                                           get_topic_de_for_pathway(des, presynaptic_topic, "GO:0099504"), #synaptic vesicle cycle
                                           get_topic_de_for_pathway(des, presynaptic_topic, "hsa04721")), # Synaptic vesicle cycle
                     "postsynaptic_genes"=c(get_topic_de_for_pathway(des, postsynaptic_topic, "GO:0099572")),
                     "ieg_genes"=c("FOS", "EGR1", "EGR2", "EGR3", "ARC", "PROX1", "FOSB", "JUNB", "HOMER1"),
                     MPhase = str_split("CLASP2/CENPC/SEM1/TPR/STAG2/NSL1/AHCTF1/ANAPC16/SMC3/ITGB3BP/BLZF1/CEP290/OFD1/CEP70/ANAPC4/CEP63/AKAP9/LEMD3/NUP205/PCM1/WAPL/CCP110", "/")[[1]],
                     NeuronalApoptotic = str_split("NEFL/MT3/APOE/GAPDH/NSMF/PTK2B/SNCB/UBE2M/PIN1/GRIN1/AGAP2/TYRO3/FIS1/NUPR1/BNIP3", "/")[[1]],
                     ResponseOxidativeStress=str_split("MT3/CST3/APOE/PRDX5/EEF2/HSPB1/PTK2B/PON2/GPX4/RPS3/AGAP3/MAPK3/GPR37L1/IPCEF1/PRDX1/TXNIP/APOD/CRYAB/CAMKK2/BNIP3/GPX3/NME2/RHOB/HBB/SELENOP", "/")[[1]],
                     HeatShockGenes=c("ATM", "TPR", "HSPH1", "HSPA4L", "NUP205"),
                     DNARepair=str_split("ATM/POLK/XPA/SMARCA5/BAZ1B/INO80D/COPS2/WRN/REV1/PPP4R2/ERCC5/RIF1", "/")[[1]]
                     )

  dput(names(signatures))

  signature_lengths <- sapply(signatures, length, USE.NAMES = TRUE) %>% paste0(names(.),"=", ., collapse = " ")
  message(glue::glue("Signatures with lengths {signature_lengths}"))

  obj <- AddModuleScore(obj, signatures)

  column_names <- paste0("Cluster", 1:length(signatures))
  names(column_names) <- names(signatures)

  obj@meta.data <- obj@meta.data %>% dplyr::rename(!!!column_names)

  obj
}

save_genes_signatures <- function(obj, cell.type=obj@project.name, h5_path=TOPICS_DATA_H5_PATH,
                                  features=c("presynaptic_genes", "postsynaptic_genes",
                                             "SST", "VIP", "NPY",
                                             "MPhase", "NeuronalApoptotic",
                                             "ResponseOxidativeStress", "HeatShockGenes", "DNARepair")) {
  df <- FetchData(object = obj, c(features, 'projid'))
  df <- df %>% rownames_to_column('cell')
  h5delete(h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/genes.signatures"))
  h5write(df, h5_path, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/genes.signatures"))
}