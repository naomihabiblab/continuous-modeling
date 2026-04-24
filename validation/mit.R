library(ggplot2)
library(Seurat)
library(dplyr)
library(withr)
library(rhdf5)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)


source("R/utils.R")
source("R/topics.R")
source("R/data_500_utils.R")
source("R/figures.funs.R")

source("validation/validation.functions.R")

project <- "mit"
cell_types <- c("inhibitory", "astrocytes", "microglia", "oligo", "opcs")

MIT_PAPER_DATA_DIR <- file.path(get_shared_dir(), "data_Single-cell atlas reveals correlates of high cognitive function, dementia, and resilience to Alzheimers disease pathology")

for (cell.type in cell_types) {
  file_name <- case_when(
    cell.type == "inhibitory" ~ "Inhibitory_neurons",
    cell.type == "excitatory" ~ "Excitatory_neurons",
    cell.type == "astrocytes" ~ "Astrocytes",
    cell.type == "microglia" ~ "Immune_cells",
    cell.type == 'opcs' ~ "OPCs",
    cell.type == 'oligo' ~ "Oligodendrocytes",
  )

  OBJECT_PATH <- file.path(MIT_PAPER_DATA_DIR, glue::glue("{file_name}.rds"))
  obj <- readRDS(OBJECT_PATH)

  if (cell.type == "microglia" && project == "mit") {
    obj <- subset(obj, subset = cell_type_high_resolution %in% c("CAMs", "Mic MKI67", "Mic P2RY12", "Mic TPT1"))
  }

  project_dir <- file.path(get_result_dir(), project)
  plots_dir <- file.path(project_dir, cell.type)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  if (identical(obj[['RNA']]@counts@x, obj[['RNA']]@data@x)) {
    message("The counts and data matrices are the same -> the object is not normalized")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  } else {
    message("The counts and data matrices are different -> the object is normalized")
  }

  donor.meta.data <- read.csv(file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv"))
  if ('subject' %in% colnames(obj@meta.data) && !('projid' %in% colnames(obj@meta.data))) {
    obj@meta.data$projid <- donor.meta.data$projid[match(obj@meta.data$subject, donor.meta.data$subject)]
  }

  genes_mapping <- get_genes_by_topic(cell.type)

  if (cell.type != "opcs") {
    pathway_to_genes <- get_pathway_signatures(cell.type)
  } else {
    pathway_to_genes <- list()
  }


  SHOULD_COMPUTE_COHERENCE <- TRUE
  if (cell.type == "inhibitory" && SHOULD_COMPUTE_COHERENCE) {
    set.seed(1)

    if (cell.type == "inhibitory" && 'k1' %in% names(genes_mapping$unique_genes_by_topic)) {
      programs_to_check <- paste0("k", c(1,7,8,10))
      gene_sets_to_test <- genes_mapping$unique_genes_by_topic[programs_to_check] %>% lapply(c)
    } else {
      gene_sets_to_test <- genes_mapping$unique_genes_by_topic %>% lapply(c)
    }
    expr_matrix <- GetAssayData(obj, slot = "data")

    coherence_results <- compute_coherence_multiple(
      expr_matrix = expr_matrix,
      gene_sets = gene_sets_to_test,
      n_permutations = 100,
      method = "pearson",
      seed = 1
    )

    saveRDS(coherence_results, file.path(plots_dir, "SignatureCoherence_results.rds"))

    h5write(coherence_results %>% as.data.frame(), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/", project, "/coherence_results"))
    h5write(unAsIs(coherence_results$bg_means_cor), TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[cell.type]], "/", project, "/coherence_results_bg_means"))
  }

  if (cell.type != "opcs") {
    SHOULD_COMPUTE_SIGNATURES <- !file.exists(file.path(plots_dir, "Signature.rds")) || !file.exists(file.path(plots_dir, "Pathway.rds"))
  } else {
    SHOULD_COMPUTE_SIGNATURES <- !file.exists(file.path(plots_dir, "Signature.rds"))
  }

  if (SHOULD_COMPUTE_SIGNATURES) {
    signature_score_prefix <- "Signature_"
    obj <- add_module_score_named(obj, genes_mapping$genes_by_topic, prefix=signature_score_prefix)

    signature_unique_score_prefix <- "SignatureUnique_"
    obj <- add_module_score_named(obj, genes_mapping$unique_genes_by_topic, prefix=signature_unique_score_prefix)

    obj@meta.data %>% dplyr::select(starts_with("Signature")) %>% saveRDS(file.path(plots_dir, "Signature.rds"))
    obj@meta.data %>% dplyr::select(starts_with("Signature")) %>% write.csv(file.path(plots_dir, "Signature.csv"))

    obj <- add_module_score_named(obj, pathway_to_genes, prefix="Pathway_")

    obj@meta.data %>% dplyr::select(starts_with("Pathway_")) %>% saveRDS(file.path(plots_dir, "Pathway.rds"))
    obj@meta.data %>% dplyr::select(starts_with("Pathway_")) %>% write.csv(file.path(plots_dir, "Pathway.csv"))
  } else {
    obj@meta.data <- obj@meta.data %>% bind_cols(readRDS(file.path(plots_dir, "Signature.rds")))
    if (file.exists(file.path(plots_dir, "Pathway.rds"))) {
      obj@meta.data <- obj@meta.data %>% bind_cols(readRDS(file.path(plots_dir, "Pathway.rds")))
    }
  }

  withr::with_pdf(file.path(plots_dir, "SignatureUniqueCorrelation.pdf"), {
    obj@meta.data %>% dplyr::select(starts_with("SignatureUnique")) %>% cor %>% ComplexHeatmap::Heatmap() %>% draw()
  })


  if (cell.type == "inhibitory") {
    obj@meta.data <- obj@meta.data %>% mutate(subtype=case_when(
      cell_type_high_resolution == "Inh ALCAM TRPM3" ~ "Inh VIP",
      cell_type_high_resolution == "Inh CUX2 MSR1" ~ "Inh SST",
      cell_type_high_resolution == "Inh ENOX2 SPHKAP" ~ "Inh SST",
      cell_type_high_resolution == "Inh FBN2 EPB41L4A" ~ "Inh SST",
      cell_type_high_resolution == "Inh GPC5 RIT2" ~ "Inh PVALB",
      cell_type_high_resolution == "Inh L1 PAX6 CA4" ~ "Inh PAX6",
      cell_type_high_resolution == "Inh L1-2 PAX6 SCGN" ~ "Inh PAX6",
      cell_type_high_resolution == "Inh L1-6 LAMP5 CA13" ~ "Inh LAMP5",
      cell_type_high_resolution == "Inh L3-5 SST MAFB" ~ "Inh SST",
      cell_type_high_resolution == "Inh L5-6 PVALB STON2" ~ "Inh PVALB",
      cell_type_high_resolution == "Inh L5-6 SST TH" ~ "Inh SST",
      cell_type_high_resolution == "Inh L6 SST NPY" ~ "Inh SST",
      cell_type_high_resolution == "Inh LAMP5 NRG1 (Rosehip)" ~ "Inh LAMP5",
      cell_type_high_resolution == "Inh LAMP5 RELN" ~ "Inh LAMP5",
      cell_type_high_resolution == "Inh PTPRK FAM19A1" ~ "Inh VIP",
      cell_type_high_resolution == "Inh PVALB CA8 (Chandelier)" ~ "Inh PVALB Chandelier",
      cell_type_high_resolution == "Inh PVALB HTR4" ~ "Inh PVALB",
      cell_type_high_resolution == "Inh PVALB SULF1" ~ "Inh PVALB",
      cell_type_high_resolution == "Inh RYR3 TSHZ2" ~ "Inh VIP",
      cell_type_high_resolution == "Inh SGCD PDE3A" ~ "Inh VIP",
      cell_type_high_resolution == "Inh SORCS1 TTN" ~ "Inh VIP",
      cell_type_high_resolution == "Inh VIP ABI3BP" ~ "Inh VIP",
      cell_type_high_resolution == "Inh VIP CLSTN2" ~ "Inh VIP",
      cell_type_high_resolution == "Inh VIP THSD7B" ~ "Inh VIP",
      cell_type_high_resolution == "Inh VIP TSHZ2" ~ "Inh VIP",
      TRUE ~ NA_character_
    ))
  }

  h5createGroup(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project))

  subtype_column <- if (cell.type == "inhibitory") "subtype" else NULL
  h5_save_mean_scores(obj@meta.data, cell.type, project, group.by='projid', subtype_column=subtype_column)

  c <- correlate_dataframes(obj@meta.data %>% dplyr::select(starts_with("SignatureUnique_")))
  h5_save_correlation(c, cell.type, paste0(project, "/signature_unique_correlation"))

  if (obj@meta.data %>% dplyr::select(starts_with("Pathway_")) %>% ncol() > 0) {
    c <- correlate_dataframes(obj@meta.data %>% dplyr::select(starts_with("Pathway_")))
    h5_save_correlation(c, cell.type, paste0(project, "/pathway_correlation"))
  }
}

# Add excitatory pathways and signatures from python notebook

cell.type <- "excitatory"
pathways_and_signatures <- read.csv(file.path(get_result_dir(), "mit", cell.type, "pathways_and_signatures.csv"), row.names=1)

if ('subject' %in% colnames(pathways_and_signatures) && !('projid' %in% colnames(pathways_and_signatures))) {
    pathways_and_signatures$projid <- donor.meta.data$projid[match(pathways_and_signatures$subject, donor.meta.data$subject)]
}
# Change pathway names with . to _
colnames(pathways_and_signatures) <- gsub("\\.", "_", colnames(pathways_and_signatures))


h5_save_mean_scores(pathways_and_signatures, cell.type, project, group.by='projid', subtype_column=NULL)

c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("SignatureUnique_")))
h5_save_correlation(c, cell.type, paste0(project, "/signature_unique_correlation"))

if (pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")) %>% ncol() > 0) {
  c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")))
  h5_save_correlation(c, cell.type, paste0(project, "/pathway_correlation"))
}

# Trait association -----
donor.meta.data <- read.csv(file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv")) %>% column_to_rownames("projid")
traits_columns <- c("sqrt.amyloid", "sqrt.amyloid_mf", "sqrt.tangles", "sqrt.tangles_mf", "cogng_demog_slope")

for (neuronal_class in c("inhibitory", "excitatory")) {
    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/mean_scores")) %>% filter(!is.na(projid)) %>% column_to_rownames("projid")

    if (neuronal_class == "inhibitory") {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"), 
        SignatureUnique_k1, SignatureUnique_k7, SignatureUnique_k8, SignatureUnique_k10)
    } else {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"))
    }

    traits <- donor.meta.data %>% dplyr::select(all_of(traits_columns))
    traits <- traits[rownames(df),]
    controls <- donor.meta.data %>% dplyr::select(age_death, msex, pmi) 
    controls <- controls[rownames(df),]
    trait.association <- associate.traits(traits, df_signature_unique, controls)
    trait.association_pathway <- associate.traits(traits, df %>% dplyr::select(starts_with("Pathway_")), controls)

    h5write(trait.association, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_signature"))
    h5write(trait.association_pathway, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_pathway"))
}

# Community analysis -----

compute_program_environment_correlation(project, identity_column="projid")
