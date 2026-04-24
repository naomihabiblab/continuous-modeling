library(ggplot2)
library(Seurat)
library(dplyr)
library(withr)
library(rhdf5)
library(ggplot2)
library(ggpointdensity)


source("R/utils.R")
source("R/topics.R")
source("R/data_500_utils.R")
source("R/figures.funs.R")

source("validation/validation.functions.R")

project <- "cuimc2"
cell_types <- c("inhibitory", "excitatory", "astrocytes", "microglia")

CUIMC2_OBJECT_DIR <- file.path(get_shared_dir(), "DLPFCSingleCellFromNatasha")

for (cell.type in cell_types) {
  OBJECT_PATH <- file.path(CUIMC2_OBJECT_DIR, project, glue::glue("{cell.type}.rds"))
  if (!file.exists(OBJECT_PATH)) next
  obj <- readRDS(OBJECT_PATH)

  project_dir <- file.path(get_result_dir(), project)
  plots_dir <- file.path(project_dir, cell.type)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


    if (identical(obj[['RNA']]@counts@x, obj[['RNA']]@data@x)) {
        message("The counts and data matrices are the same")
        obj <- NormalizeData(obj, normalization.method = "LogNormalize")
    } else {
        message("The counts and data matrices are different i.e. the object is normalized")
    }

    donor.meta.data <- read.csv(file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv"))
    if ('subject' %in% colnames(obj@meta.data) && !('projid' %in% colnames(obj@meta.data))) {
        obj@meta.data$projid <- donor.meta.data$projid[match(obj@meta.data$subject, donor.meta.data$subject)]
    }


    genes_mapping <- get_genes_by_topic(cell.type)
    pathway_to_genes <- get_pathway_signatures(cell.type)


    SHOULD_COMPUTE_COHERENCE <- TRUE
    if (cell.type %in% c("inhibitory", "excitatory") && SHOULD_COMPUTE_COHERENCE) {
        set.seed(1)

        if (cell.type == "inhibitory") {
            programs_to_check <- paste0("k", c(1,7,8,10))
            gene_sets_to_test <- genes_mapping$unique_genes_by_topic[programs_to_check] %>% lapply(c)
        } else {
            gene_sets_to_test <- genes_mapping$unique_genes_by_topic %>% lapply(c)
        }
        # Get expression matrix once (memory efficient)
        expr_matrix <- GetAssayData(obj, slot = "data")

        # Compute coherence for all programs
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

    SHOULD_COMPUTE_SIGNATURES <- !file.exists(file.path(plots_dir, "Signature.rds")) || !file.exists(file.path(plots_dir, "Pathway.rds"))
    if (SHOULD_COMPUTE_SIGNATURES) {
        signature_score_prefix <- "Signature_"
        obj <- add_module_score_named(obj, genes_mapping$genes_by_topic, prefix=signature_score_prefix)
        # obj <- add_module_score_named(obj, genes_mapping$genes_by_topic_top_50, prefix="Signature_Top50_")

        signature_unique_score_prefix <- "SignatureUnique_"
        obj <- add_module_score_named(obj, genes_mapping$unique_genes_by_topic, prefix=signature_unique_score_prefix)
        
        obj@meta.data %>% dplyr::select(starts_with("Signature")) %>% saveRDS(file.path(plots_dir, "Signature.rds"))
        obj@meta.data %>% dplyr::select(starts_with("Signature")) %>% write.csv(file.path(plots_dir, "Signature.csv"))


        obj <- add_module_score_named(obj, pathway_to_genes, prefix="Pathway_")

        obj@meta.data %>% dplyr::select(starts_with("Pathway_")) %>% saveRDS(file.path(plots_dir, "Pathway.rds"))
        obj@meta.data %>% dplyr::select(starts_with("Pathway_")) %>% write.csv(file.path(plots_dir, "Pathway.csv"))
    } else {
        obj@meta.data <- obj@meta.data %>% bind_cols(readRDS(file.path(plots_dir, "Signature.rds")))
        obj@meta.data <- obj@meta.data %>% bind_cols(readRDS(file.path(plots_dir, "Pathway.rds")))
    }

    withr::with_pdf(file.path(plots_dir, "SignatureUniqueCorrelation.pdf"), {
        obj@meta.data %>% dplyr::select(starts_with("SignatureUnique")) %>% cor %>% ComplexHeatmap::Heatmap() %>% draw()
    })


    if (cell.type == "inhibitory") {
        obj@meta.data$subtype <- case_when(
            obj@meta.data$predicted.state == "Inh.1" ~ "Inh SST",
            obj@meta.data$predicted.state == "Inh.2" ~ "Inh LAMP5",
            obj@meta.data$predicted.state == "Inh.3" ~ "Inh LAMP5",
            obj@meta.data$predicted.state == "Inh.4" ~ "Inh PAX6",
            obj@meta.data$predicted.state == "Inh.5" ~ "Inh SST",
            obj@meta.data$predicted.state == "Inh.6" ~ "Inh SST",
            obj@meta.data$predicted.state == "Inh.7" ~ "Inh SST",
            obj@meta.data$predicted.state == "Inh.8" ~ "SNCG",
            obj@meta.data$predicted.state == "Inh.9" ~ "Inh VIP",
            obj@meta.data$predicted.state == "Inh.10" ~ "Inh VIP",
            obj@meta.data$predicted.state == "Inh.11" ~ "Inh VIP",
            obj@meta.data$predicted.state == "Inh.12" ~ "Inh PVALB Chandelier",
            obj@meta.data$predicted.state == "Inh.13" ~ "Inh PVALB",
            obj@meta.data$predicted.state == "Inh.14" ~ "Inh PVALB",
            obj@meta.data$predicted.state == "Inh.15" ~ "Inh PVALB",
            obj@meta.data$predicted.state == "Inh.16" ~ "Inh PVALB",
            TRUE ~ NA_character_
        )
    }


    h5createGroup(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project))

    subtype_column <- if (cell.type == "inhibitory") "subtype" else NULL
    h5_save_mean_scores(obj@meta.data, cell.type, project, group.by='projid', subtype_column=subtype_column)

    c <- correlate_dataframes(obj@meta.data %>% dplyr::select(starts_with("SignatureUnique_")))
    h5_save_correlation(c, cell.type, paste0(project, "/signature_unique_correlation"))

    c <- correlate_dataframes(obj@meta.data %>% dplyr::select(starts_with("Pathway_")))
    h5_save_correlation(c, cell.type, paste0(project, "/pathway_correlation"))
}

# Trait association -----
project_dir <- file.path(get_result_dir(), project)
donor.meta.data <- read.csv(file.path(get_shared_dir(), "500_topics", "donor_meta_data.csv")) %>% column_to_rownames("projid")
traits_columns <- c("sqrt.amyloid", "sqrt.amyloid_mf", "sqrt.tangles", "sqrt.tangles_mf", "cogng_demog_slope")

for (neuronal_class in c("inhibitory", "excitatory")[2]) {
  df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/mean_scores")) %>% column_to_rownames("projid")

    if (neuronal_class == "inhibitory") {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"), 
        SignatureUnique_k1, SignatureUnique_k7, SignatureUnique_k8, SignatureUnique_k10)
    } else {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_"), starts_with("SignatureUnique_decreasing"))
    }

    traits <- donor.meta.data %>% dplyr::select(all_of(traits_columns))
    traits <- traits[rownames(df),]
    controls <- donor.meta.data %>% dplyr::select(age_death, msex, pmi) 
    controls <- controls[rownames(df),]
    trait.association <- associate.traits(traits, df_signature_unique, controls)
    trait.association_pathway <- associate.traits(traits, df %>% dplyr::select(starts_with("Pathway_")), controls)

  # h5write(trait.association, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_signature"))
  h5write(trait.association_pathway, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_pathway"))
}

# Community analysis -----

compute_program_environment_correlation(project, identity_column="projid")
