library(ggplot2)
library(Seurat)
library(dplyr)
library(withr)
library(rhdf5)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(stringr)

source("R/utils.R")
source("R/topics.R")
source("R/data_500_utils.R")
source("R/figures.funs.R")

source("validation/validation.functions.R")

project <- "psychAD"
cell_types <- c("inhibitory", "excitatory", "astrocytes", "microglia", "oligo", "opcs")

st_github <- read.csv(file.path(get_shared_dir(), "glia_diversity/data/PsychAD_SupplementaryTable1.csv"), check.names = FALSE)
project_dir <- file.path(get_result_dir(), project)

for (cell.type in cell_types) {
  plots_dir <- file.path(project_dir, cell.type)
  if (!file.exists(file.path(plots_dir, "pathways_and_signatures.csv.gz"))) next
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  donor.meta.data <- read.csv(file.path(plots_dir, "donor_metadata.csv"))

  pathways_and_signatures <- read.csv(file.path(plots_dir, "pathways_and_signatures.csv.gz"), row.names=1)

  pathways_and_signatures <- pathways_and_signatures %>% left_join(donor.meta.data, by="donor_id")

  # Change pathway names with . to _
  colnames(pathways_and_signatures) <- gsub("\\.", " ", colnames(pathways_and_signatures))

  # Filter out Young and non AD patients
  young_stages <- paste0(1:40, "-year-old stage")

  pathways_and_signatures <- pathways_and_signatures %>%
    filter(!(development_stage %in% young_stages)) %>%
    filter(AD_status == "Yes" | disease == "normal")

  h5createGroup(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project))
  subtype_column <- if (cell.type == "inhibitory") "subclass" else NULL
  h5_save_mean_scores(pathways_and_signatures, cell.type, project, group.by='donor_id', subtype_column=subtype_column)

  c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("SignatureUnique_")))
  h5_save_correlation(c, cell.type, paste0(project, "/signature_unique_correlation"))

  if (pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")) %>% ncol() > 0) {
    c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")))
    h5_save_correlation(c, cell.type, paste0(project, "/pathway_correlation"))
  }
}

# Trait association -----
project_dir <- file.path(get_result_dir(), project)

for (neuronal_class in c("inhibitory", "excitatory")) {
  donor.meta.data <- read.csv(file.path(project_dir, neuronal_class, "donor_metadata.csv"), check.names = FALSE)
  donor.meta.data <- donor.meta.data %>% left_join(st_github, by=c("donor_id"="DonorID"), suffix = c(".github", ".sc"))

  df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/mean_scores")) %>% filter(!is.na(donor_id)) %>% column_to_rownames("donor_id")

  if (neuronal_class == "inhibitory") {
    df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"),
      SignatureUnique_k1, SignatureUnique_k7, SignatureUnique_k8, SignatureUnique_k10)
  } else {
    df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"))
  }

  traits <- donor.meta.data %>% column_to_rownames("donor_id") %>% mutate(AD_status = if_else(AD_status == "Yes", 1, 0),
    Braak_4_5_6=if_else(Braak %in% c(4, 5, 6), 1, 0),
    HIGH_CERAD=if_else(CERAD %in% c(3, 4), 1, 0))
  traits <- traits[rownames(df),,drop=FALSE]
  controls <- donor.meta.data %>% column_to_rownames("donor_id") %>% mutate(age=as.numeric(gsub("\\+", "", Age)), sex=if_else(Sex == "M", 1, 0)) %>% dplyr::select(age, sex, PMI)
  controls <- controls[rownames(df),,drop=FALSE]

  withr::with_pdf(file.path(project_dir, neuronal_class, "TraitViolinPlotSignatures.pdf"), {
    plot_violin_by_group(df_signature_unique, traits$AD_status, ylab = "Signature score", xlab = "AD_status", title = "AD vs. Control") %>% print()
    plot_violin_by_group(df_signature_unique, traits$Braak, ylab = "Signature score", xlab = "Braak stage", title = "Braak stage vs. Control") %>% print()
    plot_violin_by_group(df_signature_unique, traits$CERAD, ylab = "Signature score", xlab = "CERAD", title = "CERAD vs. Control") %>% print()
  })

  withr::with_pdf(file.path(project_dir, neuronal_class, "TraitViolinPlotPathways.pdf"), {
    plot_violin_by_group(df %>% dplyr::select(starts_with("Pathway_")), traits$AD_status, ylab = "Pathway score", xlab = "AD_status", title = "AD vs. Control") %>% print()
    plot_violin_by_group(df %>% dplyr::select(starts_with("Pathway_")), traits$Braak, ylab = "Pathway score", xlab = "Braak stage", title = "Braak stage vs. Control") %>% print()
    plot_violin_by_group(df %>% dplyr::select(starts_with("Pathway_")), traits$CERAD, ylab = "Pathway score", xlab = "CERAD", title = "CERAD vs. Control") %>% print()
  })

  table_age <- table(traits$CERAD, traits$AD_status)
  barplot(
    t(table_age),
    beside = TRUE,
    col = c("skyblue", "salmon"),
    legend.text = c("Control", "AD"),
    args.legend = list(title = "AD status", x = "topright"),
    xlab = "Braak stage",
    ylab = "Count",
    main = "AD status by Braak stage"
  ) %>% print

  trait.association <- associate.traits(traits %>% dplyr::select(AD=AD_status, Braak, Braak_4_5_6, CERAD, HIGH_CERAD), df_signature_unique, controls)
  withr::with_pdf(file.path(project_dir, neuronal_class, "TraitAssociationSignature.pdf"), {
    plot.trait.associations.heatmap(trait.association, params = c("AD", "Braak", "Braak_4_5_6", "CERAD", "HIGH_CERAD"), row.by='covariate', show.only.significant = FALSE, column_title = glue::glue("Trait Association {neuronal_class}"))
  })

  trait.association_pathway <- associate.traits(traits %>% dplyr::select(AD=AD_status, Braak, Braak_4_5_6, CERAD, HIGH_CERAD), df %>% dplyr::select(starts_with("Pathway_")), controls)
  withr::with_pdf(file.path(project_dir, neuronal_class, "TraitAssociationPathway.pdf"), {
    plot.trait.associations.heatmap(trait.association_pathway, params = c("AD", "Braak", "Braak_4_5_6", "CERAD", "HIGH_CERAD"), row.by='covariate', show.only.significant = FALSE, column_title = glue::glue("Trait Association {neuronal_class}"))
  })

  h5write(trait.association, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_signature"))
  h5write(trait.association_pathway, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_pathway"))
}

# Community analysis -----

compute_program_environment_correlation(project, identity_column="donor_id")
