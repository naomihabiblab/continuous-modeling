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

SEA_AD_SUPP_TABLE_PATH <- file.path(get_shared_dir(), "SeaAD/paper_data/SupplementaryTables.xlsx")

project <- "seaad"
cell_types <- c("inhibitory", "excitatory", "astrocytes", "microglia", "oligo", "opcs")

for (cell.type in cell_types) {
  project_dir <- file.path(get_result_dir(), project)
  plots_dir <- file.path(project_dir, cell.type)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  donor.meta.data <- read.csv(file.path(plots_dir, "donor_metadata.csv"), check.names = FALSE)

  pathways_and_signatures <- read.csv(file.path(plots_dir, "pathways_and_signatures.csv.gz"), row.names=1, check.names = FALSE)

  pathways_and_signatures <- pathways_and_signatures %>% left_join(donor.meta.data, by="Donor ID")

  h5delete(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project))

  h5createGroup(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell.type], "/", project))
  subtype_column <- if (cell.type == "inhibitory") "Subclass" else NULL
  h5_save_mean_scores(pathways_and_signatures, cell.type, project, group.by='Donor ID', subtype_column=subtype_column)

  c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("SignatureUnique_")))
  h5_save_correlation(c, cell.type, paste0(project, "/signature_unique_correlation"))

  if (pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")) %>% ncol() > 0) {
    c <- correlate_dataframes(pathways_and_signatures %>% dplyr::select(starts_with("Pathway_")))
    h5_save_correlation(c, cell.type, paste0(project, "/pathway_correlation"))
  }
}

# Trait association -----
project_dir <- file.path(get_result_dir(), project)
pathologies_df <- readxl::read_excel(SEA_AD_SUPP_TABLE_PATH, sheet = "ST2")

for (neuronal_class in c("inhibitory", "excitatory")) {
  donor.meta.data <- read.csv(file.path(project_dir, neuronal_class, "donor_metadata.csv"), check.names = FALSE)
    donor.meta.data <- donor.meta.data %>% left_join(pathologies_df, by="Donor ID")

    df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/mean_scores")) %>% filter(!is.na(`Donor ID`)) %>% column_to_rownames("Donor ID")

    if (neuronal_class == "inhibitory") {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"), 
        SignatureUnique_k1, SignatureUnique_k7, SignatureUnique_k8, SignatureUnique_k10)
    } else {
        df_signature_unique <- df %>% dplyr::select(starts_with("SignatureUnique_increasing"), starts_with("SignatureUnique_decreasing"))
    }


    traits <- donor.meta.data %>% 
        column_to_rownames("Donor ID") %>%
        mutate(Dementia = if_else(`Cognitive Status` == "Dementia", 1, 0), 
        `Severely Affected Donor` = if_else(`Severely Affected Donor` == "Y", 1, 0)) 
    traits <- traits[rownames(df),,drop=FALSE]

    # Barplot for traits$Dementia over all Overall AD neuropathological Change categories
    dementia_oadnc_table <- table(traits$`Overall AD neuropathological Change`, traits$Dementia)
    withr::with_pdf(file.path(project_dir, neuronal_class, "DementiaOverallADNeuropathologicalChangeBarplot.pdf"), {
        barplot(
        t(dementia_oadnc_table),
        beside = TRUE,
        col = c("skyblue", "salmon"),
        legend.text = c("Control", "Dementia"),
        args.legend = list(title = "Dementia", x = "topright"),
        xlab = "Overall AD neuropathological Change",
            ylab = "Count",
            main = "Dementia Status by Overall AD neuropathological Change"
        ) %>% print
    })

    # For each column in df_signature_unique, perform a t-test comparing groups defined by AD_status
    controls <- data.frame(row.names=rownames(df), sex=traits$Sex, pmi=traits$PMI, age=traits$`Age at Death`)
    trait.association <- associate.traits(traits %>% dplyr::select(Dementia=Dementia, cps="Continuous Pseudo-progression Score", amyloid="percent 6e10 positive area_Grey matter", tangles="percent AT8 positive area_Grey matter"), df_signature_unique, controls)
    trait.association_pathway <- associate.traits(traits %>% dplyr::select(Dementia=Dementia, cps="Continuous Pseudo-progression Score", amyloid="percent 6e10 positive area_Grey matter", tangles="percent AT8 positive area_Grey matter"), df %>% dplyr::select(starts_with("Pathway_")), controls)


    trait.association <- associate.traits(traits %>% dplyr::select(Dementia=Dementia, cps="Continuous Pseudo-progression Score", amyloid="percent 6e10 positive area_Grey matter", tangles="percent AT8 positive area_Grey matter"), df_signature_unique, controls)
    trait.association_pathway <- associate.traits(traits %>% dplyr::select(Dementia=Dementia, cps="Continuous Pseudo-progression Score", amyloid="percent 6e10 positive area_Grey matter", tangles="percent AT8 positive area_Grey matter"), df %>% dplyr::select(starts_with("Pathway_")), controls)

    withr::with_pdf(file.path(project_dir, neuronal_class, "TraitAssociationSignature.pdf"), {
      plot.trait.associations.heatmap(trait.association, params = c("Dementia", "cps", "amyloid", "tangles"), row.by='covariate', show.only.significant = FALSE)
    })

    withr::with_pdf(file.path(project_dir, neuronal_class, "TraitAssociationPathway.pdf"), {
      plot.trait.associations.heatmap(trait.association_pathway, params = c("Dementia", "cps", "amyloid", "tangles"), row.by='covariate', show.only.significant = FALSE)
    })

  h5write(trait.association, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_signature"))
  h5write(trait.association_pathway, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[neuronal_class], "/", project, "/trait_association_pathway"))
}

# Community analysis -----

compute_program_environment_correlation(project, identity_column="Donor ID")
