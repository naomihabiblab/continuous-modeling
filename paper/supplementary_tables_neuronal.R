#' @title Neuronal-program sheets of the Supplementary Tables
#' @description Builds the data frames of the "program unique signature" sheet
#'              (Supplementary Table 2) and of the per-dataset "Validation" sheets
#'              (Supplementary Table 3), with the program names used in the paper.
#'              Used by paper/create.supplementary.tables.R.
#'
#' Requires R/utils.R, R/data_500_utils.R and R/topic_space_functions.R to be sourced.

library(dplyr)

NEURONAL_PROGRAM_GROUP_LABELS <- c(
  decreasing_early = "decreasing early",
  decreasing_late  = "decreasing late",
  increasing_early = "increasing early (stress/apoptosis)",
  increasing_late  = "increasing late (cell-cycle reentry)"
)

#' Paper names of the programs behind each neuronal program group.
#'
#' @return Data frame with program_group (internal key), program_group_label,
#'         inhibitory (e.g. "INH.p1") and excitatory (e.g. "EXC L2-3.p9 + EXC L4-6.p4").
neuronal_program_groups <- function() {
  do.call(rbind, lapply(names(neuronal_topics_to_annotation), function(group) {
    topics <- names(neuronal_topics_to_annotation[[group]])
    paper_names <- topic_names_mapping(topics, use_static_mapping = FALSE)
    data.frame(
      program_group = group,
      program_group_label = NEURONAL_PROGRAM_GROUP_LABELS[[group]],
      inhibitory = paper_names[grepl("^Inh", topics)],
      excitatory = paste(paper_names[grepl("^Exc", topics)], collapse = " + "),
      stringsAsFactors = FALSE
    )
  }))
}

#' Rows of the "program unique signature" sheet (Supplementary Table 2).
#'
#' @return Data frame: cell.type, program, program_group, gene.
program_signature_sheet <- function(h5_path = TOPICS_DATA_H5_PATH) {
  groups <- neuronal_program_groups()
  lapply(c("inhibitory", "excitatory"), function(ct) {
    path <- paste0(CELL_TYPE_TO_GROUP[[ct]], "/program_signature")
    df <- tryCatch(h5read(h5_path, path), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      dplyr::left_join(groups, by = c("topic" = "program_group")) %>%
      dplyr::mutate(cell.type = ct,
                    program = if (ct == "inhibitory") inhibitory else excitatory,
                    program_group = program_group_label) %>%
      dplyr::select(cell.type, program, program_group, gene)
  }) %>% dplyr::bind_rows()
}

#' Readable pathway name from a validation covariate ("Pathway_DNA_Repair" -> "DNA Repair").
pathway_label <- function(covariate) {
  x <- gsub("_", " ", sub("^Pathway_", "", covariate))
  x <- gsub("HSF1 mediated", "HSF1-mediated", x)
  sub("^antinsense$", "antisense", x)
}

#' Data frames of the per-dataset "Validation <dataset>" sheets (Supplementary Table 3).
#'
#' Adds, after `covariate`: `program` and `program_group` (signature rows; for excitatory
#' neurons the signature merges the matching EXC L2-3 and EXC L4-6 programs) and
#' `pathway` (pathway rows).
#'
#' @return Named list (sheet name -> data frame), in the order SeaAD, PsychAD, CUIMC2, MIT.
validation_association_sheets <- function(h5_path = TOPICS_DATA_H5_PATH) {
  validation_projects <- list(seaad = "SeaAD", psychAD = "PsychAD", cuimc2 = "CUIMC2", mit = "MIT")
  groups <- neuronal_program_groups()
  inh_k_to_group <- setNames(groups$program_group, sub(".*\\.p", "k", groups$inhibitory))

  sheets <- list()
  for (proj in names(validation_projects)) {
    display_name <- validation_projects[[proj]]
    parts <- list()
    for (neuronal_class in c("inhibitory", "excitatory")) {
      group_path <- paste0(CELL_TYPE_TO_GROUP[[neuronal_class]], "/", proj, "/")
      for (assoc_type in c("trait_association_signature", "trait_association_pathway")) {
        df_val <- tryCatch(h5read(h5_path, paste0(group_path, assoc_type)), error = function(e) NULL)
        if (is.null(df_val) || nrow(df_val) == 0) next
        df_val <- as.data.frame(df_val) %>%
          dplyr::mutate(dplyr::across(dplyr::any_of(c("beta", "se", "tstat", "pval", "r.sq", "n", "adj.pval")), as.numeric))
        type_label <- if (assoc_type == "trait_association_signature") "signature" else "pathway"
        df_val$validation_dataset <- display_name
        df_val$neuronal_class <- neuronal_class
        df_val$type <- type_label

        key <- sub("^SignatureUnique_", "", df_val$covariate)
        group <- if (neuronal_class == "inhibitory") unname(inh_k_to_group[key]) else key
        row <- match(group, groups$program_group)
        is_sig <- type_label == "signature"
        df_val$program <- if (is_sig) groups[[neuronal_class]][row] else NA_character_
        df_val$program_group <- if (is_sig) groups$program_group_label[row] else NA_character_
        df_val$pathway <- if (is_sig) NA_character_ else pathway_label(df_val$covariate)
        if (is_sig && anyNA(df_val$program)) {
          stop(glue::glue("Unmapped signature covariate in {display_name} {neuronal_class}: ",
                          "{paste(unique(df_val$covariate[is.na(df_val$program)]), collapse = ', ')}"))
        }
        parts[[paste(neuronal_class, type_label, sep = "_")]] <- df_val
      }
    }
    if (length(parts) > 0) {
      df <- dplyr::bind_rows(parts)
      df <- df %>% dplyr::relocate(program, program_group, pathway, .after = covariate)
      sheets[[paste0("Validation ", display_name)]] <- df
    }
  }
  sheets
}
