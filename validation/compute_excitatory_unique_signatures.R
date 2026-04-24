#' @title Compute unique gene signatures for excitatory and inhibitory and store in H5
#' @description Runs the same logic as in cuimc2.R (and other validation scripts)
#'   for each cell type: get_genes_by_topic() and stores the unique signature
#'   as a dataset "program_signature" (columns: topic, gene). Does not modify
#'   existing validation scripts (they continue to call get_genes_by_topic as before).
#' @usage Run from project root: Rscript validation/compute_excitatory_unique_signatures.R

library(rhdf5)
library(dplyr)

source("R/utils.R")
source("R/data_500_utils.R")
source("validation/validation.functions.R")

CELL_TYPES <- c("excitatory", "inhibitory")
H5_NAME <- "program_signature"

for (CELL_TYPE in CELL_TYPES) {
  message("Processing ", CELL_TYPE, " ...")
  H5_GROUP <- CELL_TYPE_TO_GROUP[[CELL_TYPE]]

  message("  Computing unique_genes_by_topic via get_genes_by_topic() ...")
  genes_mapping <- get_genes_by_topic(CELL_TYPE)
  unique_genes_by_topic <- genes_mapping$unique_genes_by_topic

  # Long-format table: topic, gene only

  df <- lapply(names(unique_genes_by_topic), function(topic) {
    data.frame(topic = topic, gene = unique_genes_by_topic[[topic]], stringsAsFactors = FALSE)
  })
  df <- bind_rows(df)

  if (nrow(df) == 0L) {
    message("  No unique genes; skipping write.")
    next
  }

  # Ensure group exists
  h5createGroup(TOPICS_DATA_H5_PATH, H5_GROUP)

  # Remove existing dataset if present so h5write can replace
  tryCatch(
    h5delete(TOPICS_DATA_H5_PATH, paste0(H5_GROUP, "/", H5_NAME)),
    error = function(e) NULL
  )

  h5write(df, TOPICS_DATA_H5_PATH, paste0(H5_GROUP, "/", H5_NAME))
  message("  Wrote ", nrow(df), " rows to ", paste0(H5_GROUP, "/", H5_NAME), " (topic, gene)")
}
