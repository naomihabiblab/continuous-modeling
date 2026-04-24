#' @title Pseudobulk Data Loading Functions
#' @description Functions for loading pseudobulk expression data aggregated by donor.

SHINY_PATH <- file.path(get_result_dir(), "../../../Shared/shiny_geneExpression_500")
SHINY_RESOURCES.PATH <- file.path(SHINY_PATH, "resources")

#' Load pseudobulk expression matrices by donor (projid)
#'
#' @description Loads pre-computed pseudobulk expression data for multiple cell types,
#'              returning matrices with donors as rows and genes as columns.
#'
#' @param cts Character vector of cell types to load (default: main cell types)
#' @param psuedobulk_path Path to directory containing pseudobulk RDS files
#'
#' @return Named list of data frames, one per cell type:
#'         - Rows: donor IDs (projid)
#'         - Columns: genes (features)
#'         - Values: average expression
#'
#' @details Files are expected to be named "dotplot_{celltype}_projid.rds" and
#'          contain columns: features.plot, id, avg.exp
load.pseudobulks.of.projid <- function(cts=c("astrocytes", "microglia", "opcs", "oligodendrocytes", "inhibitory", "cux2-", "cux2+"),
                                       psuedobulk_path=file.path(SHINY_RESOURCES.PATH, "pseudobulks")){
  sapply(cts, function(ct){
    file_name <- file.path(psuedobulk_path, paste0("dotplot_", ct, "_projid.rds"))

    if(file.exists(file_name)){
      df <- (paste0(file_name) %>% readRDS())

      return(df  %>% dplyr::select(c("features.plot", "id", "avg.exp")) %>% tidyr::spread(., features.plot, avg.exp) %>% tibble::column_to_rownames("id"))
    }
  }, simplify = F, USE.NAMES = T)
}
