SHINY_PATH <- file.path(get_result_dir(), "../../Shared/shiny_geneExpression_500")
SHINY_RESOURCES.PATH <- file.path(SHINY_PATH, "resources")

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

load.fits <- function(topics=FALSE){
  dir <- ifelse(topics, dynamics_topics_dir, dynamics_dir)
  sapply(gsub(".rds", "", list.files(dir)), function(ct){
    readRDS(file.path(dir, paste0(ct, ".rds")))
  }, USE.NAMES = T, simplify = F)
}
