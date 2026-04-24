#' @title Proteomics Analysis Functions
#' @description Functions for analyzing proteomic data and associating with traits.
#' @note Requires: org.Hs.eg.db, SummarizedExperiment

library('org.Hs.eg.db')
library("SummarizedExperiment")

PROTEOMIC_PATH <- file.path(get_result_dir(), "../../data/proteomics.rds")
proteom <- readRDS(PROTEOMIC_PATH)


#' Create UniProt to gene symbol mappings
#'
#' @description Extracts gene symbol and UniProt ID mappings from a 
#'              SummarizedExperiment proteomics object.
#'
#' @param proteom SummarizedExperiment object with rowData containing Symbol and UniProt columns
#'
#' @return Named list with two elements:
#'         - genes: Named vector mapping UniProt IDs to gene symbols
#'         - uniprot: Named vector of UniProt IDs (self-mapping for lookup)
uniprot.mapping <- function(proteom) {
  genes_mapping <- rowData(proteom)$Symbol
  names(genes_mapping) <- make.names(rowData(proteom)$UniProt)
  UniProt_mapping <- rowData(proteom)$UniProt
  names(UniProt_mapping) <- make.names(UniProt_mapping)

  return(list(genes=genes_mapping, uniprot=UniProt_mapping))
}
# uniprot.mapping <- memoize_first(uniprot.mapping)

.UniProt_mapping <- rowData(proteom)$UniProt
names(.UniProt_mapping) <- make.names(.UniProt_mapping)
.genes_mapping <- rowData(proteom)$Symbol
names(.genes_mapping) <- make.names(rowData(proteom)$UniProt)

#' Associate proteomic measurements with traits
#'
#' @description Performs trait association analysis for specified genes/proteins
#'              in proteomics data, with optional signature grouping.
#'
#' @param proteom SummarizedExperiment object with proteomics data
#' @param traits Character vector of trait column names from colData to test
#' @param genes Character vector of gene symbols to include as covariates
#' @param control_genes Character vector of gene symbols to use as controls (default: empty)
#' @param add_signature If TRUE, add aggregated signature scores for gene groups (default: FALSE)
#' @param signature_groups Named list of gene vectors for signature computation
#' @param ... Additional arguments passed to associate.traits()
#'
#' @return Data frame with trait association results including:
#'         - Standard association columns (trait, covariate, beta, pval, etc.)
#'         - gene: Gene symbol
#'         - UniProt: UniProt ID
#'         - gene_uniprot: Combined "Gene (UniProt)" label
#'
#' @details When add_signature=TRUE, each signature group is computed as the
#'          row mean of scaled expression values for genes in that group.
#'          Controls include age_death, msex, pmi, and any specified control_genes.
associate.proteomic.traits <- function (
  proteom,
  traits,
  genes,
  control_genes = c(),
  add_signature = FALSE,
  signature_groups = c(),
  ...
) {
  
  
  covariates <- assay(proteom)[rowData(proteom)$Symbol %in% genes,,drop=F] %>% t %>% as.data.frame()
  
  if (add_signature) {
    for (group in names(signature_groups)) {
      covariates[[group]] <- assay(proteom)[rowData(proteom)$Symbol %in% signature_groups[[group]],,drop=F] %>%
        t %>%
        as.data.frame() %>% 
        scale() %>%
        replace_na(0) %>%
        rowMeans()
    }
  }
  
  traits <- colData(proteom)[, traits] %>% as.data.frame()
  controls <-  data.frame(colData(proteom)[,c("age_death","msex","pmi")],
                          assay(proteom)[rowData(proteom)$Symbol %in% control_genes,,drop=F] %>% t
                          )
  
  
  trait.association <- associate.traits(
    traits = traits,
    covariates = covariates,
    controls = controls,
    ...
  )
  
  trait.association <- trait.association %>% mutate(gene=.genes_mapping[make.names(covariate)],
                                                    UniProt=.UniProt_mapping[make.names(covariate)],
                                                    gene_uniprot=glue::glue("{gene} ({UniProt})"))
  trait.association
}
