library(dplyr)
library(tibble)

source("R/utils.R")
source("R/figures.funs.R")

read.bulk <- function (path="../data/geneTpmResidualsAgeGenderUnadj") {
  bulk <- read.table(path, as.is=T, header=T, row.names = 1) %>%
  `colnames<-`(as.character(as.numeric(gsub("X","",colnames(.)))))

  bulk
}

bulk <- read.bulk()


symbol.ensembl.mapping <- function (organism = c("hsa", "mmu")) {
  organism <-  match.arg(organism)

  if (organism == "hsa") {
    suppressPackageStartupMessages(require(org.Hs.eg.db))
    df <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"), columns=c("SYMBOL","ENSEMBL")))
  } else { if (organism == "mmu") {
    suppressPackageStartupMessages(require(org.Mm.eg.db))
    df <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.eg.db, "ENTREZID"), columns=c("SYMBOL","ENSEMBL")))
  }}
  
  df
}


shared.donors <- intersect(colnames(bulk), rownames(data))
# inhibitory.pseudobulks <- inhibitory.pseudobulks[shared.donors,]
# cux2m.pseudobulks <- cux2m.pseudobulks[shared.donors,]
bulk <- bulk[, shared.donors] %>% t


shared.genes <- function(bulk, ct.psuedobulk) {
  mapping <- symbol.ensembl.mapping() %>%
    mutate(symbol_ensembl=paste(SYMBOL, ENSEMBL),
           has_psuedo_bulk=SYMBOL %in% colnames(ct.psuedobulk),
           has_bulk=ENSEMBL %in% colnames(bulk))
  
  mapping %>% filter(has_psuedo_bulk & has_bulk) %>% pull(SYMBOL)
}




correlate.genes.psuedobulk.bulk <- function(bulk, pseudobulks, genes, cell.types=NULL) {
  if (is.null(cell.types)) {
    cell.types <- c("default")
    pseudobulks <- list(default=pseudobulks)
  }
  bulk.projids <- bulk %>% rownames()
  
  cell.types.cor <- sapply(cell.types, function (cell.type) {

    ct.pseudobulk <- pseudobulks[[cell.type]]
    
    sc.projids <- rownames(ct.pseudobulk)
    projids <- intersect(sc.projids, bulk.projids)
    
    ct.pseudobulk <- ct.pseudobulk[projids, ]
    bulk <- bulk[projids,]
    
    mapping <- symbol.ensembl.mapping() %>%
      filter(SYMBOL %in% genes) %>%
      mutate(symbol_ensembl=paste(SYMBOL, ENSEMBL),
             has_bulk=ENSEMBL %in% colnames(bulk)) %>% 
      filter(has_bulk)
    
    f <- function(symbol, ensembl) {
      p <- cor.test(ct.pseudobulk[, symbol], bulk[,ensembl], use="pairwise.complete.obs", method = "spearman", exact=FALSE)
      c(corr=p$estimate, p.value=p$p.value)
    }
    
    cor <- mapping %>% bind_cols(Vectorize(f, SIMPLIFY = T, USE.NAMES = F)(.$SYMBOL, .$ENSEMBL) %>% t) %>% 
      mutate(Symbol=SYMBOL) %>% 
      mutate(adj.pval=p.adjust(p.value, method = "BH"),
             sig=cut(adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")))
    
    cor
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  if (is.null(cell.types)) {
    cell.types.cor <- cell.types.cor$default
  }
  
  cell.types.cor
}

