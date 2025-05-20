library('org.Hs.eg.db')
library("SummarizedExperiment")

PROTEOMIC_PATH <- file.path(get_result_dir(), "../../data/proteomics.rds")
proteom <- readRDS(PROTEOMIC_PATH)



#######################
# Trait Association
#######################

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



##################################
# Correlation with pseudo-bulk
##################################

corr.proteom.pseudobulk <- function(proteom, pseudobulks, genes) {
  sc.projids <- pseudobulks %>% rownames()
  protem.projids <- assay(proteom) %>% colnames
  projids <- intersect(sc.projids, protem.projids)

  proteom.df <- assay(proteom)[rowData(proteom)$Symbol %in% genes, projids, drop=F] %>% t

  ps.bulk <- pseudobulks[projids, genes, drop=FALSE]

  
  f <- function(symbol, UniProt) {
    p <- cor.test(ps.bulk[, symbol], proteom.df[,UniProt], use="pairwise.complete.obs", method = "spearman", exact=FALSE)
    c(corr=p$estimate, p.value=p$p.value)
  }
  
  
  rowData(proteom) %>%
    as.data.frame() %>% 
    filter(Symbol %in% genes) %>% 
    mutate(name=glue::glue("{Symbol} ({UniProt})")) %>% 
    arrange(Symbol) %>% 
    bind_cols(Vectorize(f, SIMPLIFY = T, USE.NAMES = F)(.$Symbol, .$UniProt) %>% t) %>% 
    mutate(adj.pval=p.adjust(p.value, method = "BH"),
           sig=cut(adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")))
}


corr.proteom.bulk <- function(proteom, bulk, genes) {
  bulk <- bulk %>% t
  
  bulk.projids <- bulk %>% rownames()
  protem.projids <- assay(proteom) %>% colnames
  projids <- intersect(bulk.projids, protem.projids)
  
  proteom.df <- assay(proteom)[rowData(proteom)$Symbol %in% genes, projids, drop=F] %>% t
  
  bulk <- bulk[projids, , drop=FALSE]
  
  bulk_mapping <- symbol.ensembl.mapping() %>%
    filter(SYMBOL %in% genes) %>%
    mutate(symbol_ensembl=paste(SYMBOL, ENSEMBL),
           has_bulk=ENSEMBL %in% colnames(bulk)) %>% 
    dplyr::rename(Symbol=SYMBOL) %>% 
    filter(has_bulk)
  
  f <- function(ENSEMBL, UniProt) {
    p <- cor.test(bulk[, ENSEMBL], proteom.df[,UniProt], use="pairwise.complete.obs", method = "spearman", exact=FALSE)
    c(corr=p$estimate, p.value=p$p.value)
  }
  
  rowData(proteom) %>%
    as.data.frame() %>% 
    filter(Symbol %in% genes) %>% 
    mutate(name=glue::glue("{Symbol} ({UniProt})")) %>% 
    dplyr::left_join(bulk_mapping) %>% 
    filter(has_bulk == TRUE) %>% 
    arrange(Symbol) %>% 
    bind_cols(Vectorize(f, SIMPLIFY = T, USE.NAMES = F)(.$ENSEMBL, .$UniProt) %>% t) %>% 
    mutate(adj.pval=p.adjust(p.value, method = "BH"),
           sig=cut(adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")))
}



plot.bulk.corr <- function(df, genes=df$Symbol %>% unique, colors=cyan2orange,
                           column_title="Proteomic - Psuedobulk correlation", ...) {
  df <- df %>% filter(Symbol %in% genes)
  rownames(df) <- df$name
  
  v <- max(abs(df$corr.rho), na.rm = T)
  col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])
  colors <- circlize::colorRamp2(col.vals, colors(length(col.vals)))
  
  hm <- Heatmap(df['corr.rho'] %>% as.matrix(), 
                cell_fun = function(j, i, x, y, w, h, fill) grid.text(df$sig[i], x,y), 
                column_title = column_title,
                show_column_names = F,
                col = colors,
                show_row_dend = F, ...)
  
  legend <- Legend(title = "adj.pval", pch = c("****","***","**","*"), type = "points", labels = c("<0.0001","<0.001","<0.01", "<0.05"))
  draw(hm, annotation_legend_list = list(legend), merge_legend=T)
}

########################
# Proteomics Dynamics
########################

fit.proteomics.dynamics <- function (proteom, data, genes) {
  sc.projids <-rownames(data)
  protem.projids <- assay(proteom) %>% colnames
  projids <- intersect(sc.projids, protem.projids)
  
  features <- assay(proteom)[rowData(proteom)$Symbol %in% genes, projids] %>% 
    as.data.frame() %>% 
    rownames_to_column("UniProt") %>% 
    mutate(gene=.genes_mapping[make.names(UniProt)],
           gene_uniprot=glue::glue("{gene} ({UniProt})")) %>% 
    dplyr::select(-gene, -UniProt) %>%
    column_to_rownames("gene_uniprot") %>% 
    t
  
  sc.only.projids <- setdiff(sc.projids, projids)
  empty_features <- matrix(ncol = ncol(features), nrow = length(sc.only.projids), dimnames = list(sc.only.projids, colnames(features))) %>% as.data.frame()
  features <- rbind(features, empty_features)
  features <- features[sc.projids, ]
  
    
  dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                           features = features,
                           trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                           trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                           evaluate.fit = T,
                           bootstrap = F)

  dynamics
}



plot.pathways.proteomic <- function(pathways_path, proteom, psuedobulk, output_dir, pathways_ids=NULL, show.only.significant=FALSE) {
  dir.create(output_dir, showWarnings = FALSE)
  
  pathways <- read.csv(pathways_path, row.names = 1)
  pathways <- pathways %>% filter(plot) %>% filter(is.null(pathways_ids) | ID %in% pathways_ids)
  for (topic in (pathways$topic %>% unique)) {
    print(topic)
    withr::with_pdf(file.path(output_dir, glue::glue("{topic}_proteom.pdf")), {
      topic_pathways <- pathways %>% filter(.data[['topic']] == .env[['topic']])
      
      lapply(1:nrow(topic_pathways), function (i){
        genes <- strsplit(topic_pathways[i, 'geneName'], "/") %>% unlist()

        if (0 == sum(rowData(proteom)$Symbol %in% genes)) {
          print(grid.text(glue::glue("No protemics for pathway {topic_pathways[i, 'ID']} - {topic_pathways[i, 'Description']}")))
          return(0)
        }
        
        print(glue::glue("Genes of pathways - {topic_pathways[i, 'ID']} - {topic_pathways[i, 'Description']}"))

        psuedobulk.bulk.corr.df <- correlate.genes.psuedobulk.bulk(bulk, psuedobulk, genes)
        
        traits_col <- c("amyloid_sqrt", "tangles_sqrt", "cogng_random_slope")
        prot.traits.df <- associate.proteomic.traits(proteom, genes=genes, traits = traits_col)
        
        corr.df <- corr.proteom.pseudobulk(proteom, psuedobulk, genes)
        
        plot.proteom(prot.traits.df, corr.df, psuedobulk.bulk.corr.df,
                     traits_col, row.by='gene_uniprot',
                     show.only.significant=show.only.significant,
                     column_title = glue::glue("Trait association - {topic_pathways[i, 'Description']}"),
                     column_title_gp = gpar(fontsize = 11))

      })
    })
  }
}


plot.proteom <- function(df, proteom.corr.df, pseudobulk.bulk.corr,
                         params=c("sqrt.amyloid_mf","sqrt.tangles_mf","cogng_demog_slope"),
                         row.by="state", 
                         col.by="trait", value.by="tstat", pval.by="adj.pval",
                         cols=window2fox,
                         column_labels=NULL,
                         show.only.significant=T, plot=T, ...) {
  if("pandas.core.frame.DataFrame" %in% class(df))
    df <- py_to_r(df)
  
  df <- df %>% filter(trait %in% params)
  
  if(show.only.significant)
    df <- df %>% group_by_at(c(row.by)) %>%
    filter(if_any(pval.by, ~sum(. < .05) > 0)) %>%
    ungroup() %>%
    data.frame()
  
  df$sig <- cut(df[[pval.by]], c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", ""))
  
  if (nrow(df) == 0) {
    print(grid::grid.text("No trait associated"))
    return()
  }
  
  vals <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = all_of(value.by), values_fill = NA_real_) %>%
    tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params)) %>% as.matrix()
  sig <- tidyr::pivot_wider(df, id_cols = all_of(row.by), names_from = all_of(col.by), values_from = "sig", values_fill = "") %>%
    tibble::column_to_rownames(row.by) %>% dplyr::select(all_of(params)) %>% as.matrix()
  
  v <- max(abs(vals), na.rm = T)
  col.vals <- c(seq(-v, 0, length.out=11), seq(0, v, length.out=11)[-1])
  
  proteom.corr.df <- proteom.corr.df %>% filter(name %in% rownames(vals)) %>%  arrange(match(name, rownames(vals)))
  
  pseudobulk.annotations <- lapply(names(pseudobulk.bulk.corr), function(cell.type) {
    psuedobulk.bulk.corr.df <- df %>% dplyr::distinct(gene, UniProt) %>% left_join(pseudobulk.bulk.corr[[cell.type]], by=c("gene"="SYMBOL")) %>% mutate(sig=unfactor(sig))
    if (nrow(psuedobulk.bulk.corr.df) != nrow(vals)) {
      warn("Bulk correaltion length is different from proteomics traits")
      psuedobulk.bulk.corr.df <- psuedobulk.bulk.corr.df %>%  group_by(gene, UniProt, ENTREZID) %>%
        summarise(sig=if_else(max(corr.rho) > 0 & min(corr.rho) < 0, '#', max(sig)), corr.rho=max(abs(corr.rho))) %>% 
        ungroup %>% 
        arrange(match(gene, df$gene))
    }
    
    annotation <- anno_simple(psuedobulk.bulk.corr.df$corr.rho, col=cyan2orange.correlation, pch=psuedobulk.bulk.corr.df$sig, na_col = 'black',
                              which = 'row')
    
    return(annotation)
  })
  names(pseudobulk.annotations) <- names(pseudobulk.bulk.corr)
  
  # args <- c(protem.correlation=anno_simple(proteom.corr.df$corr.rho, col=cyan2orange.correlation, pch=unfactor(proteom.corr.df$sig)), 
  #           pseudobulk.correlation=pseudobulk.annotations)
  # 
  # TODO: change to row_ha = do.call(rowAnnotation, args=args)
  
  if (!is.null(column_labels) & length(column_labels) != ncol(vals)) {
    column_labels <- column_labels[colnames(vals)]
  } else if (is.null(column_labels)) {
    column_labels <- colnames(vals)
  }
  
  row_ha <- rowAnnotation(protem.correlation=anno_simple(proteom.corr.df$corr.rho, col=cyan2orange.correlation, pch=unfactor(proteom.corr.df$sig))
                         ,
                         pb.corr.inhibitory=pseudobulk.annotations$inhibitory,
                         pb.corr.exc.lower=pseudobulk.annotations$`cux2-`,
                         pb.corr.exc.upper=pseudobulk.annotations$`cux2+`
                         )
  
  hm <- Heatmap(vals,
                left_annotation=row_ha,
                name = value.by,
                col = circlize::colorRamp2(col.vals, cols(length(col.vals))),
                cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
                cluster_columns = F, clustering_distance_rows = "euclidean",
                show_row_dend = F,
                column_labels = column_labels,
                ...)
  legend <- Legend(title = "adj.pval", pch = c("****","***","**","*"), type = "points", labels = c("<0.0001","<0.001","<0.01", "<0.05"))
  correlation.legend <- Legend(title = "correlation", col_fun = cyan2orange.correlation)
  if(!plot)
    return(list(hm=hm, legend=legend))
  draw(hm, annotation_legend_list = list(legend, correlation.legend), merge_legend=T)
  
}
