####################################################################################################################
##                                            #  Perform Pathways Analysis   #                                    ##
####################################################################################################################
library(igraph)
library(ggraph)
library(ggforce)

source("utils.R")

GeneIdMapping <- function(organism = c("hsa", "mmu")) {
  organism = match.arg(organism)
  
  if (organism == "hsa") {
    suppressPackageStartupMessages(require(org.Hs.eg.db))
    df <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"), columns=c("SYMBOL","ENTREZID")))
  } else { if (organism == "mmu") {
    suppressPackageStartupMessages(require(org.Mm.eg.db))
    df <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.eg.db, "ENTREZID"), columns=c("SYMBOL","ENTREZID")))
  }}
  return(list(ids=setNames(df$ENTREZID, df$SYMBOL), names=setNames(df$SYMBOL, df$ENTREZID)))
}


FindMarkersWrapper <- function(object, idents="ident", workers=NULL, ...) {
  Idents(object) <- FetchData(object, idents)
  
  if(!is.null(workers)) {
    library(future)
    plan(multicore, workers = workers)
  }
  de <- FindAllMarkers(object, ...)
  if (nrow(de) == 0)
    return(de)
  return(de %>% dplyr::mutate(id=GeneIdMapping()$ids[gene]))
}


FindMarkersSplitBy <- function(object, split.by, idents="ident", ...) {
  Idents(object) <- FetchData(object, idents)
  return(do.call(rbind, lapply(unique(Idents(object)), function(s) {
    print(paste("DE for", s))
    
    if (idents == "ident") 
      t <- object[,object@active.ident == s]
    else 
      t <- object[,object[[idents]] == as.character(s)]
    
    de <- FindMarkersWrapper(t, idents=split.by, ...)
    if (nrow(de) == 0) return(de)
    return(de %>% dplyr::mutate("{split.by}":=cluster, cluster=s)) 
  })))
}


EnrichmentAnalysis <- function(de, formula="id~cluster", fun=c("enrichKEGG","enrichPathway","enrichGO:BP","enrichGO:MF","enrichGO:CC"), species="Homo sapiens", universe=NULL, ...) {
  if(nrow(de) == 0) return(list())
  
  suppressPackageStartupMessages(require(clusterProfiler))
  f <- as.formula(formula)
  
  pathways <- list(`enrichKEGG` = function() compareCluster(f, data = de, fun="enrichKEGG", organism="hsa", universe=universe, ...),
                   `enrichPathway` = function() {
                     suppressPackageStartupMessages(require(ReactomePA))
                     compareCluster(f, data = de, fun="enrichPathway", universe=universe, ...)
                   },
                   `enrichGO:BP` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", keyType = "ENTREZID", OrgDb="org.Hs.eg.db", ont="BP", universe=universe, ...)),
                   `enrichGO:MF` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", keyType = "ENTREZID", OrgDb="org.Hs.eg.db", ont="MF", universe=universe, ...)),
                   `enrichGO:CC` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", keyType = "ENTREZID", OrgDb="org.Hs.eg.db", ont="CC", universe=universe, ...)))
  
  ea <- list()
  for(p in intersect(fun, names(pathways))) {
    tryCatch(ea[[p]] <- pathways[[p]](), error=function(e){print(e)})
  }

  
  for(category in fun[grepl("msigdbr", fun)]) {
    suppressPackageStartupMessages(require(msigdbr))
    terms <- msigdbr(species = species, category = gsub("^.*:", "", category)) %>% dplyr::select(gs_name, entrez_gene)
    
    idents <- unique(de[,gsub(".*~", "", formula)])
    ea[[category]] <- setNames(lapply(idents, function(i) as.data.frame(clusterProfiler::enricher(as.numeric(unname((de %>% dplyr::filter(cluster == i))$id)), TERM2GENE = terms))), idents)
    ea[[category]] <- new("compareClusterResult", fun = category,
                          compareClusterResult = plyr::ldply(ea[[category]], "rbind") %>% dplyr::rename(Cluster=.id) %>% dplyr::mutate(Cluster = factor(Cluster, levels=names(ea[[category]])))  )
  }

  for(i in seq_along(ea)) {
    # Patch to fix bug of duplicated descriptions but different IDs in Pathways package 
    problematic.descriptions <- ea[[i]]@compareClusterResult %>% 
      dplyr::select(ID, Description) %>% 
      dplyr::group_by_all() %>%
      dplyr::filter(length(unique(ID)) > 1) %>% 
      unique() %>%
      dplyr::mutate(desc_id = paste0(Description, " (ID: ", ID, ")"))
    
    ea[[i]]@compareClusterResult <- plyr::join(ea[[i]]@compareClusterResult, problematic.descriptions, by = c("Description","ID")) %>% 
      dplyr::mutate(Description = ifelse(is.na(desc_id), Description, desc_id))
  }
  return(ea)
}

get.pathway.genes <- function(pathway_id) {
  if (str_detect(pathway_id, "GO")) { # GO
    library(org.Hs.eg.db)
    genes <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=pathway_id, columns="SYMBOL") %>% pull(SYMBOL) %>% unique
  } else if (str_detect(pathway_id, "hsa")) { # KEGG
    ids <- clusterProfiler::download_KEGG("hsa")$KEGGPATHID2EXTID %>% filter(from %in% pathway_id) %>% pull(to)
    genes <- GeneIdMapping()$names[ids] %>% unname()
  } else if (str_detect(pathway_id, "R-HSA")) { # reactome
    library(reactome.db)
    ids <- AnnotationDbi::select(reactome.db, keytype="PATHID", keys=pathway_id, columns="ENTREZID") %>% pull(ENTREZID) %>% unique
    genes <- GeneIdMapping()$names[ids] %>% unname()
  } else {
    stop(glue::glue("unkown pathway id template {pathway_id}"))
  }

  return(genes)
}

get.pathway.genes <- memoize_first(get.pathway.genes)


PathwayGenes <- function(pathways = NULL, organism = c("hsa", "mmu"), funs = c("KEGG","Pathway","GO:BP","GO:MF","GO:CC")) {
  require(dplyr)

  organism = match.arg(organism)
  funs = match.arg(funs, several.ok = T)

  orgdbs = list(hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db")
  sets = list(
    KEGG = function(...) tryCatch({
      suppressMessages(clusterProfiler::download_KEGG(organism)$KEGGPATHID2EXTID)
    }, error=function(e) data.frame(a=1,b=1)[0,]),
    Pathway = function(...) tryCatch({
      suppressMessages(AnnotationDbi::select(reactome.db::reactome.db, keys=pathways, columns =c("ENTREZID"), keytype = "PATHID") %>% filter(!is.null(ENTREZID)))
    }, error=function(e) data.frame(a=1,b=1)[0,]),
    GO = function(domain, ...) {
      tryCatch({
        terms <- pathways
        if(is.null(terms))
          terms <- GO.db::GOTERM
        k = names(Filter(function(v) v == gsub("^.*:", "", domain), AnnotationDbi::Ontology(terms)))
        suppressMessages(AnnotationDbi::mapIds(GOSemSim::load_OrgDb("org.Hs.eg.db"),
                                               keys = k, keytype = "GOALL",
                                               column="ENTREZID", multiVals = 'list') %>% stack() %>% dplyr::select(2,1))
      }, error=function(e) data.frame(a=1,b=1)[0,])
    })

  return(do.call(rbind, lapply(funs, function(n)
    sets[[gsub(":.*$", "", n)]](n) %>% `colnames<-`(c("pathway","gene")) %>%
      filter(is.null(pathways) | pathway %in% pathways)
    )))
}



count_words <- function(terms, ngram.n=1:2) {
  # Copied and adjusted from simplifyEnrichment::count_words and https://stackoverflow.com/a/53226662/6400526
  library(tm)
  .tokenizer <- function(x)
    unlist(lapply(ngrams(stemDocument(words(x)), ngram.n), paste, collapse = "_"), use.names = FALSE)

  docs = VCorpus(VectorSource(terms))
  docs = tm_map(docs, content_transformer(tolower))
  docs = tm_map(docs, removeNumbers)
  docs = tm_map(docs, removeWords, stopwords())
  docs = tm_map(docs, removePunctuation)
  docs = tm_map(docs, stripWhitespace)
  docs = tm_map(docs, removeWords, NULL)

  tdm = TermDocumentMatrix(
    docs,
    control = list(
      wordLengths = c(1, Inf),
      tokenize = .tokenizer,
      stemming = FALSE,
      dictionary = NULL,
      tolower = FALSE,
      weighting = function(x) weightSMART(x, "ntn")
    )
  )

  return(list(terms=sort(slam::col_means(tdm) %>% `names<-`(names(terms)), decreasing = T),
              words=sort(slam::row_sums(tdm) %>% `names<-`(rownames(tdm)), decreasing = T)))
}

RankClusteredPathways <- function(membership, adjacency, method=c("betweenness","closeness","hub.score","eigen_centrality","page.rank","attributes","semantics"),
                                  attributes=NULL, descriptions=NULL, ties.method="first",
                                  n.select = function(n) as.numeric(cut(n, c(0, 5, 10, 25, 50, Inf), c(1, 2, 4, 6, 8))),
                                  ...) {
  library(igraph)

  method = match.arg(method, several.ok = T)
  if(hasArg(attributes) & !"attributes" %in% method)
    method <- c(method, "attributes")

  if(hasArg(descriptions) & !"semantics" %in% method)
    method <- c(method, "semantics")


  g <- graph_from_adjacency_matrix(adjacency, mode="undirected", weighted = T)
  V(g)$membership <- membership

  # Rank pathways in each pathway cluster independently
  scores <- lapply(unique(membership), function(cluster) {
    sub  <- igraph::induced.subgraph(g, vids = names(V(g))[V(g)$membership == cluster])
    args <- modifyList(list(), list(directed=F, graph=sub))

    # Run different pathways' scoring methods
    do.call(cbind, sapply(method, function(m)
      switch(m,
             "betweenness"      =,
             "closeness"        =R.utils::doCall(match.fun(m), args=args),
             "hub.score"        =,
             "eigen_centrality" =,
             "page.rank"        =R.utils::doCall(match.fun(m), args=args)$vector,
             "attributes"       =attributes[names(V(sub)),],
             "semantics"        ={
               descs <- descriptions[names(V(sub)),] %>% `names<-`(names(V(sub)))
               count_words(descs)$terms / sapply(strsplit(descs, " "), length)
             })
    ,simplify = F)) %>% data.frame %>% mutate(membership=as.character(cluster))
  }) %>% do.call(rbind, .)

  score.cols <- setdiff(colnames(scores), c("membership"))
  scores %>% rownames_to_column() %>%
    group_by(membership) %>%
    mutate(across(score.cols, list(r= ~ base::rank(-., ties.method = ties.method)), .names="{col}.rank")) %>%
    ungroup() %>%

    rowwise() %>%
    mutate(mean.rank = mean(c_across(paste0(score.cols,".rank")))) %>%

    group_by(membership) %>%
    mutate(selected = base::rank(mean.rank, ties.method = "first") <= n.select(n())) %>%
    ungroup() %>%

    column_to_rownames() %>%
    dplyr::select(membership, everything()) %>%
    `[`(rownames(adjacency),)
}


ClusterPathways <- function(pathways.gl,
                            adjacency.args = list(method = "kappa"),
                            clustering.args = list(method = "binary_cut"),
                            rank.pathways.args = list(method = c("page.rank","semantics")),
                            ...) {
  # If not provided set default argument values
  defaults <- list(
    adjacency.args = list(method = "kappa"),
    clustering.args = list(method = "binary_cut"),
    rank.pathways.args = list(method = c("page.rank","semantics"),
                              n.select = function(n) as.numeric(cut(n, c(0, 5, 10, 25, 50, Inf), c(1, 2, 4, 6, 8))))
  )
  for(n in names(defaults))
    assign(n, modifyList(defaults[[n]], get(n)))

  # Compute adjacency matrices
  adjacencies <- list(
    # Adjacency by DEGs
    degs = do.call(simplifyEnrichment::term_similarity,
                   modifyList(adjacency.args, list(
                     gl = pathways.gl %>% tidyr::separate_rows(geneID, sep="/") %>% unique %>% unstack(geneID~ID))))
    ,
    # "Apriori" similarity - adjacency by pathways' genes
    genes = do.call(simplifyEnrichment::term_similarity,
                    modifyList(adjacency.args, list(
                      gl = PathwayGenes(pathways.gl$ID) %>% dplyr::select(2,1) %>%
                        filter_all(Negate(is.na)) %>% unique %>% unstack() )))
  ) %>%

    # Align rows and columns of adjacency matrices
    lapply(., function(a) {
      missing <- setdiff(pathways.gl$ID, rownames(a))
      a <- rbind(a, matrix(0, nrow = length(missing), ncol = ncol(a)) %>% `rownames<-`(missing))
      a <- cbind(a, matrix(0, ncol = length(missing), nrow = nrow(a)) %>% `colnames<-`(missing))

      a[is.na(a)] <- 0
      a[pathways.gl$ID, pathways.gl$ID]
    })

  # Element-wise mean of adjacency matrices
  adjacency <- base::Reduce("+", adjacencies) / 2
  adjacency[is.na(adjacency)] <- 0

  if(nrow(pathways.gl) == 1) {
    membership <- setNames(1:nrow(pathways.gl), rownames(pathways.gl))
    ranks <- NULL
  }
  else {
    # Cluster pathways using binary cut
    membership <- setNames(rep("-", nrow(adjacency)), rownames(adjacency))
    tryCatch(membership <- do.call(simplifyEnrichment::cluster_terms,
                                   modifyList(clustering.args, list(mat = adjacency))) %>%
               `names<-`(rownames(adjacency)),
             error = function(e) warning(e))

    # Rank pathways within clusters
    ranks <- do.call(RankClusteredPathways,
                     modifyList(rank.pathways.args, list(membership = membership,
                                                         adjacency = adjacency,
                                                         descriptions = pathways.gl %>% dplyr::select(Description))))
    ranks <- ranks %>% mutate(Description = pathways.gl[rownames(.),"Description"])
  }

  return(list(
    membership  = membership,
    adjacencies = modifyList(adjacencies, list(joint=adjacency)),
    ranks       = ranks
  ))
}


pathways_clustering <- function(pts, topic, k, print_cluster_heatmap=F, dis.method="correlation") {
  pts <- pt
  pts <- pts %>%
    filter(.data[['topic']] == .env[['topic']])
  
  rownames(pts) <- pts$ID
    
  
  topic_data_tmp <- pts$geneID %>%strsplit(.,'/') %>% `names<-`(c(pts$ID)) %>% stack()
  topic_data_tmp <- table(topic_data_tmp$values,topic_data_tmp$ind) %>% as.data.frame %>%
    tidyr::spread(., key = "Var2", value = "Freq", fill = 0) %>%
    tibble::column_to_rownames("Var1")
  
  ph <- pheatmap::pheatmap(topic_data_tmp, clustering_distance_rows = dis.method, 
                           clustering_distance_cols = dis.method, )
  
  clusters <- cutree(ph$tree_col, k=k) %>% as.factor()
  clusters_tb <- as.data.frame(clusters)
  
  colnames(clusters_tb) <- c("cluster")
  
  if(print_cluster_heatmap){
    pheatmap::pheatmap(topic_data_tmp,
                       clustering_distance_rows = dis.method,
                       clustering_distance_cols = dis.method,
                       annotation_col = clusters_tb)
  }
  
  map.pID <- pts %>% dplyr::select(c("ID", "Description")) %>% unique()
  pts.by.clusters <- split(names(clusters), clusters)
  pts <- pts %>% mutate(cluster="na")
  for (i in names(pts.by.clusters)){
    pts[pts.by.clusters[[i]],]$cluster <- i
  }
  
  return(pts)
}


compute_pathways_distance_matrix <- function (pts, method='spearman') {
  df <- pts$geneID %>%
    strsplit(.,'/') %>%`names<-`(c(pts$ID)) %>% stack()

  gene_pathways_mat <- table(df$values,df$ind) %>%
    as.data.frame %>%
    tidyr::spread(., key = "Var2", value = "Freq", fill = 0) %>%
    tibble::column_to_rownames("Var1")


  correlation_matrix <- cor(gene_pathways_mat, method = method)
  dis <- 1 - correlation_matrix
  return(dis)
}

plot_pathways_overview_graph <- function(pts, k, max_distance=0.8, seed=1, named_nodes = c(),
                                         group_names=c(),
                                         rasterise_nodes_and_edges=FALSE,
                                         rasterise_edges=FALSE,
                                         label_alpha=0.6, ...) {
  pts <- pts %>% filter(topic == k) %>% filter(!is.na(group)) %>% mutate(group_name=group_names[group])

  nodes <- pts %>%
    mutate(group=as.factor(group), group_name=group_name, label=ifelse(Description %in% named_nodes, Description, "")) %>%
    dplyr::select(ID, group, log.p.adjust, label, group_name)

  corr <- compute_pathways_distance_matrix(pts)

  edges <- melt(corr, value.name = 'weight')
  colnames(edges) <- c("from", "to", "weight")
  edges <- edges %>%
    filter(from != to & weight > 0 & weight < max_distance)

  graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  isolated_nodes <- V(graph)[degree(graph) == 0]$name
  nodes_filtered <- nodes[!(nodes$ID %in% isolated_nodes), ]
  edges_filtered <- edges[edges$from %in% nodes_filtered$ID & edges$to %in% nodes_filtered$ID, ]

  g_filtered <- graph_from_data_frame(d = edges_filtered, vertices = nodes_filtered, directed = FALSE)

  custom_palette <- c("#2166ac", "#bf812d", "#dfc27d", "#f6e8c3", "#d1e5f0",
                      "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#8c510a"
                      )

  
  set.seed(seed)
  g <- ggraph(g_filtered, layout = 'fr') +
    geom_edge_link2(aes(edge_alpha = weight), show.legend = FALSE)

  if (rasterise_edges) {
    g <- ggrastr::rasterise(g)
  }

  g <- g + geom_node_point(aes(fill = group_name, size = log.p.adjust), shape = 21, stroke=0.2, color = alpha("black", 0.01), alpha=1)

  if (rasterise_nodes_and_edges) {
    g <- ggrastr::rasterise(g)
  }
  
  g + 
    geom_label_repel(aes(x=x, y=y, label = label), segment.colour='black', 
                     min.segment.length = 0,
                     fill='white', label.size = NA,
                     alpha=label_alpha,
                     segment.alpha=1,
                     point.padding=unit(1, "lines"),
                     box.padding=unit(2, "lines"),
                    max.overlaps = Inf
                     ) +
    scale_fill_manual(values = custom_palette) +
    scale_size_continuous(range = c(5, 13))+
    theme_graph(base_family = 'Helvetica')

}
