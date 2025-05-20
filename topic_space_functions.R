library(dplyr)
library(reticulate)
library(dendextend)

source("utils.R")

AD.traits        <- c(sqrt.amyloid_mf="neocortical amyloid", sqrt.tangles_mf="neocortical tangles", cogng_demog_slope="cognitive decline rate")
AD.traits.colors <- c(sqrt.amyloid_mf="midnightblue",     sqrt.tangles_mf="olivedrab4",       cogng_demog_slope="firebrick3")

cyan2orange.correlation <- circlize::colorRamp2(c(-1,0,1), c("#4db6ac","white","#ef6c00"))

# Topic Space ####

CELL_TYPE_TO_K <- list(
  "inhibitory"="17",
  "oligo"="8",
  "microglia"="15",
  "opcs"="7",
  "cux2p"="14",
  "cux2m"="16",
  "astrocytes"="10"
)

CELL_TYPE_DIRECTORY <- list(
  inhibitory='InhibitoryAll',
  astrocytes='AstrocytesAll',
  oligo='OligoAllv1.1',
  opcs='OPCAllv1.1',
  microglia='MicrogliaAll',
  cux2p='ExcitatoryCux2+All',
  cux2m='ExcitatoryCux2-All'
)

read_topic_space <- function(join_by = 'projid', result_dir = get_result_dir()) {
  # TODO: read from h5 instead
  oligo <- read.csv1and2(file.path(result_dir, "OligoAll/8/mean_topic_scores.csv"))
  astrocytes <- read.csv1and2(file.path(result_dir, "AstrocytesAll/10/mean_topic_scores.csv"))
  microglia <- read.csv1and2(file.path(result_dir, "MicrogliaAll/15//mean_topic_scores.csv"))
  opcs <- read.csv1and2(file.path(result_dir, "OPCAll/7/mean_topic_scores.csv"))
  excitatory_cuxp <- read.csv1and2(file.path(result_dir, "ExcitatoryCux2+All", "14", "mean_topic_scores.csv"))
  excitatory_cuxm <- read.csv1and2(file.path(result_dir, "ExcitatoryCux2-All/16/mean_topic_scores.csv"))
  inhibitory_all <- read.csv1and2(file.path(result_dir, "InhibitoryAll", "17", "mean_topic_scores.csv"))

  oligo %>%
    dplyr::left_join(astrocytes, by = join_by) %>%
    dplyr::left_join(microglia, by = join_by) %>%
    dplyr::left_join(opcs, by = join_by) %>%
    dplyr::left_join(inhibitory_all, by = join_by) %>%
    dplyr::left_join(excitatory_cuxp, by = join_by) %>%
    dplyr::left_join(excitatory_cuxm, by = join_by) %>%
    mutate_if(is.numeric,coalesce,0)
}

# DELETE
save_topics_pathways <- function (h5_path=TOPICS_DATA_H5_PATH, topics_k=CELL_TYPE_TO_K) {
  for (cell_type in names(topics_k)) {
    cell_type_dir <- CELL_TYPE_DIRECTORY[[cell_type]]
    k <- topics_k[[cell_type]]
    dataset_path <- paste0(v, "/pathways")
    upper <- ifelse(str_detect(dataset_path, "glia"), "080", "095")
    pathways_path <- file.path(cell_type_dir, k, glue::glue("filtered_de/de_gene_ontology_{upper}.csv"))
    pt <- read.csv1and2(file.path(get_result_dir(), pathways_path))
    pt$cell.type <- cell_type

    attr(pt, "upper_quantile") <- upper
    attr(pt, "path") <- pathways_path
    
    print(glue::glue("Writing pathways of {cell_type} ({pathways_path}) to h5 {dataset_path}"))
    h5delete(h5_path, dataset_path)
    h5write(pt, h5_path, dataset_path, write.attributes=TRUE)
  }
}


save_topics_de <- function (h5_path=TOPICS_DATA_H5_PATH, topics_k=CELL_TYPE_TO_K) {
  for (cell_type in names(topics_k)) {
    cell_type_dir <- CELL_TYPE_DIRECTORY[[cell_type]]
    k <- topics_k[[cell_type]]
    dataset_path <- paste0(CELL_TYPE_TO_GROUP[[cell_type]], "/de")
    upper <- ifelse(str_detect(dataset_path, "glia"), "080", "095")
    de_path <- file.path(cell_type_dir, k, glue::glue("filtered_de/filtered_de_{upper}.csv"))
    de <- read.csv1and2(file.path(get_result_dir(), de_path), row.names=1)
    de$cell.type <- cell_type

    attr(de, "upper_quantile") <- upper
    attr(de, "path") <- de_path

    print(glue::glue("Writing de of {cell_type} ({de_path}) to h5 {dataset_path}"))
    tryCatch({
        h5delete(h5_path, dataset_path)
      }, error = function(e) {
    })

    h5write(de, h5_path, dataset_path, write.attributes=TRUE)
  }
}


save_topic_fit <- function(fit, dataset_path, h5_path=TOPICS_DATA_H5_PATH) {
  groups <- h5ls(h5_path)$group
  if (!(dataset_path %in% groups) & !(paste0("/", dataset_path) %in% groups)) {
    h5createGroup(h5_path, dataset_path)
  }

  loading <- fit[['L']] %>% as.data.frame() %>% rownames_to_column("cell")

  tryCatch({
    h5delete(h5_path, paste0(dataset_path, "/L"))
  }, error = function(e) {
  })

  h5write(loading, h5_path, paste0(dataset_path, "/L"))
  print("writing F")
  f <- fit[['F']] %>% as.data.frame() %>% rownames_to_column("gene")

  tryCatch({
    h5delete(h5_path, paste0(dataset_path, "/F"))
  }, error = function(e) {
  })

  h5write(f, h5_path, paste0(dataset_path, "/F"))
}


save_topics_score <- function(h5_path=TOPICS_DATA_H5_PATH, topics_k=CELL_TYPE_TO_K, cell_types=names(topics_k)) {
  FITS_PATH <- list(
    inhibitory="inhibitory_topics/all_inhibitory_topic_fit",
    astrocytes="astrocytes_topics/all_astrocytes_topic_fit",
    oligo="oligo_topics/all_oligo_topic_fit",
    opcs="opcs_topics/all_opcs_topic_fit",
    microglia="microglia_topics/all_microglia_topic_fit",
    cux2p="excitatory_topics/all_cux2_plus_excitatory_topic_fit",
    cux2m="excitatory_topics/all_cux2_minus_excitatory_topic_fit"
  )

  SUBSET_FITS_PATH <- list(
    inhibitory="inhibitory_topics/inhibitory_topic_fit",
    astrocytes="astrocytes_topics/astrocytes_topic_fit",
    microglia="microglia_topics/microglia_topic_fit"
  )

  ks_range <- c(-1, 1)
  for (cell_type in cell_types) {
    for (k_diff in ks_range) {
      for (all_cells in c(TRUE, FALSE)) {
        # Do we have partial fits for this cell type?
        if  (!all_cells && !cell_type %in% names(SUBSET_FITS_PATH)) {
          next
        }
        
        # Save different ks only for the all fits
        if (!all_cells && 0 != k_diff) {
          next
        }
  
        fits_path <- ifelse(all_cells, FITS_PATH[[cell_type]], SUBSET_FITS_PATH[[cell_type]])
        k <- topics_k[[cell_type]]
        dataset_path <- CELL_TYPE_TO_GROUP[[cell_type]]
  
        if (!all_cells) {
          dataset_path <- paste0(dataset_path, "/random_subset")
          h5createGroup(h5_path, dataset_path)
        } else if (0 != k_diff) {
          k <- as.character(as.integer(k) + k_diff)
          dataset_path <- paste0(dataset_path, "/fit_", k)
          h5createGroup(h5_path, dataset_path)
        }
  
        print(glue::glue("Writing loading and F of {cell_type} ({fits_path}) to h5 {dataset_path}"))
        
        topic_fits <- read_topic_models(file.path(get_result_dir(), fits_path), k=k)
        stopifnot(k %in% names(topic_fits))

        loading <- topic_fits[[k]][['L']] %>% as.data.frame() %>% rownames_to_column("cell")
  
        tryCatch({
          h5delete(h5_path, paste0(dataset_path, "/L"))
        }, error = function(e) {
        })
  
        h5write(loading, h5_path, paste0(dataset_path, "/L"))
        print("writing F")
        f <- topic_fits[[k]][['F']] %>% as.data.frame() %>% rownames_to_column("gene")
  
        tryCatch({
          h5delete(h5_path, paste0(dataset_path, "/F"))
        }, error = function(e) {
        })
  
        h5write(f, h5_path, paste0(dataset_path, "/F"))
  
      }
    }
  }
}


# Correlation ####

correlate_prevalence <- function (X) {
  cor <- list(
    names = colnames(X),
    corr = stats::cor(X, use = "pairwise.complete.obs", method = "spearman"),
    pval = outer(1:ncol(X), 1:ncol(X), Vectorize(function(i,j)
      cor.test(X[,i], X[,j], use="pairwise.complete.obs", method = "spearman")[["p.value"]]))
  )

  cor$adj.pval <- matrix(p.adjust(cor$pval, method = "BH"), nrow=nrow(cor$pval))
  cor$sig <- matrix(cut(cor$adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")), nrow=nrow(cor$pval))

  cor$params <- list(cor.method = "spearman",
                                 cor.use = "pairwise.complete.obs",
                                 p.adjust.method = "BH")

  cor
}

plot_landscape_correaltion <- function(cor, title, names=cor$names, add_significance=TRUE,
                                       col = cyan2orange.correlation,
                                       show_column_names = T,
                                       ...) {
  corr <- cor$corr %>% `dimnames<-`(list(cor$names, cor$names))
  corr.sig <- cor$sig %>% `dimnames<-`(list(cor$names, cor$names))

  corr <- corr[names, names]
  corr.sig <- corr.sig[names, names]
  
  if (add_significance) {
    cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y, gp = grid::gpar(fontsize=6, fontface="bold"))
  } else {
    cell_fun <- NULL
  }

  cors.heatmap <- ComplexHeatmap::Heatmap(corr,
                                          name = title,
                                          column_title_gp = gpar(fontsize=10, fontface="bold"),
                                          heatmap_legend_param = list(title=title, direction="horizontal"),
                                          col = col,
                                          clustering_distance_rows = "spearman",
                                          clustering_distance_columns = "spearman",
                                          column_names_rot = 90,
                                          show_row_names = T,
                                          show_column_names = show_column_names,
                                          cell_fun = cell_fun, ...)

  cors.heatmap
}


# Communities #######

build_communities <- function(data, features, n.sub.communities = NULL) {
  features <- features[order(match(features, data$var_names))]
  
  dynamics <- py_to_r(data$uns$trajectories$palantir$dynamics$pred.vals) %>% filter(feature %in% features)
  correlations <- data$uns$tt.cor$corr %>% `dimnames<-`(list(data$uns$tt.cor$names, data$uns$tt.cor$names))
  correlations <- correlations[features, features]
  
  communities <- construct.communities(
    dynamics     = dynamics, 
    correlations = correlations,
    k=10,
    weights = c(1,-1,1),
    partition.type = "RBERVertexPartition")
  
  communities$membership <- communities$membership[order(match(names(communities$membership), data$var_names))]
  communities$membership <- recode(communities$membership, "0"="C2", "1"="C1", "2"="C3")
  
  dends <- split(features, communities$membership) %>%
    lapply(., function(states) {
      lapply(c("dynamics.adjacency", "correlation.mtx"), function(n)
        communities[[n]][states, states]) %>%
        base::Reduce(`+`, .) %>%
        dist %>%
        hclust %>%
        dendsort::dendsort(.) %>%
        as.dendrogram() %>%
        ladderize %>%
        set("labels_to_character")
      
    })
  
  n.communities <- communities$membership %>% unique %>% length
  if (n.communities > 3) {
    warn("Number of communities is larger than 3")
  }
  
  
  if (is.null(n.sub.communities)) {
    n.sub.communities <- rep(2,n.communities)
  } else if (length(n.sub.communities) == 1) {
    n.sub.communities <- rep(n.sub.communities ,n.communities)
  } # Else use the supplied n.sub.communities
  
  sub.dends <- lapply(1:n.communities, function(i, comm, nsplit)
    cutree(dends[[comm[[i]]]], nsplit[[i]]) %>%
      split(names(.), .) %>%
      lapply(., function(s) prune(dends[[comm[[i]]]], leaves = setdiff(labels(dends[[comm[[i]]]]), s))) %>%
      `names<-`(paste0(comm[[i]], ".", names(.))),
    comm = paste0("C", 1:n.communities), nsplit = n.sub.communities) %>% unlist(., recursive=F)
  
  membership <- 
    data.frame(community = communities$membership, row.names = features) %>%
    merge(., lapply(sub.dends, labels) %>% stack, by.x="row.names", by.y="values", all.x=TRUE) %>%
    column_to_rownames("Row.names") %>% 
    dplyr::select(community, sub.community=ind) %>%
    `[`(match(data$var_names,rownames(.)),)

  data$uns$communities <- list(
    similarities = communities %>% `[`(c("dynamics.mtx","dynamics.adjacency","correlation.mtx")) %>% 
      `names<-`(c("dynamics","dynamics.adjacency","correlation")),
    names = rownames(communities$dynamics.mtx),
    dynamics.colnames = colnames(communities$dynamics.mtx),
    resolution.parameters = communities$resolution.parameters
  )
  
  data$var <- membership %>% mutate(community=factor(community)) %>% `rownames<-`(colnames(data))
  
  data$obsm$communities <- donor.community.proportion(data$X, data$var %>% dplyr::select(community))
  data$obsm$sub.communities <- donor.community.proportion(data$X, data$var %>% dplyr::select(sub.community))
  
  
  data$uns$communities$dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                                                features = cbind(data$obsm$communities %>% dplyr::select(starts_with("C")),
                                                                 data$obsm$sub.communities %>% dplyr::select(starts_with("C"))),
                                                trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                                                trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                                                bootstrap=F)
  
  data$uns$communities$trait.association <- associate.traits(
    traits = data$obsm$meta.data[,c("sqrt.amyloid_mf","sqrt.tangles_mf","cogng_demog_slope")],
    covariates = sqrt(cbind(data$obsm$communities %>% dplyr::select(starts_with("C")),
                            data$obsm$sub.communities %>% dplyr::select(starts_with("C")))),
    controls = data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")]))
  
  data
}


plot_communities <- function(data, communities_plots_dir) {
  dir.create(communities_plots_dir, showWarnings = FALSE)
  
  
  withr::with_pdf(file.path(communities_plots_dir, "TraitAssociationCommunities.pdf"), {
    plot.trait.associations.heatmap(data$uns$communities$trait.association, row.by = 'covariate')  
  })
  
  communities <- data$var %>% filter(!is.na(community)) %>% pull(community) %>% unique
  all.sub.communities <- data$var %>% filter(!is.na(sub.community)) %>% pull(sub.community) %>% unique
  
  withr::with_pdf(file.path(communities_plots_dir, "CommunitiesTrajectories.pdf"), {
    print(plot.dynamics.wrapper(communities, dynamics = data$uns$communities$dynamics, include.points=F,
                                labels=topics_to_states_annoation,
                                legend.position = "none", overlap.pseudotime = .3) + labs(x=NULL, y=NULL, title='Communities'))
    
    print(plot.dynamics.wrapper(all.sub.communities %>% unique, dynamics = data$uns$communities$dynamics, include.points=F,
                                labels=topics_to_states_annoation,
                                legend.position = "none", overlap.pseudotime = .3) + labs(x=NULL, y=NULL, title="Sub Communities"))
    
    for (c in communities) {
      message(glue::glue("Plotting dynamics for community {c}"))
      sub.communities <- data$var %>% dplyr::filter(community == c) %>% pull(sub.community) %>% unique %>% S4Vectors::unfactor(.)
      community_and_sub_communities <- c(c, sub.communities)
      print(plot.dynamics.wrapper(community_and_sub_communities, dynamics = data$uns$communities$dynamics, include.points=F,
                                  labels=topics_to_states_annoation,
                                  legend.position = "none", overlap.pseudotime = .3) + labs(x=NULL, y=NULL, title=glue::glue("Community {c} vs sub communities")))
      
      print(plot.dynamics.wrapper(data$var %>% filter(community == c) %>% rownames(), 
                                  dynamics =  data$uns$trajectories$palantir$dynamics, include.points=F,
                                  labels=topics_to_states_annoation,
                                  legend.position = "none", overlap.pseudotime = .3) + labs(x=NULL, y=NULL, title=glue::glue("Community {c} topics")))
      
      for (sc in sub.communities) {
        print(plot.dynamics.wrapper(data$var %>% filter(sub.community == sc) %>% rownames(), 
                                    dynamics =  data$uns$trajectories$palantir$dynamics, include.points=F,
                                    labels=topics_to_states_annoation,
                                    legend.position = "none", overlap.pseudotime = .3) + labs(x=NULL, y=NULL, title=glue::glue("Sub Community {sc} topics")))
      }
      
    }
    
    
  })
  
  write.csv(data$var %>% dplyr::select(community, sub.community) %>% rownames_to_column("topic"), file.path(communities_plots_dir, "communities.csv"), row.names = F)
}


topic_quantiles_healthy <- function(df, topic, q = c(0.5,0.8,0.9,0.95, 0.99)) {
  q_start <- quantile(df %>% filter(group == 'start') %>% pull(topic), q) %>% data.frame(threshold=.)
  q_start$prAD_quantile <- sapply(q_start %>% pull(threshold), function (threshold) {
    below_threshold <- mean((df %>% filter(group == 'prAD') %>% pull(topic)) < threshold)

    return(as.integer(below_threshold * 100))
  })
  q_start
}


topics_to_states_annoation <- c(
  c("Ast.10.k6"="AST.p6 (Disease)", "Ast.10.k2"="AST.p2 (ABA)", "Ast.10.k5"="AST.p5 (Homestatic)"),
  c("Micro.15.k9"="MIC.p9 (Homestatic)", "Micro.15.k15"="MIC.p15 (Reactive)",
   "Micro.15.k4"="MIC.p4 (Monocytes/Macrophages)",
    "Micro.15.k3"="MIC.p3 (Proliferating)"),
  c("Oligo.8.k8"="OLI.p8 (Disease)"),
  c("OPC.7.k4"="OPC.p4 (ABA)", "OPC.7.k2"="OPC.p2 (COP)")

)

topic_names_mapping <- function(names, use_static_mapping=FALSE) {
  if (use_static_mapping) {
    static_mapping <- c(unlist(unname(neuronal_topics_to_annotation)),
                        topics_to_states_annoation)
  } else {
    static_mapping <- c()
  }


  sapply(names, function(topic_name) {
    k <- str_match(topic_name, "k(.*)")[[2]]
    ct <- case_when(str_detect(topic_name, "Inh") ~ "INH",
                    str_detect(topic_name, "Exc.Cux2P") ~ "EXC L2-3",
                    str_detect(topic_name, "Exc.Cux2M") ~ "EXC L4-6",
                    str_detect(topic_name, "Micro") ~ "MIC",
                    str_detect(topic_name, "Oligo") ~ "OLI",
                    TRUE ~ str_to_upper(str_match(topic_name, "(.*?)\\.")[[2]]))

    ifelse(
      topic_name %in% names(static_mapping),
      static_mapping[[topic_name]],
      glue::glue("{ct}.p{k}")
    )
  }, simplify = TRUE, USE.NAMES = F)
}
