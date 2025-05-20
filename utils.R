load.expr <- function(cell_types = NULL, data_path="data/genes_exp_by_donor_and_gene"){
  if (is.null(cell_types)) {
    files <- list.files(data_path)
    files <- files[grepl(".csv", files)]
  } else {
    files <- paste0(cell_types, ".csv")
  }

  objs <- lapply(files, function(f_){
    print(f_)
    tidyr::spread(read.table(file.path(data_path, f_)), gene, mean.exp)
  })
  names(objs) <- substr(files, 1, nchar(files) - 4)

  return(format.obj(objs))
}

format.obj <- function(objs){
  return(lapply(objs, function(obj){
    # Remove donor that are NA & Set the rownames to be the donor column
    obj[!(apply(obj, 1, function(x){any(is.na(x))})),]%>% tibble::column_to_rownames("donor")
  }))
}

read.csv1and2 <- function(path, ...) {
  line <- readLines(path , n = 1)
  if (strsplit(line, ";")[[1]] %>% length > 1) {
    data <- read.csv2(path, ...)
  } else {
    data <- read.csv(path, ...)
  }

  data
}

number_of_expressing_cells <- function(object) {
  gene_counts <- tabulate(object[['SCT']]@counts@i + 1)
  names(gene_counts) <- rownames(object[['SCT']]@counts)
  gene_counts
}

fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]

  factor(dplyr::case_when(...), levels=levels)
}

get_result_dir <- function () {
  CLUSTER_PLOTS_DIR <- "/ems/elsc-labs/habib-n/roi.meir/500/Plots/"
  LOCAL_CLUSTER_MOUNT_DIR <- "/Volumes/habib-lab/roi.meir/500/Plots"

  if (dir.exists(CLUSTER_PLOTS_DIR)) {
    CLUSTER_PLOTS_DIR
  } else {
    LOCAL_CLUSTER_MOUNT_DIR
  }
}

get_shared_dir <- function () {
  CLUSTER_SHARED_DIR <- "/ems/elsc-labs/habib-n/Shared/"
  LOCAL_CLUSTER_MOUNT_DIR <- "/Volumes/habib-lab/Shared/"

  if (dir.exists(CLUSTER_SHARED_DIR)) {
    CLUSTER_SHARED_DIR
  } else {
    LOCAL_CLUSTER_MOUNT_DIR
  }
}

SEURAT_OBJECT_DATA_PATH <- file.path(get_shared_dir(), "NextSeq/500/v1.1.objects")
SEURAT_SUBSET_OBJECT_PATH <- file.path(get_shared_dir(), "500/subset")


###########
# Jaccard #
###########

jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard_correlation <- function(a, b=a) {
  outer(a, b, Vectorize(jaccard))
}

overlap_coefficient <- function (a, b) {
  intersection <- length(intersect(a, b))
  return (intersection/ min(length(a), length(b)))
}

overlap_coefficient_correlation <- function (a, b=a) {
  outer(a, b, Vectorize(overlap_coefficient))
}

tibble_jaccard <- function(a, b) {
  jaccard(a[[1]], b[[1]])
}

tibble_overlap_coefficient <- function(a, b) {
  overlap_coefficient(a[[1]], b[[1]])
}


compare_jaccard <- function (df1, df2=NULL, group_column_1, id_column,
                             group_column_2 = group_column_1, method=tibble_jaccard) {
  if (is.null(df2)) {
    df2 <- df1
  }

  df1_nested <- df1 %>% group_by({ {  group_column_1 } }) %>% dplyr::select({{id_column}}, {{group_column_1}}) %>% group_nest()
  df2_nested <- df2 %>% group_by({ {  group_column_2 } }) %>% dplyr::select({{id_column}}, {{group_column_2}}) %>% group_nest()

  mat <- outer(df1_nested %>% pull(data),
               df2_nested %>% pull(data),
               Vectorize(method))

  rownames(mat) <- df1_nested %>% pull({ { group_column_1 } })
  colnames(mat) <- df2_nested %>% pull({ { group_column_2 } })

  mat
}

plot_jaccard_comparison <- function (df1, df2=NULL,
                                     group_column_1, id_column_1,
                                     group_column_2 = group_column_1,
                                     column_title="Jaccard comparison",
                                     cluster_rows = FALSE, cluster_columns = FALSE, ...) {
    mat <- compare_jaccard(df1, df2, group_column, id_column)

    ComplexHeatmap::Heatmap(mat,
                        col = colorRampPalette(c("white", "firebrick4"))(30),
                        show_column_dend = FALSE, show_row_dend = FALSE,
                        show_row_names = T, cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                        column_names_gp = gpar(fontsize=7),
                        row_names_gp = gpar(fontsize=7),
                        column_title=column_title, ...)
}

memoize_first <- function(fun) {
  fun
  cache <- list()
  dec <- function(arg, ...) {
    arg_hash <- rlang::hash(arg)
    if (is.null(cache[[arg_hash]])) cache[[arg_hash]] <<- fun(arg, ...)
    cache[[arg_hash]]
  }
  dec
}

correlate_dataframes <- function (X, Y=NA, use="pairwise.complete.obs", method = "spearman", p.adjust.method = "BH") {
  if (is.na(Y)) Y <- X
  cor <- list(
    names = list(colnames(X), colnames(Y)),
    corr = stats::cor(X, Y, use = use, method = method),
    pval = outer(1:ncol(X), 1:ncol(Y), Vectorize(function(i,j)
      cor.test(X[,i], Y[,j], use=use, method = method)[["p.value"]]))
  )

  cor$adj.pval <- matrix(p.adjust(cor$pval, method = p.adjust.method), nrow=nrow(cor$pval))
  cor$sig <- matrix(cut(cor$adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")), nrow=nrow(cor$pval))

  cor$params <- list(cor.method = method,
                     cor.use = use,
                     p.adjust.method = p.adjust.method)

  cor
}


neuronal_topics_to_annotation <- list(
  decreasing_early=c("Inh.17.k1"="Inh decreasing early","Exc.Cux2P.14.k9"="Exc L2-3 decreasing early","Exc.Cux2M.16.k4"="Exc L4-6 decreasing early"),
  decreasing_late=c("Inh.17.k7"="Inh decreasing late","Exc.Cux2P.14.k10"="Exc L2-3 decreasing late","Exc.Cux2M.16.k8"="Exc L4-6 decreasing late"),
  increasing_early=c("Inh.17.k10"="Inh increasing early","Exc.Cux2P.14.k11"="Exc L2-3 increasing early","Exc.Cux2M.16.k3"="Exc L4-6 increasing early"),
  increasing_late=c("Inh.17.k8"="Inh increasing late","Exc.Cux2P.14.k4"="Exc L2-3 increasing late","Exc.Cux2M.16.k9"="Exc L4-6 increasing late")
)

TOPIC_OVERLLAPING_PSEUDOTIME <- 0.375
GREEN_ET_AL_OVERLLAPING_PSEUDOTIME <- 0.112
