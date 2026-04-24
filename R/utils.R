#' @title Utility Functions
#' @description General utility functions for data loading, path management,
#'              similarity metrics, and data manipulation.

library(dplyr)


#' Read CSV file with automatic delimiter detection
#'
#' @description Reads CSV files, automatically detecting whether the delimiter
#'              is comma (read.csv) or semicolon (read.csv2).
#'
#' @param path Path to CSV file
#' @param ... Additional arguments passed to read.csv/read.csv2
#'
#' @return Data frame
read.csv1and2 <- function(path, ...) {
  line <- readLines(path , n = 1)
  if (strsplit(line, ";")[[1]] %>% length > 1) {
    data <- read.csv2(path, ...)
  } else {
    data <- read.csv(path, ...)
  }

  data
}

#' Count number of cells expressing each gene
#'
#' @description Counts the number of cells with non-zero expression for each gene
#'              in a Seurat object's SCT assay.
#'
#' @param object Seurat object with SCT assay
#'
#' @return Named numeric vector of cell counts per gene
number_of_expressing_cells <- function(object) {
  gene_counts <- tabulate(object[['SCT']]@counts@i + 1)
  names(gene_counts) <- rownames(object[['SCT']]@counts)
  gene_counts
}

#' Factor-preserving case_when
#'
#' @description Wrapper around dplyr::case_when that returns a factor with
#'              levels ordered by their appearance in the formula RHS values.
#'
#' @param ... Formula pairs passed to case_when
#'
#' @return Factor with levels in order of appearance
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]

  factor(dplyr::case_when(...), levels=levels)
}

#' Get results directory path
#'
#' @description Returns the path to the results/plots directory, automatically
#'              detecting whether running on cluster or local mount.
#'
#' @return Character string path to results directory
get_result_dir <- function () {
  CLUSTER_PLOTS_DIR <- "/ems/elsc-labs/habib-n/roi.meir/500/Plots/"
  LOCAL_CLUSTER_MOUNT_DIR <- "/Volumes/habib-lab/roi.meir/500/Plots"

  if (dir.exists(CLUSTER_PLOTS_DIR)) {
    CLUSTER_PLOTS_DIR
  } else {
    LOCAL_CLUSTER_MOUNT_DIR
  }
}

#' Get shared data directory path
#'
#' @description Returns the path to the shared data directory, automatically
#'              detecting whether running on cluster or local mount.
#'
#' @return Character string path to shared directory
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

#' Compute Jaccard similarity between two sets
#'
#' @param a First set (vector)
#' @param b Second set (vector)
#'
#' @return Numeric Jaccard index (intersection / union)
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

#' Compute pairwise Jaccard similarity matrix
#'
#' @param a List of sets
#' @param b List of sets (default: same as a)
#'
#' @return Matrix of Jaccard similarities
jaccard_correlation <- function(a, b=a) {
  outer(a, b, Vectorize(jaccard))
}

#' Compute overlap coefficient between two sets
#'
#' @description Computes the overlap coefficient (Szymkiewicz-Simpson coefficient):
#'              intersection / min(size_a, size_b)
#'
#' @param a First set (vector)
#' @param b Second set (vector)
#'
#' @return Numeric overlap coefficient (0-1)
overlap_coefficient <- function (a, b) {
  intersection <- length(intersect(a, b))
  return (intersection/ min(length(a), length(b)))
}

#' Compute pairwise overlap coefficient matrix
#'
#' @param a List of sets
#' @param b List of sets (default: same as a)
#'
#' @return Matrix of overlap coefficients
overlap_coefficient_correlation <- function (a, b=a) {
  outer(a, b, Vectorize(overlap_coefficient))
}

#' Jaccard for nested tibble columns
#'
#' @param a Tibble with single list column
#' @param b Tibble with single list column
#'
#' @return Numeric Jaccard index
tibble_jaccard <- function(a, b) {
  jaccard(a[[1]], b[[1]])
}

#' Overlap coefficient for nested tibble columns
#'
#' @param a Tibble with single list column
#' @param b Tibble with single list column
#'
#' @return Numeric overlap coefficient
tibble_overlap_coefficient <- function(a, b) {
  overlap_coefficient(a[[1]], b[[1]])
}


#' Compare groups using Jaccard or overlap coefficient
#'
#' @description Computes pairwise similarity between groups defined in data frames.
#'
#' @param df1 First data frame
#' @param df2 Second data frame (default: NULL = same as df1)
#' @param group_column_1 Column name for grouping in df1
#' @param id_column Column name for IDs to compare
#' @param group_column_2 Column name for grouping in df2 (default: same as group_column_1)
#' @param method Comparison function (default: tibble_jaccard)
#'
#' @return Matrix of similarities with group names as row/column names
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

#' Plot Jaccard comparison as heatmap
#'
#' @description Creates a ComplexHeatmap visualization of Jaccard comparisons.
#'
#' @param df1 First data frame
#' @param df2 Second data frame (default: NULL)
#' @param group_column_1 Column for grouping in df1
#' @param id_column_1 Column for IDs
#' @param group_column_2 Column for grouping in df2
#' @param column_title Heatmap title (default: "Jaccard comparison")
#' @param cluster_rows Cluster rows (default: FALSE)
#' @param cluster_columns Cluster columns (default: FALSE)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return ComplexHeatmap object
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

#' Memoize a function based on first argument
#'
#' @description Creates a memoized version of a function that caches results
#'              based on the hash of the first argument.
#'
#' @param fun Function to memoize
#'
#' @return Memoized function
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

#' Compute correlations between two data frames with p-values
#'
#' @description Computes all pairwise correlations between columns of X and Y
#'              with significance testing and multiple testing correction.
#'
#' @param X First data frame or matrix
#' @param Y Second data frame or matrix (default: NA = same as X)
#' @param use Method for handling missing values (default: "pairwise.complete.obs")
#' @param method Correlation method: "spearman" or "pearson" (default: "spearman")
#' @param p.adjust.method Method for p-value adjustment (default: "BH")
#'
#' @return List with:
#'         - names: list of column names for X and Y
#'         - corr: correlation matrix
#'         - pval: p-value matrix
#'         - adj.pval: adjusted p-value matrix
#'         - sig: significance symbols (*** < 0.001, ** < 0.01, * < 0.05)
#'         - params: list of method parameters
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


#' @section Neuronal Topic Annotations:
#' Named list mapping trajectory groups to topic-annotation pairs.
#' Groups: decreasing_early, decreasing_late, increasing_early, increasing_late
neuronal_topics_to_annotation <- list(
  decreasing_early=c("Inh.17.k1"="Inh decreasing early","Exc.Cux2P.14.k9"="Exc L2-3 decreasing early","Exc.Cux2M.16.k4"="Exc L4-6 decreasing early"),
  decreasing_late=c("Inh.17.k7"="Inh decreasing late","Exc.Cux2P.14.k10"="Exc L2-3 decreasing late","Exc.Cux2M.16.k8"="Exc L4-6 decreasing late"),
  increasing_early=c("Inh.17.k10"="Inh increasing early","Exc.Cux2P.14.k11"="Exc L2-3 increasing early","Exc.Cux2M.16.k3"="Exc L4-6 increasing early"),
  increasing_late=c("Inh.17.k8"="Inh increasing late","Exc.Cux2P.14.k4"="Exc L2-3 increasing late","Exc.Cux2M.16.k9"="Exc L4-6 increasing late")
)

TOPIC_OVERLLAPING_PSEUDOTIME <- 0.375
GREEN_ET_AL_OVERLLAPING_PSEUDOTIME <- 0.112


#' Add module scores with custom column names
#'
#' @description Wrapper around Seurat::AddModuleScore that renames the output
#'              columns using the feature list names instead of "Cluster1", etc.
#'
#' @param obj Seurat object
#' @param features Named list of gene vectors (names will become column names)
#' @param prefix Prefix to add to column names (default: "")
#'
#' @return Seurat object with renamed module score columns in meta.data
add_module_score_named <- function(obj, features, prefix="") {
  stopifnot(!is.null(names(features)))
  mapping <- setNames(paste0("Cluster", seq_along(features)), paste0(prefix, names(features)))
  obj <- Seurat::AddModuleScore(obj, features = features)
  obj@meta.data <- obj@meta.data %>%
    dplyr::rename(!!!mapping)

  obj
}

#' Remove AsIs class from object
#'
#' @description Removes the "AsIs" class attribute that can be added by I() function,
#'              useful for cleaning data before HDF5 serialization.
#'
#' @param X Object potentially with AsIs class
#'
#' @return Object with AsIs class removed (if present)
unAsIs <- function(X) {
    if("AsIs" %in% class(X)) {
        class(X) <- class(X)[-match("AsIs", class(X))]
    }
    X
}

#' Check if a list is nested
#'
#' @description Checks if a list contains nested lists.
#'
#' @param x List to check
#'
#' @return Logical indicating if list is nested
is_nested <- function(x) {
  if (class(x) != "list") {
    stop("Expecting 'x' to be a list")
  }
  out <- any(sapply(x, is.list))
  return(out)
}

#' Convert inner names of a named list to their root (outer) names
#'
#' @description
#' Given a named list where each element is itself a named vector or list,
#' this function returns a named character vector mapping each inner name
#' to its root (outer) name. The resulting vector has names corresponding
#' to the inner elements and values corresponding to their respective root.
#'
#' @param x A named list, where each element is itself a named vector or list.
#'
#' @return A named character vector. Names are the inner names of each element, 
#' values are the corresponding root (outer) names.
#'
#' @examples
#' x <- list(A = c(gene1=1, gene2=2), B = c(gene3=3))
#' inner_to_root(x)
#' # Returns: 
#' #  gene1  gene2  gene3 
#' #   "A"    "A"    "B"
inner_to_root <- function(x) {
  unlist(
    lapply(names(x), function(root) {
      setNames(rep(root, length(x[[root]])), names(x[[root]]))
    }),
    use.names = TRUE
  )
}
