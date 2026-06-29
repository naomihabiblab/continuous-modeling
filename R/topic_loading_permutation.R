#' @title Topic-loading within-cell permutation null
#' @description Reusable functions for testing whether two NMF topic-loading
#'              programs are more correlated (across cells) than expected under a
#'              within-cell topic-label permutation null. Reads loadings (L) from
#'              the topics HDF5 so the test is self-contained and can be run
#'              independently of any Seurat object or fit object.
#'
#' Typical usage (any cell type, after L was written to H5 by save_topics_score):
#'   source("R/topic_loading_permutation.R")
#'   run_increasing_programs_permutation_null(
#'     cell.type  = "inhibitory",
#'     output_dir = file.path(get_result_dir(), "InhibitoryAll", "permutation_null")
#'   )

library(dplyr)
library(ggplot2)

source("R/utils.R")
source("R/data_500_utils.R")

NEURONAL_TYPE_LABELS <- c(
  inhibitory = "Inhibitory",
  cux2p = "Exc upper layers (CUX2+)",
  cux2m = "Exc lower layers (CUX2-)"
)

#' Human-readable label for a neuronal cell type, falling back to the id itself.
neuronal_type_label <- function(cell.type) {
  if (cell.type %in% names(NEURONAL_TYPE_LABELS)) {
    NEURONAL_TYPE_LABELS[[cell.type]]
  } else {
    cell.type
  }
}

#' Load the topic-loading matrix (cells x programs) for a cell type from H5.
#'
#' @param cell.type Cell type key into CELL_TYPE_TO_GROUP (e.g. "inhibitory").
#' @param h5_path Path to the topics HDF5 file.
#'
#' @return Numeric matrix with cells as rows and only the k* program columns.
topic_loading_matrix_from_h5 <- function(cell.type, h5_path = TOPICS_DATA_H5_PATH) {
  h5_group <- CELL_TYPE_TO_GROUP[[cell.type]]
  L <- rhdf5::h5read(h5_path, paste0(h5_group, "/L")) %>%
    tibble::column_to_rownames("cell")
  as.matrix(L[, grep("^k", colnames(L)), drop = FALSE])
}

#' Load the topic-loading matrix from an in-memory fit (helper for tests / offline use).
#'
#' @param fit A fastTopics-style fit with element L (cells x programs).
#'
#' @return Numeric matrix with cells as rows and only the k* program columns.
topic_loading_matrix_from_fit <- function(fit) {
  L <- as.data.frame(fit[["L"]])
  as.matrix(L[, grep("^k", colnames(L)), drop = FALSE])
}

#' Resolve the increasing-program topic pair (stress vs cell cycle) for a cell type.
#'
#' @param cell.type Neuronal cell type ("inhibitory", "cux2p", "cux2m").
#'
#' @return List with elements stress (increasing_early) and cell_cycle (increasing_late).
get_increasing_topic_pair <- function(cell.type) {
  idx <- neuronal_cell_type_to_index(cell.type)
  list(
    stress = gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_early)[[idx]]),
    cell_cycle = gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_late)[[idx]])
  )
}

#' One within-cell permutation replicate for two topic slots.
#'
#' Null model: within each cell, randomly reassign which original topic column
#' contributes the loading to each slot (topic_a / topic_b), breaking topic-specific
#' coupling while keeping each cell's loading values on the simplex.
#'
#' Full row permutation would shuffle all p columns per cell; for correlating only
#' two slots, it suffices to sample two distinct source columns per cell:
#'   - src_a[i] ~ Uniform(1..p)             which column feeds slot topic_a in cell i
#'   - src_b[i] ~ Uniform(1..p \\ {src_a})   distinct column for slot topic_b
#' The pair (src_a, src_b) is uniform over all p*(p-1) ordered distinct pairs,
#' matching the marginal of a random within-row permutation at those two slots.
#'
#' replace=TRUE on sample.int draws indices independently across cells (not within a cell).
#' The src_b adjustment maps k in 1..(p-1) to 1..p while skipping src_a:
#'   if k < src_a then src_b = k; else src_b = k + 1
permute_within_cell_pair_cor <- function(L_mat, topic_a, topic_b, method = "spearman", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ia <- match(topic_a, colnames(L_mat))
  ib <- match(topic_b, colnames(L_mat))
  if (is.na(ia) || is.na(ib)) {
    stop(glue::glue("Topics not found in L: {topic_a}, {topic_b}"))
  }
  n <- nrow(L_mat)
  p <- ncol(L_mat)

  # Per cell: pick source column for slot topic_a (1..p, independent across cells)
  src_a <- sample.int(p, n, replace = TRUE)

  # Per cell: pick a different source column for slot topic_b
  src_b <- sample.int(p - 1L, n, replace = TRUE)
  src_b <- src_b + (src_b >= src_a)  # skip src_a: shift indices at/above the hole

  x <- L_mat[cbind(seq_len(n), src_a)]
  y <- L_mat[cbind(seq_len(n), src_b)]
  stats::cor(x, y, method = method)
}

#' Run the within-cell permutation null for one topic pair.
#'
#' Pure compute + summary; performs no file I/O.
#'
#' @param L_mat Numeric matrix (cells x programs); column names must include topic_a/topic_b.
#' @param topic_a,topic_b Program column names to correlate (e.g. "k10", "k8").
#' @param n_perm Number of permutations.
#' @param seed Base RNG seed; permutation b uses seed + b for reproducibility.
#' @param method Correlation method.
#' @param n_cores Cores for parallel permutations (sequential with progress bar otherwise).
#' @param cell.type Optional cell-type id stored in the output.
#' @param cell_type_label Optional human-readable label stored in the output.
#' @param pair_label Optional label for the pair (defaults to "topic_a-topic_b").
#'
#' @return List with null_df and summary_df.
run_topic_pair_permutation_null <- function(L_mat, topic_a, topic_b,
                                            n_perm = 200,
                                            seed = 42,
                                            method = "spearman",
                                            n_cores = max(1L, parallel::detectCores() - 1L),
                                            cell.type = NA_character_,
                                            cell_type_label = NA_character_,
                                            pair_label = NULL) {
  L_mat <- as.matrix(L_mat)
  ia <- match(topic_a, colnames(L_mat))
  ib <- match(topic_b, colnames(L_mat))
  if (is.na(ia) || is.na(ib)) {
    stop(glue::glue("Topics not found in L: {topic_a}, {topic_b}"))
  }
  n <- nrow(L_mat)
  if (is.null(pair_label)) pair_label <- paste(topic_a, topic_b, sep = "-")

  message(glue::glue(
    "Permutation null{if (is.na(cell.type)) '' else paste0(' for ', cell.type)}: ",
    "{topic_a} vs {topic_b} ({n_perm} permutations, seed = {seed})"
  ))

  r_obs <- stats::cor(L_mat[, ia], L_mat[, ib], method = method)
  perm_one <- function(b) {
    # seed + b: distinct, reproducible RNG state per permutation replicate
    permute_within_cell_pair_cor(L_mat, topic_a, topic_b, method = method, seed = seed + b)
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    message(glue::glue("Running {n_perm} permutations in parallel ({n_cores} cores)..."))
    r_perm <- unlist(parallel::mclapply(seq_len(n_perm), perm_one, mc.cores = min(n_cores, n_perm)))
  } else {
    message(glue::glue("Running {n_perm} permutations sequentially..."))
    r_perm <- numeric(n_perm)
    pb <- utils::txtProgressBar(min = 0, max = n_perm, style = 3)
    on.exit(close(pb), add = TRUE)
    for (b in seq_len(n_perm)) {
      r_perm[b] <- perm_one(b)
      utils::setTxtProgressBar(pb, b)
    }
  }

  p_one_sided <- (sum(r_perm <= r_obs) + 1L) / (length(r_perm) + 1L)

  null_df <- data.frame(
    cell_type = cell.type,
    perm_id = c(0L, seq_len(n_perm)),
    is_observed = c(TRUE, rep(FALSE, n_perm)),
    correlation = c(r_obs, r_perm),
    stringsAsFactors = FALSE
  )

  summary_df <- data.frame(
    cell_type = cell.type,
    cell_type_label = cell_type_label,
    topic_a = topic_a,
    topic_b = topic_b,
    pair_label = pair_label,
    correlation_observed = r_obs,
    correlation_null_mean = mean(r_perm),
    correlation_null_sd = stats::sd(r_perm),
    correlation_null_median = stats::median(r_perm),
    correlation_null_min = min(r_perm),
    correlation_null_max = max(r_perm),
    p_one_sided = p_one_sided,
    n_cells = n,
    n_perm = n_perm,
    cor_method = method,
    seed = seed,
    stringsAsFactors = FALSE
  )

  message(glue::glue(
    "r_obs = {round(r_obs, 4)}, null mean = {round(summary_df$correlation_null_mean, 4)}, ",
    "p (one-sided) = {signif(p_one_sided, 3)}"
  ))

  list(null_df = null_df, summary_df = summary_df)
}

#' Write permutation null result (null distribution + summary) to CSV.
#'
#' @param result List from run_topic_pair_permutation_null.
#' @param output_dir Directory to write into (created if needed).
#' @param file_prefix Prefix for the two CSV files.
#'
#' @return List with null_path and summary_path.
save_topic_pair_permutation_null <- function(result, output_dir, file_prefix) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  null_path <- file.path(output_dir, glue::glue("{file_prefix}_null_distribution.csv"))
  summary_path <- file.path(output_dir, glue::glue("{file_prefix}_summary.csv"))
  write.csv(result$null_df, null_path, row.names = FALSE)
  write.csv(result$summary_df, summary_path, row.names = FALSE)
  invisible(list(null_path = null_path, summary_path = summary_path))
}

#' Read a previously saved permutation null result back into a result list.
#'
#' @param output_dir Directory containing the CSVs.
#' @param file_prefix Prefix used when saving.
#'
#' @return List with null_df and summary_df, or NULL if files are missing.
read_topic_pair_permutation_null <- function(output_dir, file_prefix) {
  null_path <- file.path(output_dir, glue::glue("{file_prefix}_null_distribution.csv"))
  summary_path <- file.path(output_dir, glue::glue("{file_prefix}_summary.csv"))
  if (!file.exists(null_path) || !file.exists(summary_path)) {
    message(glue::glue("Missing {null_path} or {summary_path}"))
    return(NULL)
  }
  list(
    null_df = read.csv(null_path, stringsAsFactors = FALSE),
    summary_df = read.csv(summary_path, stringsAsFactors = FALSE)
  )
}

#' Plot the permutation null distribution with the observed correlation marked.
#'
#' The x-axis always spans zero so the observed value's sign is visible.
#'
#' @param result List with null_df and summary_df (e.g. from run_/read_ functions).
#' @param output_path Output PDF path.
#' @param title Optional plot title (defaults to a cell-type-aware title).
#'
#' @return Invisibly, the output path.
plot_topic_pair_permutation_null <- function(result, output_path, title = NULL) {
  null_df <- result$null_df
  summary_df <- result$summary_df
  r_obs <- summary_df$correlation_observed[[1]]
  pair_label <- summary_df$pair_label[[1]]
  p_one_sided <- summary_df$p_one_sided[[1]]
  cell_label <- summary_df$cell_type_label[[1]]
  if (is.null(title)) {
    title <- glue::glue("Within-cell permutation null: {cell_label}")
  }

  perm_cor <- null_df$correlation[!null_df$is_observed]
  x_min <- min(c(perm_cor, r_obs, 0))
  x_max <- max(c(perm_cor, r_obs, 0))
  x_pad <- 0.02 * (x_max - x_min)

  p <- ggplot(null_df %>% dplyr::filter(!is_observed), aes(x = correlation)) +
    geom_histogram(bins = 30, fill = "grey80", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = r_obs, color = "#c7522a", linewidth = 0.9) +
    coord_cartesian(xlim = c(x_min - x_pad, x_max + x_pad)) +
    annotate(
      "text", x = r_obs, y = Inf, hjust = -0.05, vjust = 1.5,
      label = glue::glue("{pair_label}: r = {round(r_obs, 3)}, p = {signif(p_one_sided, 3)}"),
      color = "#c7522a", size = 3.5
    ) +
    labs(
      title = title,
      subtitle = glue::glue(
        "{summary_df$n_perm[[1]]} permutations (seed = {summary_df$seed[[1]]}); ",
        "one-sided p = {signif(p_one_sided, 3)}"
      ),
      x = glue::glue("{summary_df$cor_method[[1]]} correlation ({pair_label})"),
      y = "Count"
    ) +
    theme_classic(base_size = 11)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  withr::with_pdf(output_path, width = 6.5, height = 4.5, print(p))
  message(glue::glue("Permutation null plot saved to {output_path}"))
  invisible(output_path)
}

#' End-to-end increasing-programs permutation null for a cell type (H5-first).
#'
#' Loads L from H5, resolves the stress (increasing_early) vs cell-cycle
#' (increasing_late) topic pair, runs the permutation null, and writes CSVs + PDF.
#' Self-contained: requires only that L exists in H5 for this cell type.
#'
#' @param cell.type Neuronal cell type ("inhibitory", "cux2p", "cux2m").
#' @param output_dir Directory for CSV + PDF outputs.
#' @param h5_path Path to the topics HDF5 file.
#' @param file_prefix Prefix for output files (defaults to cell.type).
#' @param n_perm,seed,method,n_cores Passed to run_topic_pair_permutation_null.
#'
#' @return List with null_df and summary_df (schema uses topic_stress/topic_cell_cycle).
run_increasing_programs_permutation_null <- function(cell.type, output_dir,
                                                     h5_path = TOPICS_DATA_H5_PATH,
                                                     file_prefix = cell.type,
                                                     n_perm = 200,
                                                     seed = 42,
                                                     method = "spearman",
                                                     n_cores = max(1L, parallel::detectCores() - 1L)) {
  pair <- get_increasing_topic_pair(cell.type)
  L_mat <- topic_loading_matrix_from_h5(cell.type, h5_path = h5_path)

  result <- run_topic_pair_permutation_null(
    L_mat,
    topic_a = pair$stress,
    topic_b = pair$cell_cycle,
    n_perm = n_perm,
    seed = seed,
    method = method,
    n_cores = n_cores,
    cell.type = cell.type,
    cell_type_label = neuronal_type_label(cell.type),
    pair_label = paste(pair$stress, pair$cell_cycle, sep = "-")
  )

  # Backward-compatible summary schema: stress = topic_a, cell cycle = topic_b.
  result$summary_df <- result$summary_df %>%
    dplyr::rename(topic_stress = topic_a, topic_cell_cycle = topic_b)

  save_topic_pair_permutation_null(result, output_dir, file_prefix)
  plot_topic_pair_permutation_null(
    result,
    file.path(output_dir, glue::glue("{file_prefix}_null_distribution.pdf"))
  )
  result
}
