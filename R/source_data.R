# Source Data helpers for paper/figures.R
#
# Every sd_* call writes one gzipped CSV right away to
#   {SD_DIR}/staging/<fig>/<panel>[__<part>].csv.gz
# and overwrites any earlier file with the same key. Loops, re-runs and running
# blocks by hand therefore need no special handling. paper/source_data_compile.R
# later turns the staging folders into one xlsx (+ zip) per figure.
#
# The helpers never change what they are given and never draw random numbers,
# so adding them to a plotting block cannot change the plot.
#
# Switches (environment variables, read when this file is sourced):
#   SD_ENABLED=FALSE   turn every sd_* call into a no-op
#   SD_DIR=<path>      output root (default: {get_result_dir()}/source_data)
#
# A fig/panel of NA or NULL means "not in the paper" and the call does nothing,
# so a lookup table can route loop iterations (see sd_lookup()).

SD_ENABLED <- toupper(Sys.getenv("SD_ENABLED", "TRUE")) != "FALSE"
SD_DIR     <- Sys.getenv("SD_DIR", file.path(get_result_dir(), "source_data"))
SD_STAGING <- file.path(SD_DIR, "staging")

.sd_skip <- function(fig, panel) {
  !SD_ENABLED || is.null(fig) || is.null(panel) || length(fig) == 0 ||
    is.na(fig) || is.na(panel)
}

.sd_key <- function(panel, part) {
  part <- part[!is.na(part) & nzchar(part)]
  if (length(part) == 0) return(panel)
  paste0(panel, "__", gsub("[^A-Za-z0-9.+-]+", "_", paste(part, collapse = "_")))
}

.sd_as_df <- function(data) {
  if (is.matrix(data)) {
    rn <- rownames(data)
    df <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(rn)) df <- cbind(row = rn, df, stringsAsFactors = FALSE)
    return(df)
  }
  df <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  # keep meaningful rownames (not 1..n)
  rn <- rownames(df)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df)))) && !"row" %in% colnames(df)) {
    df <- cbind(row = rn, df, stringsAsFactors = FALSE)
  }
  # list columns (e.g. from ggplot layer data) cannot go to CSV
  is_list <- vapply(df, is.list, logical(1))
  df[is_list] <- lapply(df[is_list], function(x) vapply(x, function(v) paste(format(v), collapse = ";"), ""))
  df
}

#' Write one table now.
#'
#' @param data data.frame, tibble or matrix (matrix rownames go to column "row")
#' @param fig  figure id used as folder and later as file name, e.g. "Fig1", "ED2"
#' @param panel panel letter(s), e.g. "c" or "a-d"
#' @param part optional sub-key (loop iteration, cell type, table name); several are joined by "_"
#' @param description optional one-line text for the README sheet
#' @return data, invisibly
sd_add <- function(data, fig, panel, part = NULL, description = NULL) {
  if (.sd_skip(fig, panel)) return(invisible(data))
  dir <- file.path(SD_STAGING, fig)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  key <- .sd_key(panel, part)
  df <- .sd_as_df(data)
  data.table::fwrite(df, file.path(dir, paste0(key, ".csv.gz")))
  message(sprintf("[source data] %s/%s: %d rows x %d cols", fig, key, nrow(df), ncol(df)))
  if (!is.null(description)) sd_note(fig, panel, description, part)
  invisible(data)
}

#' Pass-through: save the data behind a plot object and return the object unchanged.
#'
#' Works for ggplot (plot data + layers with their own data, ggraph edges),
#' drawn ComplexHeatmap HeatmapList (matrices in displayed order), undrawn
#' Heatmap/HeatmapList without clustering, data.frame and matrix.
sd_tee <- function(x, fig, panel, part = NULL, description = NULL) {
  if (.sd_skip(fig, panel)) return(x)
  tabs <- sd_extract(x)
  if (length(tabs) == 0) warning("sd_tee: nothing to save for ", fig, " ", panel)
  multi <- length(tabs) > 1
  for (i in seq_along(tabs)) {
    nm <- names(tabs)[i]
    sub <- if (multi && !is.null(nm) && nzchar(nm)) nm else NULL
    sd_add(tabs[[i]], fig, panel, c(part, sub))
  }
  if (!is.null(description)) sd_note(fig, panel, description, part)
  x
}

sd_extract <- function(x) UseMethod("sd_extract")
sd_extract.default <- function(x) {
  warning("sd_extract: no method for class ", paste(class(x), collapse = "/"))
  list()
}
sd_extract.data.frame <- function(x) list(x)
sd_extract.matrix <- function(x) list(x)

sd_extract.gg <- function(x) {
  out <- list()
  if (is.data.frame(x$data) && nrow(x$data) > 0) out[["data"]] <- x$data
  for (i in seq_along(x$layers)) {
    d <- x$layers[[i]]$data
    if (is.data.frame(d) && nrow(d) > 0) out[[paste0("layer", i)]] <- d
  }
  g <- attr(x$data, "graph")              # ggraph: add the edge list
  if (!is.null(g)) out[["edges"]] <- igraph::as_data_frame(g, what = "edges")
  out
}

# ComplexHeatmap -----------------------------------------------------------

.sd_heatmaps <- function(x) {
  if (methods::is(x, "Heatmap")) return(list(x))
  Filter(function(h) methods::is(h, "Heatmap"), x@ht_list)
}

# One heatmap's matrix in displayed order, with the row slice (row_split) as a column.
.sd_ordered_matrix <- function(h, ro, co) {
  m <- h@matrix
  if (is.null(rownames(m))) rownames(m) <- seq_len(nrow(m))
  ro <- if (is.list(ro)) ro else list(ro)
  co <- if (is.list(co)) co else list(co)
  r <- unlist(ro, use.names = FALSE)
  cc <- unlist(co, use.names = FALSE)
  if (is.null(cc) || length(cc) == 0) cc <- seq_len(ncol(m))
  df <- .sd_as_df(m[r, cc, drop = FALSE])
  # simple row annotations (e.g. LFC, group) as extra columns, in displayed row order
  for (side in c("left_annotation", "right_annotation")) {
    ha <- tryCatch(methods::slot(h, side), error = function(e) NULL)
    if (is.null(ha)) next
    for (nm in names(ha@anno_list)) {
      v <- tryCatch({
        val <- ha@anno_list[[nm]]@fun@var_env$value   # set for simple annotations; NULL for anno_mark/anno_text
        if (is.matrix(val)) NULL else val
      }, error = function(e) NULL)
      if (!is.null(v) && length(v) == nrow(m)) df[[paste0("annotation_", nm)]] <- v[r]
    }
  }
  if (length(ro) > 1) {
    slice <- rep(if (is.null(names(ro))) seq_along(ro) else names(ro), lengths(ro))
    df <- cbind(row_slice = slice, df, stringsAsFactors = FALSE)
  }
  df
}

sd_extract.HeatmapList <- function(x) {
  hts <- .sd_heatmaps(x)
  drawn <- isTRUE(x@layout$initialized)
  if (!drawn) return(.sd_undrawn(hts))
  ro <- ComplexHeatmap::row_order(x)
  co <- ComplexHeatmap::column_order(x)
  out <- lapply(hts, function(h) {
    co_h <- if (length(hts) > 1) co[[h@name]] else co
    .sd_ordered_matrix(h, ro, co_h)
  })
  names(out) <- vapply(hts, function(h) h@name, "")
  out
}

# Undrawn: the order is only known when nothing is clustered or reordered.
sd_extract.Heatmap <- function(x) .sd_undrawn(list(x))

.sd_undrawn <- function(hts) {
  out <- lapply(hts, function(h) {
    if (isTRUE(h@row_dend_param$cluster) || isTRUE(h@column_dend_param$cluster)) {
      warning("sd_tee on an undrawn, clustered heatmap: matrix saved in input order, ",
              "not displayed order. Call sd_tee after draw() instead.")
    }
    ro <- h@row_order; co <- h@column_order
    .sd_ordered_matrix(h, if (length(ro)) ro else seq_len(nrow(h@matrix)),
                       if (length(co)) co else seq_len(ncol(h@matrix)))
  })
  names(out) <- vapply(hts, function(h) h@name, "")
  out
}

# Statistics ----------------------------------------------------------------

#' correlate_dataframes() result -> long table (row, column, correlation, adjusted p).
sd_corr <- function(c, fig, panel, part = NULL, description = NULL) {
  if (.sd_skip(fig, panel)) return(invisible(c))
  long <- function(m, value) {
    m <- as.matrix(m)
    df <- data.frame(row = rep(rownames(m), times = ncol(m)),
                     column = rep(colnames(m), each = nrow(m)),
                     value = as.vector(m), stringsAsFactors = FALSE)
    stats::setNames(df, c("row", "column", value))
  }
  df <- long(c$corr, "correlation")
  if (!is.null(c$adj.pval)) df$adj_pval <- as.vector(as.matrix(c$adj.pval))
  if (!is.null(c$pval)) df$pval <- as.vector(as.matrix(c$pval))
  if (!is.null(c$n)) df$n <- as.vector(as.matrix(c$n))
  sd_add(df, fig, panel, c(part, "stats"), description)
  invisible(c)
}

#' Dynamics panels: the two tables plot.dynamics() draws from.
#'
#' Same feature/trajectory filter as plot.dynamics.wrapper().
#' _individuals : x = pseudotime, y = observed per-individual value, X.weights. = trajectory probability
#' _fitted_curve: x, fit, se.fit (the ribbon is fit +/- 2 * se.fit)
sd_dynamics <- function(dynamics, features, fig, panel, part = NULL, trajectories = NULL,
                        labels = NULL, description = NULL) {
  if (.sd_skip(fig, panel)) return(invisible(NULL))
  as_r <- function(d) if (inherits(d, "python.builtin.object")) reticulate::py_to_r(d) else d
  fv <- as_r(dynamics$fitted.vals)
  pv <- as_r(dynamics$pred.vals)
  if (is.null(trajectories)) trajectories <- unique(fv$trajectory)
  keep <- function(d) {
    d <- d[d$feature %in% features & d$trajectory %in% trajectories, , drop = FALSE]
    if (!is.null(labels)) d$label <- unname(labels[as.character(d$feature)])
    d
  }
  sd_add(keep(fv), fig, panel, c(part, "individuals"))
  sd_add(keep(pv), fig, panel, c(part, "fitted_curve"))
  if (!is.null(description)) sd_note(fig, panel, description, part)
  invisible(NULL)
}

# Helpers ------------------------------------------------------------------

#' Route loop iterations to figures: returns c(fig, panel) or c(NA, NA).
#' lookup <- list(astrocytes = c("Fig1", "c"), microglia = c("ED2", "b"))
sd_lookup <- function(lookup, key) {
  v <- lookup[[key]]
  if (is.null(v)) c(NA_character_, NA_character_) else v
}

#' One README line per panel/part; collected into the README sheet.
sd_note <- function(fig, panel, text, part = NULL) {
  if (.sd_skip(fig, panel)) return(invisible(NULL))
  dir <- file.path(SD_STAGING, fig)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, "_README.csv")
  key <- .sd_key(panel, part)
  old <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else
    data.frame(key = character(), panel = character(), description = character())
  old <- old[old$key != key, , drop = FALSE]
  new <- rbind(old, data.frame(key = key, panel = panel, description = text))
  utils::write.csv(new[order(new$key), ], f, row.names = FALSE)
  invisible(NULL)
}

#' Trait-association table behind plot.trait.associations.heatmap(): the long
#' table (covariate/state, trait, t statistic, p, adjusted p) for the plotted traits.
sd_trait <- function(df, params, fig, panel, part = NULL, description = NULL) {
  if (.sd_skip(fig, panel)) return(invisible(df))
  d <- if (inherits(df, "python.builtin.object")) reticulate::py_to_r(df) else df
  d <- as.data.frame(d, check.names = FALSE)
  if (!is.null(params) && "trait" %in% colnames(d)) d <- d[d$trait %in% params, , drop = FALSE]
  sd_add(d, fig, panel, c(part, "stats"), description)
  invisible(df)
}
