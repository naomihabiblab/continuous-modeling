# Run selected figure blocks of paper/figures.R without running the whole file.
#
# figures.R is a list of top-level blocks of the form (function() { ... })().
# This script parses the file (nothing is run by parsing) and then:
#   1. runs every top-level statement that is not a figure block (library(),
#      source(), colours, constants), in file order. Statements that read an
#      h5ad/rds into a variable (e.g. `data <- anndata::read_h5ad(...)`) are
#      turned into delayedAssign(), so the file is only read if a selected
#      block actually uses that variable;
#   2. runs only the selected blocks, each inside tryCatch, logging time and errors.
#
# Plots go to FIGURES_DIR (default {SD_DIR}/figures_rerun), never to the
# published figures directory. Source data go to {SD_DIR}/staging.
#
# Usage (from the project root, on a compute node / via submit_job.py):
#   Rscript paper/run_figure_blocks.R --list
#   Rscript paper/run_figure_blocks.R --blocks "Figure 1c,Figure 1d"
#   Rscript paper/run_figure_blocks.R --blocks 123,192
#
# --blocks takes comma-separated tokens. A number selects the block that
# contains that line of figures.R; anything else is a case-insensitive regex
# matched against the section comment above the block (e.g. "Figure 1c",
# "Extended Data Figure 2c"). Blocks run in file order.
# Blocks that figures.R defines but never calls (`function() {...}` without
# the trailing `()`) are called when selected; the log says so.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
# Every non-flag argument must be the value of a flag. A stray one usually means
# a block name with a space was split by a wrapper (submit_job.py passes
# `--blocks "Figure 4f"` as `--blocks Figure 4f`, which would select every block).
.vals <- which(!startsWith(args, "--"))
.stray <- .vals[.vals == 1 | !startsWith(args[pmax(.vals - 1, 1)], "--")]
if (length(.stray)) stop("Unexpected argument(s): ", paste(args[.stray], collapse = " "),
                         ". Use line numbers or a regex without spaces (e.g. 'Figure.4f').")
FIG_FILE  <- get_arg("--file", "paper/figures.R")
LIST_ONLY <- isTRUE(get_arg("--list", FALSE))
BLOCKS    <- get_arg("--blocks", NULL)

# Parse ---------------------------------------------------------------------

exprs <- parse(FIG_FILE, keep.source = TRUE)
srcrefs <- attr(exprs, "srcref")
src_lines <- readLines(FIG_FILE, warn = FALSE)
start_line <- vapply(srcrefs, function(s) s[1], 1L)
end_line   <- vapply(srcrefs, function(s) s[3], 1L)

is_fun <- function(e) is.call(e) && identical(e[[1]], as.name("function"))
kind <- vapply(exprs, function(e) {
  if (is_fun(e)) return("uncalled")                                        # function() {...}
  if (is.call(e) && identical(e[[1]], as.name("(")) && is_fun(e[[2]])) return("uncalled")  # (function() {...})
  if (is.call(e) && is.call(e[[1]]) && identical(e[[1]][[1]], as.name("(")) &&
      is_fun(e[[1]][[2]])) return("block")                                  # (function() {...})()
  "setup"
}, "")

# Section header of a block: the last comment line mentioning a figure, between
# the previous top-level statement and this one.
header_of <- function(i) {
  from <- if (i > 1) end_line[i - 1] + 1 else 1
  to <- start_line[i] - 1
  if (to < from) return("")
  cand <- grep("^\\s*#.*(Fig|Figure|ED)", src_lines[from:to], value = TRUE)
  cand <- trimws(gsub("[#+-]+\\s*$", "", gsub("^[#+ ]+", "", cand)))
  cand <- cand[nzchar(cand)]
  if (length(cand) == 0) "" else cand[length(cand)]
}
blocks <- which(kind != "setup")
block_info <- data.frame(
  index = blocks, start = start_line[blocks], end = end_line[blocks],
  kind = kind[blocks], header = vapply(blocks, header_of, ""),
  stringsAsFactors = FALSE)

if (LIST_ONLY) {
  for (r in seq_len(nrow(block_info))) with(block_info[r, ], cat(sprintf(
    "%5d-%-5d %-8s %s\n", start, end, if (kind == "uncalled") "UNCALLED" else "", header)))
  quit(save = "no")
}
if (is.null(BLOCKS) || isTRUE(BLOCKS)) stop("Pass --blocks (or --list)")

# Select ----------------------------------------------------------------------

tokens <- trimws(strsplit(BLOCKS, ",")[[1]])
sel <- integer()
for (tok in tokens) {
  if (grepl("^[0-9]+$", tok)) {
    ln <- as.integer(tok)
    hit <- block_info$index[block_info$start <= ln & block_info$end >= ln]
  } else {
    hit <- block_info$index[grepl(tok, block_info$header, ignore.case = TRUE)]
  }
  if (length(hit) == 0) stop("No block matches '", tok, "'. Use --list to see blocks.")
  sel <- c(sel, hit)
}
sel <- sort(unique(sel))
cat("Selected blocks:\n")
for (.blk in sel) cat(sprintf("  lines %d-%d  %s%s\n", start_line[.blk], end_line[.blk],
                           block_info$header[block_info$index == .blk],
                           if (kind[.blk] == "uncalled") "  [never called in figures.R; calling it here]" else ""))

# Output directories -------------------------------------------------------------

source("R/utils.R")   # get_result_dir()
if (Sys.getenv("SD_DIR") == "") Sys.setenv(SD_DIR = file.path(get_result_dir(), "source_data"))
if (Sys.getenv("FIGURES_DIR") == "") Sys.setenv(FIGURES_DIR = file.path(Sys.getenv("SD_DIR"), "figures_rerun"))
dir.create(Sys.getenv("FIGURES_DIR"), recursive = TRUE, showWarnings = FALSE)
cat("FIGURES_DIR:", Sys.getenv("FIGURES_DIR"), "\nSD_DIR:     ", Sys.getenv("SD_DIR"), "\n")
published_dir <- normalizePath(file.path(get_result_dir(), "figures"), mustWork = FALSE)

# Setup statements --------------------------------------------------------------

# `x <- <read call>` at top level -> delayedAssign
lazy_readers <- c("read_h5ad", "readRDS", "readr::read_rds")
lazy_target <- function(e) {
  if (!(is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]))) return(NULL)
  fn <- e[[3]]
  if (!is.call(fn)) return(NULL)
  name <- paste(deparse(fn[[1]]), collapse = "")
  if (sub("^.*::", "", name) %in% sub("^.*::", "", lazy_readers)) as.character(e[[2]]) else NULL
}

t0 <- Sys.time()
for (.stp in which(kind == "setup")) {
  e <- exprs[[.stp]]
  var <- lazy_target(e)
  if (!is.null(var)) {
    cat(sprintf("[setup] line %d: `%s` loads on first use\n", start_line[.stp], var))
    do.call(delayedAssign, list(var, e[[3]], eval.env = globalenv(), assign.env = globalenv()))
    next
  }
  ts <- Sys.time()
  eval(e, envir = globalenv())
  dt <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
  if (dt > 5) cat(sprintf("[setup] line %d took %.0fs: %s\n", start_line[.stp], dt,
                          substr(paste(deparse(e)[1], collapse = ""), 1, 80)))
}
cat(sprintf("[setup] done in %.0fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Packages loaded by the sourced files (AnnotationDbi, S4Vectors, ...) mask
# dplyr::select/rename/filter etc. The figures were made in an interactive
# session where dplyr won; re-attach dplyr so it is first on the search path.
for (pkg in c("tidyr", "dplyr")) {
  suppressWarnings(detach(paste0("package:", pkg), character.only = TRUE, unload = FALSE))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("[setup] re-attached dplyr/tidyr; select() comes from:", environmentName(environment(select)), "\n")

if (!exists("figures_dir") ||
    normalizePath(figures_dir, mustWork = FALSE) == published_dir) {
  stop("figures_dir points to the published figures directory; refusing to run.")
}
if (!exists("sd_add")) stop("R/source_data.R was not sourced by figures.R")
cat("figures_dir:", figures_dir, "\n")

# Blocks -------------------------------------------------------------------
# Loop variables are named .stp/.blk (not i) so a block that uses an undefined
# global `i` fails instead of silently picking up the runner's value.

results <- data.frame(lines = character(), header = character(), status = character(),
                      seconds = numeric(), message = character(), stringsAsFactors = FALSE)
for (.blk in sel) {
  hdr <- block_info$header[block_info$index == .blk]
  cat(sprintf("\n==== lines %d-%d: %s ====\n", start_line[.blk], end_line[.blk], hdr))
  e <- exprs[[.blk]]
  if (kind[.blk] == "uncalled") e <- as.call(list(e))
  grDevices::pdf(NULL)   # default device: stray plots outside with_pdf() go nowhere (no Rplots.pdf)
  ts <- Sys.time()
  err <- NULL
  n_warn <- 0
  tryCatch(
    withCallingHandlers(
      eval(e, envir = globalenv()),
      warning = function(w) {
        n_warn <<- n_warn + 1
        message("warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      error = function(cnd) {
        calls <- sys.calls()
        cat("Traceback (innermost last):\n")
        for (cl in utils::tail(calls, 15)) cat("  ", substr(paste(deparse(cl)[1], collapse = ""), 1, 120), "\n")
      }),
    error = function(cnd) err <<- conditionMessage(cnd))
  dt <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
  while (grDevices::dev.cur() > 1) grDevices::dev.off()   # close a device left open by an error
  status <- if (is.null(err)) "OK" else "ERROR"
  cat(sprintf("==== %s in %.0fs (%d warnings)%s\n", status, dt, n_warn,
              if (is.null(err)) "" else paste0(": ", err)))
  results <- rbind(results, data.frame(lines = sprintf("%d-%d", start_line[.blk], end_line[.blk]),
                                       header = hdr, status = status, seconds = round(dt),
                                       message = if (is.null(err)) "" else err))
}

cat("\n==== Summary ====\n")
print(results, right = FALSE, row.names = FALSE)
cat(sprintf("Total %.0f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (any(results$status != "OK")) quit(save = "no", status = 1)
