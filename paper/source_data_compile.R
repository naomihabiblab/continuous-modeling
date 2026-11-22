# Build the Nature Source Data files from the staging tables written by the
# sd_* helpers (R/source_data.R).
#
#   {SD_DIR}/staging/<Fig>/<panel>__<part>.csv.gz   (written by figures.R)
#   -> {SD_DIR}/final/Source_Data_<Fig>.xlsx        (one sheet per table + README)
#   -> {SD_DIR}/final/Source_Data_<Fig>.zip         (only for tables too big for the xlsx)
#
# Reads only the staging folder, so it is cheap and can run locally against the
# read-only mount:
#   SD_DIR=<result_dir>/source_data OUT_DIR=<local dir> Rscript paper/source_data_compile.R
# OUT_DIR defaults to {SD_DIR}/final.

suppressPackageStartupMessages({
  library(openxlsx)
  library(data.table)
})

source("R/utils.R")
SD_DIR  <- Sys.getenv("SD_DIR", file.path(get_result_dir(), "source_data"))
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(SD_DIR, "final"))
STAGING <- file.path(SD_DIR, "staging")
MAX_XLSX_ROWS <- as.integer(Sys.getenv("MAX_XLSX_ROWS", "100000"))   # bigger tables go to the zip
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ROSMAP participants are shown by their ROSMAP individualID (as in the Supplementary Tables,
# paper/create.supplementary.tables.R), not by the internal projid.
ROSMAP_CLINICAL <- file.path(get_shared_dir(), "NextSeq/500/v1.1.objects.for.synapse/ROSMAP_clinical.csv")
rosmap <- fread(ROSMAP_CLINICAL, select = c("projid", "individualID"))
projid_to_individualID <- setNames(rosmap$individualID, as.character(as.integer(rosmap$projid)))
as_projid <- function(v) suppressWarnings(as.character(as.integer(as.character(v))))

# Replace ROSMAP projid columns with individualID, and drop the leftover integer row index
# ("row" = 0,1,2,... from pandas/rownames) that carries no information.
clean_ids <- function(dt, key, log) {
  for (cn in names(dt)) {
    v <- dt[[cn]]
    ok <- !is.na(v) & as.character(v) != ""
    if (!any(ok)) next
    is_projid <- !is.na(projid_to_individualID[as_projid(v[ok])])
    if (mean(is_projid) > 0.9 && (cn %in% c("projid", "row", "individual", "donor", "donor_id") || mean(is_projid) == 1)) {
      if (!all(is_projid)) stop(key, ": column ", cn, " has projids without a ROSMAP individualID")
      new <- unname(projid_to_individualID[as_projid(v)]); new[!ok] <- NA
      dt[[cn]] <- new
      setnames(dt, cn, "individualID")
      log$ids <- c(log$ids, paste0(key, ":", cn))
    } else if (cn == "row" && is.numeric(v) && all(v[ok] == round(v[ok]))) {
      dt[[cn]] <- NULL
      log$rows <- c(log$rows, key)
    }
  }
  for (cn in names(dt)) {
    if ((tolower(cn) %in% c("cell", "barcode", "cell_id", "row")) && is_cell_id(dt[[cn]])) {
      dt[[cn]] <- seq_len(nrow(dt))
      setnames(dt, cn, "cell")
      log$cells <- c(log$cells, paste0(key, ":", cn))
    }
  }
  # per-cell tables do not need the donor id (no test is run per donor on them)
  if ("cell" %in% names(dt) && "individualID" %in% names(dt)) {
    dt[["individualID"]] <- NULL
    log$donor_dropped <- c(log$donor_dropped, key)
  }
  # matrices with participants as columns
  pc <- as_projid(names(dt)); hit <- !is.na(projid_to_individualID[pc])
  if (sum(hit) > 1) { setnames(dt, names(dt)[hit], unname(projid_to_individualID[pc[hit]])); log$ids <- c(log$ids, paste0(key, ":<columns>")) }
  dt
}
id_log <- new.env(); id_log$ids <- character(); id_log$rows <- character(); id_log$cells <- character(); id_log$donor_dropped <- character()

# Structure-plot tables (Fig 1b) are written long (sample, topic, prop); show them as one row per
# plotted cell (column order of the plot) with one column per program.
widen_structure <- function(dt) {
  if (!setequal(names(dt), c("sample", "topic", "prop"))) return(NULL)
  w <- dcast(dt, sample ~ factor(topic, levels = unique(topic)), value.var = "prop")
  setnames(w, "sample", "cell")
  w
}

# Cell barcodes (snRNA-seq, e.g. XXXXX-XX-A_AACCTGACACGTCATA-1) and spatial cell ids are replaced by the
# row number within the table (1..n); the full ids are in the deposited raw data.
is_cell_id <- function(v) {
  v <- as.character(v[!is.na(v)])
  length(v) > 0 && mean(grepl("[ACGT]{10,}", v) | grepl("^[0-9]{12,}$", v)) > 0.9
}

fig_label <- function(fig) {
  if (grepl("^ED", fig)) paste("Extended Data Fig.", sub("^ED", "", fig)) else paste("Fig.", sub("^Fig", "", fig))
}
out_name <- function(fig) {
  if (grepl("^ED", fig)) paste0("Source_Data_ED_Fig", sub("^ED", "", fig)) else paste0("Source_Data_Fig", sub("^Fig", "", fig))
}

# Excel sheet names: <= 31 chars, no []:*?/\ , unique within the workbook
sheet_name <- function(key, used) {
  nm <- gsub("[]:*?/\\\\[]", "_", key)
  nm <- gsub("__", " ", nm)
  if (nchar(nm) > 31) nm <- substr(nm, 1, 31)
  base <- nm; i <- 1
  while (tolower(nm) %in% tolower(used)) { i <- i + 1; nm <- paste0(substr(base, 1, 28), "~", i) }
  nm
}

figs <- list.dirs(STAGING, full.names = FALSE, recursive = FALSE)
figs <- figs[order(grepl("^ED", figs), as.integer(gsub("\\D", "", figs)))]

for (fig in figs) {
  files <- sort(list.files(file.path(STAGING, fig), pattern = "\\.csv\\.gz$", full.names = TRUE))
  if (length(files) == 0) next
  keys <- sub("\\.csv\\.gz$", "", basename(files))
  notes <- if (file.exists(file.path(STAGING, fig, "_README.csv")))
    read.csv(file.path(STAGING, fig, "_README.csv"), stringsAsFactors = FALSE) else
    data.frame(key = character(), panel = character(), description = character())

  wb <- createWorkbook()
  readme_rows <- list()
  used <- "README"
  zip_files <- character()
  zip_dir <- file.path(tempdir(), out_name(fig))
  unlink(zip_dir, recursive = TRUE); dir.create(zip_dir)

  for (k in seq_along(files)) {
    key <- keys[k]
    panel <- sub("__.*$", "", key)
    dt <- clean_ids(fread(files[k]), file.path(fig, key), id_log)
    wide <- widen_structure(dt)
    if (!is.null(wide)) dt <- wide
    # description: exact key, else the longest README key that prefixes this key (e.g. part tables).
    # Statistics tables (key ends in _stats) without their own note say so, instead of inheriting
    # the description of the plotted data.
    desc <- notes$description[notes$key == key]
    if (length(desc) == 0) {
      parent <- sub("_+stats$", "", key)
      cand <- notes[startsWith(parent, notes$key), , drop = FALSE]
      desc <- if (nrow(cand)) cand$description[which.max(nchar(cand$key))] else ""
      if (grepl("_stats$", key))
        desc <- paste0("Statistics behind the stars/P values shown in this panel (test, P value and adjusted P value per comparison)",
                       if (nzchar(desc)) paste0(", for: ", desc) else "")
    }
    if (!is.null(wide)) desc <- sub("\\(long format\\)", "(one row per cell, in the plotted order; one column per program)", desc)
    if (nrow(dt) > MAX_XLSX_ROWS) {
      fname <- paste0(key, ".csv.gz")
      fwrite(dt, file.path(zip_dir, fname))   # the cleaned table (individualID, no row index)
      zip_files <- c(zip_files, fname)
      location <- paste0(out_name(fig), ".zip / ", fname)
    } else {
      sh <- sheet_name(key, used); used <- c(used, sh)
      addWorksheet(wb, sh)
      writeData(wb, sh, dt)
      location <- paste0("sheet '", sh, "'")
    }
    readme_rows[[k]] <- data.frame(Panel = paste0(fig_label(fig), panel), Table = key, Location = location,
                                   Rows = nrow(dt), Columns = ncol(dt), Description = desc)
  }

  readme <- do.call(rbind, readme_rows)
  addWorksheet(wb, "README")
  writeData(wb, "README", data.frame(x = c(paste("Source Data for", fig_label(fig)),
                                           "One sheet per plotted table. Heatmap matrices are in the displayed row/column order.",
                                           "Dynamics panels: '_individuals' = per-individual values; '_fitted_curve' = fit and se.fit (ribbon = fit +/- 2 se.fit).",
                                           "Tables too large for Excel are provided as csv.gz in the accompanying zip file.", "")),
            colNames = FALSE)
  writeData(wb, "README", readme, startRow = 7)
  worksheetOrder(wb) <- c(length(names(wb)), seq_len(length(names(wb)) - 1))
  saveWorkbook(wb, file.path(OUT_DIR, paste0(out_name(fig), ".xlsx")), overwrite = TRUE)

  if (length(zip_files)) {
    zip::zip(file.path(OUT_DIR, paste0(out_name(fig), ".zip")), files = zip_files, root = zip_dir)
  }
  sz <- file.size(file.path(OUT_DIR, paste0(out_name(fig), c(".xlsx", ".zip"))))
  cat(sprintf("%-6s %3d tables -> xlsx %.1f MB%s\n", fig, length(files), sz[1] / 1e6,
              if (length(zip_files)) sprintf(", zip %.1f MB (%d tables)", sz[2] / 1e6, length(zip_files)) else ""))
  if (any(sz > 30e6, na.rm = TRUE)) cat("   WARNING: file over 30 MB; consider a repository for the large tables\n")
}

cat("\nROSMAP projid -> individualID:", length(id_log$ids), "columns\n  ", paste(id_log$ids, collapse = "\n   "), "\n")
cat("Dropped integer row index:", length(id_log$rows), "tables\n  ", paste(id_log$rows, collapse = "\n   "), "\n")
cat("Donor id dropped from per-cell tables:", length(id_log$donor_dropped), "tables\n  ", paste(id_log$donor_dropped, collapse = "\n   "), "\n")
cat("Cell ids replaced by numbers:", length(id_log$cells), "columns\n  ", paste(id_log$cells, collapse = "\n   "), "\n")
