# Run spatial/spatial_figures.R for the Source Data, writing plots to
# {SD_DIR}/figures_rerun/spatial instead of the published Plots/spatial.
#   submit_job.py ... paper/run_spatial_source_data.R
source("R/utils.R")
if (Sys.getenv("SD_DIR") == "") Sys.setenv(SD_DIR = file.path(get_result_dir(), "source_data"))
Sys.setenv(SPATIAL_PLOTS_DIR = file.path(Sys.getenv("SD_DIR"), "figures_rerun", "spatial"))
stopifnot(normalizePath(Sys.getenv("SPATIAL_PLOTS_DIR"), mustWork = FALSE) !=
          normalizePath(file.path(get_result_dir(), "spatial"), mustWork = FALSE))
cat("SPATIAL_PLOTS_DIR:", Sys.getenv("SPATIAL_PLOTS_DIR"), "\n")
grDevices::pdf(NULL)   # stray plots go nowhere
source("spatial/spatial_figures.R", echo = FALSE)
stopifnot(SPATIAL_PLOTS_DIR == Sys.getenv("SPATIAL_PLOTS_DIR"))
cat("spatial source data done\n")
