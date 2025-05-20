neuronal_markers <- list(
    GABAergic=c("MEG3", "PVALB","SST","VIP", "KIT","GAD2", "GAD1"),
    Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2", "NXPH3", "PCDH17"),
    GeneralMarkers=c("NPY", "SST", "VIP", "PVALB", "KIT", "ADARB2", "PAX6", "LAMP5", "LHX6", "RELN", "NDNF"),
    Neurogliaform=c("RELN", "NDNF", "ACTN2", "NR2F2", "NPY", "NOS1"),
    Inhibitory=c("PVALB", "VIP", "TAC3", "CALB2", "SST", "TRHDE", "RALYL", "RELN", "CXCL14", "CNR1", "FBXL7", "NRG1", "KIT", "FGF13", "TRPS1", "CAB", "PTPRK", "NOS1", "TOX", "NTNG1", "CACNA2D1"),
    Excitatory=c("SLC17A7", "RORB", "HS3ST4", "PTPN3", "RTKN2", "ESR1", "COL22A1", "FILIP1L", "SEMA3E", "FOXP2", "DCC",  "DPP4", "THEMIS", "C1QL3",  "NPSR1", "LAMP5", "GLRA3", "PLCH1", "RMST", "OTOGL", "FEZF2", "SCN4B", "IL26"),
    # https://www.nature.com/articles/ncomms11022
    IEG=c("FOS", "EGR1", "EGR2", "EGR3", "ARC", "PROX1", "FOSB", "JUNB", "HOMER1")
  )

marker_analysis <- function(obj, output_path, markers, group.by = NULL, markers_name) {
  number_of_markers <- sapply(markers, length) %>% sum
  number_of_rows <- ceiling(sapply(markers,length) / 2)

  pdf(file = file.path(output_path, paste0(markers_name, "Markers.pdf")), width = 30, height = sum(number_of_rows) * 3.3)

  lapply(names(markers), function (type) {
    print(type)
    genes <- markers[[type]]
    FeaturePlot(obj, genes, ncol=2, cols = viridisLite::turbo(20))
  }) %>% plot_grid(plotlist = ., ncol = 1, labels=paste(markers_name, names(markers), sep = ' - ')) %>% print()
  dev.off()

  pdf(file = file.path(output_path, paste0(markers_name, "MarkersDensity.pdf")), width = 30, height = sum(number_of_rows) * 3.3)
  lapply(names(markers), function (type) {
    print(type)
    genes <- markers[[type]]
    FeatureDensity(obj, features=genes, ncol=2, cols = viridisLite::turbo(20))
  }) %>% plot_grid(plotlist = ., ncol = 1, labels=paste(markers_name, names(markers), sep = ' - ')) %>% print()
  dev.off()

  if (length(group.by) == 1 ) {
    group.by <- c(group.by)
  }
  for (gb in group.by) {
    pdf(file = file.path(output_path, sprintf("%sMarkersDotPlotGroupedBy%s.pdf", markers_name, gb)), height = 1.8 * number_of_markers, width=number_of_markers*0.4)
    lapply(names(markers), function (type) {
      print(type)
      genes <- markers[[type]]

      DotPlot(obj, features=genes, group.by = gb)
    }) %>% plot_grid(plotlist = ., ncol = 1, labels=paste(markers_name, names(markers), sep = ' - ')) %>% print()
    dev.off()
  }
}
