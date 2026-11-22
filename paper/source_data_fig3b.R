# Partial Source Data for Fig. 3b (inhibitory pathway networks, INH.p10 and INH.p8).
#
# The published networks were drawn from a pathway table that had a `group`
# (pathway cluster) column. The current topics_data.h5 table
# (neuronal/inhibitory/pathways) no longer has it, so the cluster colours and
# the exact node subset (only clustered pathways were drawn) cannot be
# reproduced. This script saves what can be reproduced from the stored data:
#   - the enrichment table of each program (node size = -log10 FDR);
#   - the pathway-pathway distances computed exactly as in
#     plot_pathways_overview_graph() (R/pathways.analysis.R), with the same
#     max_distance as in paper/figures.R (k10: 0.3, k8: 0.7), over all pathways
#     of the program (the published graph used only the clustered subset).
# It also reports which of the named pathways in the published figure are present.
#
# Reads only the h5 (no single-cell data):  Rscript paper/source_data_fig3b.R

suppressPackageStartupMessages({ library(rhdf5); library(dplyr); library(reshape2) })
source("R/utils.R")
source("R/data_500_utils.R")
source("R/pathways.analysis.R")
source("R/source_data.R")

pts_all <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[["inhibitory"]], "/pathways"))

# named nodes in paper/figures.R (Figure 3b block)
named <- list(k10=c("Metallothioneins bind metals", "actin binding", "CaM pathway", "microtubule",
                    "synaptic vesicle", "Synaptic vesicle cycle", "Translation",
                    "oxidative phosphorylation", "regulation of cell morphogenesis",
                    "response to oxidative stress", "neuron apoptotic process", "regulation of transporter activity",
                    "rRNA binding", "neuron spine"),
              k8=c("M Phase", "Centrosome maturation", "DNA Repair", "chromatin remodeling", "sister chromatid segregation",
                   "RNA degradation", "RNA splicing"))
programs <- list(k10=list(name="INH.p10", max_distance=0.3), k8=list(name="INH.p8", max_distance=0.7))

for (k in names(programs)) {
  pts <- pts_all %>% filter(topic == k)
  nodes <- pts %>% mutate(log10_FDR = -log10(p.adjust), named_in_figure = Description %in% named[[k]]) %>%
    dplyr::select(ID, Description, db, p.adjust, log10_FDR, Count, GeneRatio, BgRatio, geneID, named_in_figure)
  sd_add(nodes, "Fig3", "b", c(programs[[k]]$name, "nodes"),
         description = glue::glue("{programs[[k]]$name}: enriched pathways (node size = -log10 FDR). PARTIAL: cluster colours (pathway groups) and the exact clustered node subset are no longer stored and cannot be reproduced"))

  dis <- compute_pathways_distance_matrix(pts)          # same function as the figure
  edges <- melt(dis, value.name = "weight")
  colnames(edges) <- c("from", "to", "weight")
  edges <- edges %>% filter(as.character(from) < as.character(to) & weight > 0 & weight < programs[[k]]$max_distance)
  sd_add(edges, "Fig3", "b", c(programs[[k]]$name, "edges"),
         description = glue::glue("{programs[[k]]$name}: pathway pairs with distance (1 - Spearman correlation of gene membership) < {programs[[k]]$max_distance}, computed over all pathways of the program"))

  missing <- setdiff(named[[k]], pts$Description)
  cat(sprintf("%s: %d pathways, %d edges; named pathways missing from current table: %s\n",
              programs[[k]]$name, nrow(pts), nrow(edges), if (length(missing)) paste(missing, collapse = "; ") else "none"))
}
