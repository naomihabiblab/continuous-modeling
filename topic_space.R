library(dplyr)
library(ComplexHeatmap)
library(anndata)
library(stringr)

source("topic_space_functions.R")
source("utils.R")
source("data_500_utils.R")
source("landscape.analysis.funs.R")
source("figures.funs.R")
source("topics.R")

result_dir <- get_result_dir()
topic_space_plots_dir <- file.path(result_dir, 'topic_space')
dir.create(topic_space_plots_dir, showWarnings = FALSE)

topics_h5ad_path <- file.path(topic_space_plots_dir, '500_topics_all.h5ad')

#++++++++++++++++++++++++++++++++#
# Create topic space object ######
#++++++++++++++++++++++++++++++++#



topic_space <- read_topic_space()
data.da <- anndata::read_h5ad(AD_500_ANDATA_FILE)
topic_space <- topic_space %>% filter(!is.na(projid)) %>% tibble::column_to_rownames("projid")

# Shared donors
data.da <- data.da[intersect(rownames(data.da), rownames(topic_space)),]
topic_space <- topic_space[rownames(data.da),]


data <- AnnData(X=data$X[,!grepl("PeriFibSMC|Endothelial", colnames(data$X))], obsm = list(meta.data=data.da$obsm$meta.data, QCs=data.da$obsm$QCs))
data$layers['sqrt.prev'] <- sqrt(data$X)

anndata::write_h5ad(data, topics_h5ad_path)
rm(data)

#+++++++++++++++++++++++++++++++++#
# Topic to Topic correlation  #####
#+++++++++++++++++++++++++++++++++#

data <- anndata::read_h5ad(topics_h5ad_path)

data$uns$tt.cor <- correlate_prevalence(data$X)
data$uns$ts.cor <- correlate_prevalence(cbind(data$X, data.da$X))

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsCorrelation.pdf"), height = 50, width = 40, {
  plot_landscape_correaltion(data$uns$tt.cor,
                             title = "Topics correlations") %>% draw
})

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsCorrelationNeuronal.pdf"), height = 50, width = 40, {
  plot_landscape_correaltion(data$uns$tt.cor,
                             title = "Neuronal topics correlations",
                             names = grepl("Inh|Exc", data$uns$tt.cor$names)) %>% draw
})

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsToStateCorrelation.pdf"), height = 70, width = 60, {
  plot_landscape_correaltion(data$uns$ts.cor,
                             title = "Topics to state correlations",
                             add_significance = FALSE) %>% draw
})

anndata::write_h5ad(data, topics_h5ad_path)
rm(data)

#+++++++++++++++++++++++#
# Traits association ####
#+++++++++++++++++++++++#

data <- anndata::read_h5ad(topics_h5ad_path)

traits_columns <- c("sqrt.amyloid", "sqrt.amyloid_mf", "sqrt.tangles", "sqrt.tangles_mf", "cogng_demog_slope")
traits <- data$obsm$meta.data[,c(traits_columns)] %>% mutate(across(all_of(traits_columns), ~as.numeric(as.character(.))))

controls <- data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")])
controls_with_qc <- data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")],
                       data$obsm$QCs[,c("Total_Genes_Detected","Estimated_Number_of_Cells")])

data$uns$trait.analysis <- associate.traits(traits,
                                            data$layers[['sqrt.prev']],
                                            controls) %>% dplyr::mutate("topic"="covariate")

data$uns$trait.analysis.with.qc <- associate.traits(traits,
                                                    data$layers[['sqrt.prev']],
                                                    controls_with_qc) %>% dplyr::mutate("topic"="covariate")

withr::with_pdf(file.path(topic_space_plots_dir, "TraitAssocationHeatmap.pdf"), {
  plot.trait.associations.heatmap(data$uns$trait.analysis , params = traits_columns, row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitAssocationNeuronsHeatmap.pdf"), {
  plot.trait.associations.heatmap(data$uns$trait.analysis %>% py_to_r %>% filter(grepl("Inh|Exc", covariate)) ,
                                  params = c("sqrt.amyloid", "sqrt.tangles", "sqrt.amyloid_mf",  "sqrt.tangles_mf", "cogng_demog_slope"), row.by='covariate')
})


withr::with_pdf(file.path(topic_space_plots_dir, "QCAssocationHeatmap.pdf"), height = 15, {
  associate.traits(data$obsm$QCs[,c("Total_Genes_Detected","Estimated_Number_of_Cells")],
                   data$layers[['sqrt.prev']],
                   controls) %>%
    dplyr::mutate("topic"="covariate") %>%
    arrange(adj.pval) %>%
    plot.trait.associations.heatmap(params=c("Total_Genes_Detected", "Estimated_Number_of_Cells"),
                                    row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitAssocationHeatmapWithQC.pdf"), {
  plot.trait.associations.heatmap(data$uns$trait.analysis.with.qc , params = traits_columns, row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitAssocationNeuronsHeatmapWithQC.pdf"), {
  plot.trait.associations.heatmap(data$uns$trait.analysis.with.qc %>% py_to_r %>% filter(grepl("Inh|Exc", covariate)) ,
                                  params = c("sqrt.amyloid", "sqrt.tangles", "sqrt.amyloid_mf",  "sqrt.tangles_mf", "cogng_demog_slope"), row.by='covariate')
})

syn_traits_col <- c("zcapture_syn_3cort",
                "synap_3cort_vamp",
                "synap_3cort_syntaxin",
                "synap_3cort_synaptophys",
                "synap_3cort_stagmin",
                "synap_3cort_snap25",
                "synap_3cort_sept5",
                "synap_3cort_complex2",
                "synap_3cort_complex1",
                "synap_3cort_capture4",
                "synap_3cort_capture3",
                "synap_3cort_capture2",
                "synap_3cort_capture1")
syn_traits <- data$obsm$meta.data[,c(syn_traits_col)] %>% mutate(across(all_of(syn_traits_col), ~as.numeric(as.character(.))))

data$uns$syn.trait.analysis <- associate.traits(syn_traits,
                                                data$layers[['sqrt.prev']],
                                                controls) %>% dplyr::mutate("topic"="covariate")


data$uns$disease.syn.trait.analysis <- associate.traits(traits,
                                                        syn_traits,
                                                controls)

data$uns$disease.syn.trait.analysis.with.qc <- associate.traits(traits,
                                                        syn_traits,
                                                controls_with_qc)

data$uns$syn.trait.analysis.with.qc <- associate.traits(syn_traits,
                                                        data$layers[['sqrt.prev']],
                                                        controls_with_qc)


data.da$uns$syn.trait.analysis <- associate.traits(syn_traits,
                                            data.da$layers[['sqrt.prev']],
                                            controls) %>% dplyr::mutate("state"="covariate")

withr::with_pdf(file.path(topic_space_plots_dir, "TraitSynAssocationHeatmap.pdf"), {
   plot.trait.associations.heatmap(data$uns$syn.trait.analysis , params = syn_traits_col, row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitSynAssocationHeatmapWithQC.pdf"), {
  plot.trait.associations.heatmap(data$uns$syn.trait.analysis.with.qc %>% py_to_r %>%
                                    filter(grepl("Inh|Exc", covariate))
                                  , params = syn_traits_col, row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitSynAssocationHeatmapState.pdf"), {
  plot.trait.associations.heatmap(data.da$uns$syn.trait.analysis , params = syn_traits_col, row.by='covariate')
})

withr::with_pdf(file.path(topic_space_plots_dir, "TraitDieaseSynAssocationHeatmap.pdf"), {
  plot.trait.associations.heatmap(data$uns$disease.syn.trait.analysis , params = traits_columns, row.by='covariate', show.only.significant = F)
})



global_features <- cbind(wc.prev=data.da$uns$cell.types$wc.prev %>% py_to_r(), 
                         total.prev=data.da$uns$cell.types$prev %>% py_to_r(),
                         data.da$obsm$QCs %>% dplyr::select(Total_Genes_Detected, Estimated_Number_of_Cells),
                         total_number_of_cells=data.da$uns$cell.types$counts %>% py_to_r() %>% rowSums())
colnames(global_features) <- make.names(colnames(global_features))




global_features.trait.analysis <- associate.traits(traits,
                                                   global_features,
                                                   controls) %>% data.frame()

withr::with_pdf(file.path(topic_space_plots_dir, "TraitAssocationHeatmapGlobalFeatures.pdf"), {
  plot.trait.associations.heatmap(global_features.trait.analysis, params = traits_columns, row.by='covariate', show.only.significant = F)
})


anndata::write_h5ad(data, topics_h5ad_path)
rm(data)


#+++++++++++++++++++++++++++#
# Fit State Dynamics  #######
#+++++++++++++++++++++++++++#

data <- anndata::read_h5ad(topics_h5ad_path)

data <- fit_west_south_dnynamics(data = data, data.da = data.da, bootstrap = FALSE)

anndata::write_h5ad(data, topics_h5ad_path)


# Plots dynamics
withr::with_pdf(file.path(topic_space_plots_dir, "TopicDynamicsStatePsuedoTime.pdf"), {
   for (topic in colnames(data)) {
      print(plot.dynamics.wrapper(data$uns$trajectories$da_fits_sqrt, features = c(topic), overlap.pseudotime=.11, ncol=1, strip.position="left", scales="free_x", label=TRUE, include.points=FALSE) +
                        theme(strip.text = element_blank(),
                              legend.position="none",
                              axis.line = element_line(),
                              axis.text.y = element_text()) +
                        labs(x=NULL, y=NULL, title=NULL))
     grid::grid.text(topic,x = (0.2), y = (0.95))
    }
})




#+++++++++++++++++++++++++++++++++++++++++
# Topics over states landscape phates ####
#+++++++++++++++++++++++++++++++++++++++++

df <- data.frame(data$X)

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsOverStatePhate.pdf"), width=20, height=5 * ncol(df) / 4, {
  plot.landscape(df, data=data.da, embedding='X_core_phate', enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=4) %>% print(.)
})



#++++++++++++++++++++++++++++++++++++++++++++++++++++
## Topics space clusters and 2/3d representation ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++

sc <- reticulate::import("scanpy")


sc$pp$neighbors(data, n_neighbors = as.integer(10), use_rep = "X", metric = "cosine")
data <- sc$tl$leiden(data, resolution =.25, copy=TRUE)
sc$tl$leiden(data, resolution = .75, restrict_to=reticulate::tuple("leiden", reticulate::np_array(c("0"))))
data$obs["clusters"] = plyr::mapvalues(data$obs$leiden_R, levels(data$obs$leiden_R), 1:length(levels(data$obs$leiden_R)))
data$obs$leiden <- data$obs$leiden_R <- NULL


sc$tl$tsne(data, n_pcs = 0, use_rep = "X", learning_rate = 100)
sc$tl$umap(data, maxiter = as.integer(3000), spread = 3)


sc$external$tl$phate(data,
                     n_components = as.integer(3),
                     k = as.integer(10), a = as.integer(40),
                     knn_dist =  "euclidean", mds_dist = "euclidean",
                     mds_solver = "smacof", verbose = F)
data$obsm$X_all_3d_phate <- data$obsm$X_phate
data$obsm$X_phate <- NULL

sc$external$tl$phate(data, n_components = as.integer(2),
                     k = as.integer(15), a = as.integer(100),
                     knn_dist =  "euclidean", mds_dist = "correlation",
                     mds_solver = "smacof", verbose = F)



for(e in c("X_all_3d_phate","X_phate","X_umap","X_tsne")) {
  data$obsp[[paste0("similarity_", e)]] <- embedding.similatity(data$obsm[[e]], knn = 5)
}


anndata::write_h5ad(data, topics_h5ad_path)

rm(e, sc)


## Plots over phate and umap ####

df <- data.frame(data$obs["clusters"]) 

withr::with_pdf(file.path(topic_space_plots_dir, "TopicSpaceClusterOverStatePhate.pdf"), {
  plot.landscape.3D(df, theta = 135, phi = 30, data. = data.da, cols = scales::hue_pal(), legend.position = "left") %>% print()
})

withr::with_pdf(file.path(topic_space_plots_dir, "StateSpaceClusterOverStatePhate.pdf"), {
  plot.landscape.3D(data.da$obs["clusters"], theta = 135, phi = 30, data. = data.da, cols = scales::hue_pal(), legend.position = "left") %>% print()
})


df <- data.frame(data$X,
                 data.da$X,
                 data.da$uns$trajectories$palantir$branch.probs %>% py_to_r,
                 data.da$uns$trajectories$palantir$pseudotime)

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsOverTopicPhate.pdf"), width=20, height=5 * ncol(df) / 4, {
  plot.landscape(df, data=data, embedding='X_phate', enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=4) %>% print(.)
})

withr::with_pdf(file.path(topic_space_plots_dir, "TopicsOverTopicUmap.pdf"), width=20, height=5 * ncol(df) / 4, {
  plot.landscape(df, data=data, embedding='X_umap', enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=4) %>% print(.)
})

withr::with_pdf(file.path(topic_space_plots_dir, "TopicSpaceClusterTopicUmap.pdf"), width=5, height=5, {
  plot.landscape(data$obs["clusters"], data=data, embedding='X_umap', cols = scales::hue_pal(), legend.position = "left") %>% print(.)
})

rm(data)

#++++++++++++++++++++++++
# Fits trajectories #####
#++++++++++++++++++++++++

data <- anndata::read_h5ad(topics_h5ad_path)

# Run Palantir algorithm
data$uns$trajectories$palantir <- fit.trajectories.palantir(data, dm=5, dm.k=30,
                                                            palantir.k = 30, scale=FALSE,
                                                            use_early_cell_as_start=FALSE,
                                                            early_cell = "4330337")




data$uns$trajectories$palantir_core <- fit.trajectories.palantir(data, dm=5, dm.k=30,
                                                                 palantir.k = 15, scale=FALSE,
                                                                 use_early_cell_as_start=FALSE,
                                                                 early_cell = "4330337",
                                                                 include=!is.na(data.da$uns$trajectories$palantir$pseudotime))

change_trajectoies_names <- function(trajectories) {
  if (all(trajectories$branch.probs$columns == c('X20779834', 'X3380931'))) {
    trajectories$branch.probs$columns <- c("ABA", "prAD")
    trajectories$terminals$index <- c("ABA", "prAD")
  } else if ((all(trajectories$branch.probs$columns == c('X20779834', 'X76647134')))) {
    trajectories$branch.probs$columns <- c("ABA", "prAD")
    trajectories$terminals$index <- c("ABA", "prAD")
  } else {
    stop("Unknown terminal points")
  }
  
  trajectories
}

# change trajectories names
change_trajectoies_names(data$uns$trajectories$palantir_core)
change_trajectoies_names(data$uns$trajectories$palantir)


df <- data.frame(state=data.da$uns$trajectories$palantir$branch.probs %>% py_to_r,
                 topic=data$uns$trajectories$palantir$branch.probs %>% py_to_r,
                 state.pseudotime=data.da$uns$trajectories$palantir$pseudotime, 
                 topic.pseudotime=data$uns$trajectories$palantir$pseudotime)



plot.landscape(df, data=data, embedding='X_phate', enforce.same.color.scale = FALSE, smoothened = TRUE,
               ncol=2, cols = colorRampPalette(c("#4db6ac","#ef6c00"))) %>% print(.)

withr::with_jpeg(file.path(topic_space_plots_dir, "TopicsTrajectoriesOverPhate.jpg"), {
  plot.landscape(df, data=data, embedding='X_phate', enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=2, cols = colorRampPalette(c("#4db6ac", "white", "#ef6c00"))) %>% print(.)
})


withr::with_jpeg(file.path(topic_space_plots_dir, "TopicsTrajectoriesOverStatePhate.jpg"), {
  plot.landscape(df, data=data.da, embedding='X_core_phate', enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=2) %>% print(.)
})

withr::with_pdf(file.path(topic_space_plots_dir, 'PseudotimeVsBranchProb.pdf'), {
  plot_grid(
    ggplot(df)  + geom_point(aes(x=topic.pseudotime, y=topic.prAD, alpha=is.na(state.pseudotime))),
    ggplot(df)  + geom_point(aes(x=state.pseudotime, y=state.prAD, alpha=is.na(state.pseudotime))),
    ncol = 1
  )  %>% print
})


anndata::write_h5ad(data, topics_h5ad_path)
rm(data)

#+++++++++++++++
# Dynamics  ####
#+++++++++++++++

data <- anndata::read_h5ad(topics_h5ad_path)

features <- data.frame(data$layers[["sqrt.prev"]],
                       data.da$layers[["sqrt.prev"]],
                       data$obsm$meta.data[,c("cogng_demog_slope","sqrt.tangles_mf","sqrt.amyloid_mf", "sqrt.tangles", "sqrt.amyloid",
                                              "age_bl","age_death","msex")], 
                       row.names = data$obs_names)

data$uns$trajectories$palantir$dynamics <-
  fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
               features = features,
               trajectory.probs = data$uns$trajectories$palantir$branch.probs %>% py_to_r(),
               trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
               evaluate.fit = T,
               bootstrap = F)

data$uns$trajectories$palantir$dynamics_core <-
  fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime[data.da$obs$core],
               features = features[data.da$obs$core,],
               trajectory.probs = (data$uns$trajectories$palantir$branch.probs %>% py_to_r())[data.da$obs$core,],
               trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
               evaluate.fit = T,
               bootstrap = F)

anndata::write_h5ad(data, topics_h5ad_path)


syn_features <- data.frame(data$obsm$meta.data[,syn_traits_col], row.names = data$obs_names)
syn_and_all_features <- cbind(features, syn_features)



data$uns$trajectories$syn_dynamics_state <-
    fit.dynamics(pseudotime = data.da$uns$trajectories$palantir$pseudotime,
                 features = syn_features,
                 trajectory.probs = py_to_r(data.da$uns$trajectories$palantir$branch.probs),
                 trajectory.terminal.pseudotime = setNames(py_to_r(data.da$uns$trajectories$palantir$terminals)[,1], data.da$uns$trajectories$palantir$terminals$index),
                 evaluate.fit = T,
                 bootstrap = F)

data$uns$trajectories$palantir$syn_dynamics <-
  fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
               features = syn_features,
               trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
               trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
               evaluate.fit = T,
               bootstrap = F)


data.da$uns$trajectories$palantir$dynamics_with_topics <-
  fit.dynamics(pseudotime = data.da$uns$trajectories$palantir$pseudotime,
               features = features,
               trajectory.probs = data.da$uns$trajectories$palantir$branch.probs %>% py_to_r(),
               trajectory.terminal.pseudotime = setNames(py_to_r(data.da$uns$trajectories$palantir$terminals)[,1], data.da$uns$trajectories$palantir$terminals$index),
               evaluate.fit = T,
               bootstrap = F)



data$uns$trajectories$palantir$syn_and_all_dynamics <-
  fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
               features = syn_and_all_features,
               trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
               trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
               evaluate.fit = T,
               bootstrap = F)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Construct communities based on state correlations and dynamics ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


topics_snr <- py_to_r(data$uns$trajectories$palantir$dynamics$pred.vals) %>%
  filter(feature %in% colnames(data)) %>% 
  dplyr::group_by(feature, trajectory) %>%
  dplyr::summarise(max_fit=max(fit), min_fit=min(fit), mean_se=mean(se.fit), 
                   max_se=max(se.fit),
                   diff_fit=max_fit-min_fit, snr=diff_fit/max_se, snr_mean=diff_fit/mean_se)


dynamical_topics_10 <- topics_snr %>%
  filter(snr > 10) %>%
  pull(feature) %>%
  unique

assoicated_topics <- data$uns$trait.analysis %>%
  py_to_r %>%
  filter(sig != '') %>%
  pull(covariate) %>%
  unique %>% 
  sort

data <- build_communities(data, c(assoicated_topics, dynamical_topics_10) %>% unique,
                                                    n.sub.communities = c(2,3, 2)
                          )
plot_communities(data=data, communities_plots_dir = file.path(topic_space_plots_dir, "CommunitiesDynamical10AndAssociatedTopics"))

anndata::write_h5ad(data, topics_h5ad_path)


# Cell cycle vs apoptosis summary data

(function () {
  START_MAX_PSEUDOTIME <- 0.25
  PR_AD_END_MIN_PROBABILITY <- 0.9

  cts <- c("inhibitory", "cux2p", "cux2m")
  for (ct in cts) {
    ct_number <- which(cts == ct)

    L <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/L"))

    signatures <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/genes.signatures"))
    df <- L %>% left_join(signatures, by='cell')

    df <- df %>%
      left_join(data.frame(py_to_r(data$uns$trajectories$palantir$branch.probs),
                           pseudotime=data$uns$trajectories$palantir$pseudotime) %>% rownames_to_column('projid'))

    df <- df %>% mutate(group=case_when(prAD > PR_AD_END_MIN_PROBABILITY ~ 'prAD',
                                        pseudotime < START_MAX_PSEUDOTIME ~ 'start'))

    allen <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/allen.annotation"))
    df <- df %>% left_join(allen, by=c('cell'))

    apoptosis_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_early)[[ct_number]])
    cell_cycle_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$increasing_late)[[ct_number]])
    decreasing_early_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$decreasing_early)[[ct_number]])
    decreasing_late_topic <- gsub("\\w+\\.", "", names(neuronal_topics_to_annotation$decreasing_late)[[ct_number]])

    df <- df %>% mutate(decreasing_early_topic_score=.data[[decreasing_early_topic]],
                        decreasing_late_topic_score=.data[[decreasing_late_topic]],
                        cell_cycle_topic_score=.data[[cell_cycle_topic]],
                        apoptosis_topic_score=.data[[apoptosis_topic]])


    quantiles_apoptosis <- topic_quantiles_healthy(df, apoptosis_topic)
    quantiles_cell_cycle <- topic_quantiles_healthy(df, cell_cycle_topic)

    df <- df %>% mutate(
      `cell_cycle_95`=as.integer(cell_cycle_topic_score > quantiles_cell_cycle['95%', 'threshold']),
      `cell_cycle_90`=as.integer(cell_cycle_topic_score > quantiles_cell_cycle['90%', 'threshold']),
      `cell_cycle_80`=as.integer(cell_cycle_topic_score > quantiles_cell_cycle['80%', 'threshold']),
      `apoptosis_95`=as.integer(apoptosis_topic_score > quantiles_apoptosis['95%', 'threshold']),
      `apoptosis_90`=as.integer(apoptosis_topic_score > quantiles_apoptosis['90%', 'threshold']),
      `apoptosis_80`=as.integer(apoptosis_topic_score > quantiles_apoptosis['80%', 'threshold']),
    )

    quantiles_cell_cycle <- quantiles_cell_cycle %>% rownames_to_column("start_quantile")
    quantiles_apoptosis <- quantiles_apoptosis %>% rownames_to_column("start_quantile")
    
    df <- df %>% dplyr::select(-starts_with('k'))

    
    tryCatch({
    h5delete(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))
    }, error=function(e){})
    h5write(df, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_apoptosis"))

    attr(quantiles_cell_cycle, 'start_max_pseudotime') <- START_MAX_PSEUDOTIME
    attr(quantiles_cell_cycle, 'pr_ad_end_min_probability') <- PR_AD_END_MIN_PROBABILITY
    tryCatch({
    h5delete(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_quantiles"))
  }, error=function(e){})
    h5write(quantiles_cell_cycle, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/cell_cycle_quantiles"),  write.attributes=TRUE)

    attr(quantiles_apoptosis, 'start_max_pseudotime') <- START_MAX_PSEUDOTIME
    attr(quantiles_apoptosis, 'pr_ad_end_min_probability') <- PR_AD_END_MIN_PROBABILITY
    tryCatch({
    h5delete(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/apoptosis_quantiles"))
}, error=function(e){})
    h5write(quantiles_apoptosis, TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/apoptosis_quantiles"),  write.attributes=TRUE)
  }
})()
