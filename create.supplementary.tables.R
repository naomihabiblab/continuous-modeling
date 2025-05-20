library(openxlsx)

source("proteomics.R")
source("topic_space_functions.R")

creator <- "Roi Meir, Habib Lab, HUJI"

SUPPLEMENTARY_DATA_DIR <- file.path(get_result_dir(), "supplementary_data")
dir.create(SUPPLEMENTARY_DATA_DIR, showWarnings = FALSE)

topics_h5ad_path <- file.path(get_result_dir(), 'topic_space', '500_topics_all.h5ad')
data <- anndata::read_h5ad(topics_h5ad_path)

ROSMAP_CLINICAL_DATA_PATH <- file.path(get_result_dir(), "../../..//Shared/NextSeq/500/v1.1.objects.for.synapse/ROSMAP_clinical.csv")

donor.mapping <- read.csv(ROSMAP_CLINICAL_DATA_PATH) %>% dplyr::select(projid, individualID)
donor.mapping <- setNames(donor.mapping$individualID, donor.mapping$projid)

replace_projid_with_individualID <- function(df, column=NA) {
  df <- df %>% as.data.frame()
  if (!is.na(column)) {
    df[column] <- donor.mapping[df[[column]]]
  } else {
    rownames(df) <- donor.mapping[rownames(df)]
  }
  
  df
}

header.style <- function() {
  createStyle(fontSize = 10, halign = "left", fgFill = "#808080", textDecoration = "bold")
}

add.header.style <- function(wb, sheet, ncols) {
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncols+1, gridExpand = TRUE)
}

format.sheet <- function(wb, sheet, df, sets, horizontal.border.after=NULL) {
  columns <- setNames(LETTERS[0:ncol(df)+1], colnames(df))
  
  # identify rows of same state in order to apply different colors in sheet
  sets  <- cumsum(tidyr::replace_na(sets, FALSE))
  
  add.header.style(wb, sheet, ncol(df))
  
  # Even groups style
  addStyle(wb, sheet = sheet, gridExpand = TRUE,
           createStyle(fontSize = 11, fgFill = "#e0dede", borderStyle = "thin"),
           rows = which(as.logical(sets %% 2)) + 1, 
           cols = 0:ncol(df)+1)
  
  if(!is.null(horizontal.border.after)) {
    addStyle(wb, sheet, createStyle(border="Left", borderStyle ="thick"), rows = 0:nrow(df)+1, 
             cols = columns[c(horizontal.border.after)], 
             stack = TRUE, gridExpand = TRUE)
  }
  return(wb)
}


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Supplementary Table 1 - Cohorts Description.              ########
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(function (){
  wb <- createWorkbook(creator = creator, title = "Cohorts description")
  addWorksheet(wb, "Description", gridLines = TRUE)

  # -------------------------------------------------------- #
  # Create sheet for characteristics of snRNA-seq cohort     #
  # -------------------------------------------------------- #
  df <- data$obsm$meta.data %>%
    replace_projid_with_individualID %>%
    rownames_to_column("individualID") %>%
    dplyr::select(individualID, study, msex, age_death, pmi, cogng_demog_slope, cogdx, ceradsc, braaksc, niareagansc,
                  sqrt.amyloid, sqrt.amyloid_mf, sqrt.tangles, sqrt.tangles_mf, mglia3_mf, mglia123_mf,
                  zcapture_syn_3cort,  synap_3cort_vamp,  synap_3cort_syntaxin,  synap_3cort_synaptophys,  synap_3cort_stagmin,
                  synap_3cort_snap25,  synap_3cort_sept5,  synap_3cort_complex2,  synap_3cort_complex1,
                  synap_3cort_capture4,  synap_3cort_capture3,  synap_3cort_capture2,  synap_3cort_capture1,
                  )

  addWorksheet(wb, (sheet <- "snRNA-seq cohort"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)



  # -------------------------------------------------------- #
  # Create sheet for characteristics of proteomics cohort    #
  # -------------------------------------------------------- #
  df <- colData(proteom) %>%
    as.data.frame %>%
    dplyr::select(individualID=IndividualID, batch=Batch, sampleID=SampleID,
                  study, msex, age_death, pmi, cogng_demog_slope=cogng_random_slope, cogdx, ceradsc, braaksc, niareagansc,
                  amyloid_sqrt, sqrt_amyloid_mf, tangles_sqrt, sqrt_tangles_mf)

  addWorksheet(wb, (sheet <- "proteomics cohort"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)



  saveWorkbook(wb, file.path(SUPPLEMENTARY_DATA_DIR, "Supplementary Table 1 - Cohorts description.xlsx"),
               overwrite = TRUE)
})()


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Supplementary Table 2 - Programs Atlas Characterization  ########
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(function () {
  wb <- createWorkbook(creator = creator, title = "Atlas Characterization")
  addWorksheet(wb, "Description", gridLines = TRUE)
  
  cts <- c("inhibitory", "cux2p", "cux2m", "microglia","astrocytes", "oligo", "opcs")
  
  # -------------------------------------------------------- #
  # Create sheet for DEGs of all programs                    #
  # -------------------------------------------------------- #
  df <- lapply(cts, function(ct) {
    h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/de")) %>% 
      mutate(program=paste0(cell.type.to.prefix[[ct]], gsub("^k", ".p", topic))) %>% 
      mutate(cell.type=case_when(cell.type == 'cux2p' ~ "Excitatory L2-3",
                                 cell.type == 'cux2m' ~ "Excitatory L4-6", 
                                 cell.type == 'oligo' ~ 'oligodendrocytes', 
                                 TRUE ~ cell.type)) %>% 
      dplyr::arrange(as.integer(gsub('^k', '', topic))) %>% 
      dplyr::select(cell.type, program, gene, gene_ID=id, gene_score,
                    fastTopics_postmean=postmean, fastTopics_lfsr=lfsr,
                    reweighted_gene_score, log_z_score_reweighted_gene_score=z_score_log)
      
  }) %>%  do.call(rbind, .)
  
  addWorksheet(wb, "DEGs", gridLines = TRUE)
  writeData(wb, sheet = "DEGs", df, rowNames = FALSE)
  wb <- format.sheet(wb, "DEGs", df, df %>% mutate(v = program != lag(program)) %>% pull(v))
  
  
  
  # -------------------------------------------------------- #
  # Create sheet for pathways of all programs.               #
  # -------------------------------------------------------- #
  df <- lapply(cts, function(ct)
    h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[[ct]], "/pathways")) %>% 
      mutate(program=paste0(cell.type.to.prefix[[ct]], gsub("^k", ".p", topic))) %>% 
      mutate(cell.type=case_when(cell.type == 'cux2p' ~ "Excitatory L2-3",
                                 cell.type == 'cux2m' ~ "Excitatory L4-6", 
                                 cell.type == 'oligo' ~ 'oligodendrocytes', 
                                 TRUE ~ cell.type)) %>% 
      dplyr::arrange(as.integer(gsub('^k', '', topic))) %>% 
      dplyr::select(cell.type, program, pathway=Description, gene=geneName, dataset=db, pathwayID=ID, 
                    Count, GeneRatio, BgRatio, pvalue, p.adjust, qvalue)
    ) %>% do.call(rbind, .)
  
  addWorksheet(wb, "pathways", gridLines = TRUE)
  writeData(wb, sheet = "pathways", df, rowNames = FALSE)
  wb <- format.sheet(wb, "pathways", df, df %>% mutate(v = program != lag(program)) %>% pull(v))
  
  saveWorkbook(wb, file.path(SUPPLEMENTARY_DATA_DIR, "Supplementary Table 2 - Atlas Characterization.xlsx"), overwrite = TRUE)
})()


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Supplementary Table 3 - Abundance and association.        ########
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(function (){
  wb <- createWorkbook(creator = creator, title = "Programs Abundance and Associations")
  addWorksheet(wb, "Description", gridLines = TRUE)

  # -------------------------------------------------------- #
  # Program abundance matrix                                 #
  # -------------------------------------------------------- #
  
  df <- replace_projid_with_individualID(data$X) %>% dplyr::rename_with(~topic_names_mapping(.x, use_static_mapping=FALSE))
  
  addWorksheet(wb, (sheet <- "programs abundance"), gridLines = TRUE)
  writeData(wb, sheet = sheet, 
            df %>% rownames_to_column("donor.id"), 
            startRow = 1, rowNames = FALSE)
  
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(data$X)+1, gridExpand = TRUE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1:nrow(data$X)+1, cols = 1, gridExpand = TRUE)
  
  
  # -------------------------------------------------------- #
  # Programs-endophenotype associations                      #
  # -------------------------------------------------------- #
  df <- data$uns$trait.analysis %>% py_to_r %>%
    mutate(program=topic_names_mapping(covariate, use_static_mapping=FALSE)) %>% 
    dplyr::select(trait, program, beta, se, tstat, pval, r.sq, n, formula, adj.pval, sig)

  addWorksheet(wb, (sheet <- "Endophenoptye associations"), gridLines = TRUE)
  writeData(wb, sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)
  
  # -------------------------------------------------------- #
  # Programs-ELISA associations                              #
  # -------------------------------------------------------- #
  df <- data$uns$syn.trait.analysis %>% py_to_r %>%
    mutate(program=topic_names_mapping(covariate, use_static_mapping=FALSE)) %>% 
    dplyr::select(trait, program, beta, se, tstat, pval, r.sq, n, formula, adj.pval, sig)
  
  addWorksheet(wb, (sheet <- "ELISA associations"), gridLines = TRUE)
  writeData(wb, sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)
  
  # -------------------------------------------------------- #
  # ELISA-Endophenoptye associations                         #
  # -------------------------------------------------------- #
  df <- data$uns$disease.syn.trait.analysis %>% py_to_r %>%
    dplyr::select(trait, covariate, beta, se, tstat, pval, r.sq, n, formula, adj.pval, sig)
  
  addWorksheet(wb, (sheet <- "ELISA-phenoptye associations"), gridLines = TRUE)
  writeData(wb, sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)
  
  saveWorkbook(wb, file.path(SUPPLEMENTARY_DATA_DIR, "Supplementary Table 3 - Abundance and Endophenotypes Associations.xlsx"),
               overwrite = TRUE)
})()

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Supplementary Table 4 - BEYOND analysis                   ########
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(function (){
  wb <- createWorkbook(creator = creator, title = "BEYOND analysis results")
  addWorksheet(wb, "Description", gridLines = TRUE)
  
  # -------------------------------------------------------- #
  # Cellular landscape 2d embedding                          #
  # -------------------------------------------------------- #
  df <- data.frame(donor.id = donor.mapping[data$obs_names],
                   data$obsm$X_phate %>% `colnames<-`(paste0("PHATE_", 1:2)),
                   data$obsm$X_umap %>% `colnames<-`(paste0("UMAP_", 1:2)))
  
  addWorksheet(wb, (sheet <- "2D Landscape embedding"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  
  # -------------------------------------------------------- #
  # Donors' pseudotime and trajectory probability            #
  # -------------------------------------------------------- #
  trajectories <- data$uns$trajectories$palantir
  terminals <- trajectories$terminals %>% py_to_r %>% dplyr::select(x=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("x")
  root      <- setNames(c(T), trajectories$user.root)
  
  df <- data.frame(pseudotime = trajectories$pseudotime,
                   trajectories$branch.probs %>% py_to_r) %>%
    rownames_to_column("individualID") %>%
    mutate(terminal = terminals[individualID,"terminal"],
           root = root[individualID],
           individualID = donor.mapping[individualID])
  
  
  addWorksheet(wb, (sheet <-"Trajectories"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  
  # -------------------------------------------------------- #
  # Subpopulation community assignment                       #
  # -------------------------------------------------------- #
  df <- data$var %>% dplyr::select(community, sub.community) %>%
    rownames_to_column("program") %>%
    mutate(program=topic_names_mapping(program, use_static_mapping=FALSE))

  addWorksheet(wb, (sheet <- "Program communities"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


  # -------------------------------------------------------- #
  # Donor community proportions                              #
  # -------------------------------------------------------- #
  df <- cbind(data$obsm$communities %>% dplyr::rename(C.NA='NA.'),
              data$obsm$sub.communities %>% dplyr::rename(SUB.NA='NA.')) %>%
    replace_projid_with_individualID() %>%
    rownames_to_column("donor.id") %>% 
    dplyr::rename()

  addWorksheet(wb, (sheet <- "Participant comm. prop."), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


  # -------------------------------------------------------------- #
  # Community-endophenotype associations                           #
  # -------------------------------------------------------------- #
  df <- data$uns$communities$trait.association %>% py_to_r 
  addWorksheet(wb, (sheet <- "Community endophe. assoc."), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


  # ------------------------------------------------------------ #
  # Program, community, ELISA and endophenotype dynamics - fitted values  #
  # ------------------------------------------------------------ #
  df <- rbind(data$uns$trajectories$palantir$dynamics$fitted.vals %>% py_to_r,
              data$uns$communities$dynamics$fitted.vals %>% py_to_r,
              data$uns$trajectories$palantir$syn_dynamics$fitted.vals %>% py_to_r) %>%
    dplyr::select(feature, trajectory, pseudotime=x, value=y, weight=X.weights., fit, se.fit) %>%
    mutate(trajectory = recode(trajectory, "ABA"="Alternative Aging", "prAD" = "prAD"))
  
  programs <- df %>% filter(str_detect(feature, "\\.k\\d+")) %>% pull(feature) %>% unique
  programs_names_mapping <- topic_names_mapping(programs, use_static_mapping = FALSE)
  names(programs_names_mapping) <- programs
  
  df <- df %>% mutate(feature=ifelse(feature %in% programs, programs_names_mapping[feature], feature))

  addWorksheet(wb, (sheet <- "Dynamics fitted values"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


  # -------------------------------------------------------------- #
  # Program, community, ELISA and endophenotype dynamics - predicted values #
  # -------------------------------------------------------------- #
  df <- rbind(data$uns$trajectories$palantir$dynamics$pred.vals %>% py_to_r,
              data$uns$communities$dynamics$pred.vals %>% py_to_r,
              data$uns$trajectories$palantir$syn_dynamics$pred.vals %>% py_to_r) %>%
    dplyr::select(feature, trajectory, pseudotime=x, fit, se.fit) %>%
    mutate(trajectory = recode(trajectory, "ABA"="Alternative Aging", "prAD" = "prAD"))
  
  programs <- df %>% filter(str_detect(feature, "\\.k\\d+")) %>% pull(feature) %>% unique
  programs_names_mapping <- topic_names_mapping(programs, use_static_mapping = FALSE)
  names(programs_names_mapping) <- programs
  
  df <- df %>% mutate(feature=ifelse(feature %in% programs, programs_names_mapping[feature], feature))

  addWorksheet(wb, (sheet <- "Dynamics predicted values"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


  # ------------------------------------------------------------ #
  # Program dynamics over Green et al. pseudotime - fitted values.   #
  # ------------------------------------------------------------ #
  df <- data$uns$trajectories$da_fits_sqrt$fitted.vals %>% 
    py_to_r %>% 
    dplyr::select(feature, trajectory, pseudotime=x, value=y, weight=X.weights., fit, se.fit) %>%
    mutate(trajectory = recode(trajectory, "ABA"="Green et al. Alternative Aging", "prAD" = "Green et al. prAD"))
  
  programs <- df %>% filter(str_detect(feature, "\\.k\\d+")) %>% pull(feature) %>% unique
  programs_names_mapping <- topic_names_mapping(programs, use_static_mapping = FALSE)
  names(programs_names_mapping) <- programs
  
  df <- df %>% mutate(feature=ifelse(feature %in% programs, programs_names_mapping[feature], feature))
  
  addWorksheet(wb, (sheet <- "Dyn. fit. Green pseudotime"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  
  # -------------------------------------------------------------- #
  # Program dynamics over Green et al. pseudotime - predicted values.  #
  # -------------------------------------------------------------- #
  df <- data$uns$trajectories$da_fits_sqrt$pred.vals %>% 
    py_to_r %>% 
    dplyr::select(feature, trajectory, pseudotime=x, fit, se.fit) %>%
    mutate(trajectory = recode(trajectory, "ABA"="Alternative Aging", "Green et al. prAD" = "Green et al. prAD"))
  
  programs <- df %>% filter(str_detect(feature, "\\.k\\d+")) %>% pull(feature) %>% unique
  programs_names_mapping <- topic_names_mapping(programs, use_static_mapping = FALSE)
  names(programs_names_mapping) <- programs
  
  df <- df %>% mutate(feature=ifelse(feature %in% programs, programs_names_mapping[feature], feature))
  
  addWorksheet(wb, (sheet <- "Dyn. pred. Green pseudotime"), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  saveWorkbook(wb, file.path(SUPPLEMENTARY_DATA_DIR, "Supplementary Table 4 - BEYOND analysis results.xlsx"),
               overwrite = TRUE)
})()


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Supplementary Table 5 - AD-associated programs            ########
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(function (){
  wb <- createWorkbook(creator = creator, title = "AD-associated programs")
  addWorksheet(wb, "Description", gridLines = TRUE)
  
  # -------------------------------------------------------- #
  # Proportions of cell cycle reentry and apoptosis neurons   #
  # -------------------------------------------------------- #

  cts <- c(inhibitory="inhibitory", "excitatory L2-3"="cux2p", "excitatory L4-6"="cux2m")
  for (ct in names(cts)) {
    prop_cells_df <- get_AD_programs_proportion(cts[[ct]])
    
    
    df <- prop_cells_df %>% 
      filter(projid %in% rownames(data)) %>% 
      mutate(group=recode_factor(group, start='healthy')) %>% 
      replace_projid_with_individualID(column = 'projid') %>% 
      dplyr::rename_with(~gsub("cell_cycle", "cell_cycle_reentry", .x)) %>% 
      dplyr::rename_with(~gsub("apoptosis", "stress_apoptosis", .x)) %>% 
      dplyr::rename(individualID=projid)
      
    addWorksheet(wb, (sheet <- glue::glue("Neuronal prop. {ct}")), gridLines = TRUE)
    writeData(wb, sheet = sheet, df, rowNames = FALSE)
    addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  }

  
  saveWorkbook(wb, file.path(SUPPLEMENTARY_DATA_DIR, "Supplementary Table 5 - AD-associated programs.xlsx"),
               overwrite = TRUE)
})()
 

