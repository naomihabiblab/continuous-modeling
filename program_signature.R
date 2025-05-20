library(Seurat)
library(fastTopics)
library(Matrix)
library(ggcorrplot)
library(gplots)
library(ltm)
library(ggplotify)

###############################################################
################# create reweighted F matrix  #################
###############################################################
# Returns a reweighted gene scores F' matrix, based on the new method.
# Defines a scaled gene score which prioritizes for high original gene score and unique genes in each topic.
# f_matrix - the fitted_model$F, the nmf gene score matrix (genes*topics).
# upper_quantile - the upper quantile to use in the formula.
# lower_quantile - the lower quantile to use in the formula.
topic_reweight_f<-function(f_matrix, upper_quantile=0.8, lower_quantile=0.5){
  topics_amount<-ncol(f_matrix)
  reweighted_f<-f_matrix
  for(topic_idx in 1:topics_amount){
    # create a vector with FALSE in each idx that isn't the topic_idx (put there TRUE)
    logical_vec_curr_idx<-logical(topics_amount)
    logical_vec_curr_idx[topic_idx]<-TRUE
    cur_topic_as_mat<-f_matrix[, logical_vec_curr_idx, drop=FALSE]
    all_other_topics_mat<-f_matrix[, !logical_vec_curr_idx, drop=FALSE]
    # find quantile per each row (hence MARGIN=1) => across all other topics
    upper_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, upper_quantile)})
    lower_quantile_other_topics<-apply(all_other_topics_mat, 1, function(x){quantile(x, lower_quantile)})
    log_ratio_ours_vs_upper<-log(cur_topic_as_mat+1e-15) - log(upper_quantile_other_topics+1e-15)
    log_ratio_ours_vs_lower<-log(cur_topic_as_mat+1e-15) - log(lower_quantile_other_topics+1e-15)
    reweighted_f[, topic_idx]<- log_ratio_ours_vs_upper * 
      ifelse(log_ratio_ours_vs_lower<0, Inf, (f_matrix[, logical_vec_curr_idx]*log_ratio_ours_vs_lower)^sign(log_ratio_ours_vs_upper)) 
  } 
  return(reweighted_f)
}

##############################################################################
########## calculate percent expression of each gene in each topic  ##########
##############################################################################
# For each topic (i) and gene (j), calculates the percent of cells in topic i 
# that express gene j. Returns a matrix (genes, topics) of those values.
# count_mat - obj[['RNA']]@counts, matrix with the expression of genes across cells. Must have - rows = genes, cols = cells.
# l_mat - fitted_model$L, the L matrix of the topic modeling. Must have rows = cells, cols = topics.
# number_of_chunks - in case of heavy count matrix, the matrix multiplication can be done by chunks. If the function's
# running time is too long or it collapses, increase this parameter.
percent_cells_exp_gene_in_topic_mat<-function(count_mat, l_mat, number_of_chunks=2){
  bool_mat=(count_mat>0)
  datalist = list()
  row_ind<-seq(1, nrow(bool_mat), ceiling(nrow(bool_mat)/number_of_chunks))
  j<-1
  for(i in row_ind){
    end_ind<-min(i + ceiling(nrow(bool_mat)/number_of_chunks) - 1, nrow(bool_mat))
    datalist[[j]]<-(bool_mat[i:end_ind, ] %*% l_mat)
    j<-j+1
  }
  cell_topics_mat = do.call(rbind, datalist)
  cell_topics_mat<-sweep(as.matrix(cell_topics_mat), 2, colSums(l_mat), "/")
  return(as.matrix(cell_topics_mat))
}

#########################################################################
################# create DEs df with all relevant info  #################
#########################################################################
# Returns a df of positive (reweighted_gene_score>0) DE genes with their name,
# original score, new reweighted score, percent of cells expression that gene in the 
# topic, gene_id and z score across the log of the reweighted scores.
# obj - the object relevant to the topics.
# reweighted_f - the new F', returned from the function - topic_reweight_f().
# fitted_model - the fit of the topics, which includes L, F, etc.
# assay - object assay to take the count from
# organism - for gene name conversion
df_pos_des<-function(obj, reweighted_f, fitted_model, assay="RNA", organism="mmu"){
  df_pos_des<-df_gene_scores(obj, reweighted_f, fitted_model, assay, organism) %>% filter(reweighted_gene_score>0)
  df_pos_des<-df_pos_des %>% group_by(topic) %>% 
    mutate(z_score_log=((log(reweighted_gene_score) - mean(log(reweighted_gene_score))) / sd(log(reweighted_gene_score)))) %>%
    arrange(topic, desc(z_score_log))
  return(df_pos_des)
}

# Returns a df of all genes names, original score, new reweighted score and percent 
# of cells expression that gene in the topic, gene_id.
# obj - the object relevant to the topics.
# reweighted_f - the new F', returned from the function - topic_reweight_f().
# fitted_model - the fit of the topics, which includes L, F, etc.
df_gene_scores<-function(obj, reweighted_f, fitted_model, assay="RNA", organism="mmu"){
  percent_cells<-percent_cells_exp_gene_in_topic_mat(obj[[assay]]@counts[,rownames(fitted_model$L)], fitted_model$L)
  percent_cells<-percent_cells %>%
    melt(variable.factor = FALSE, varnames=c("gene", "topic"), value.name = 'percent_cells')
  reweighted_gene_score<-reweighted_f %>%
    melt(variable.factor = FALSE, varnames=c("gene", "topic"), value.name = 'reweighted_gene_score')

  gene_scores <- fitted_model$F %>%
    as.data.frame %>%
    rownames_to_column('gene') %>%
    pivot_longer(!gene, names_to='topic', values_to='gene_score') %>%
    dplyr::select(gene, topic, gene_score)
  # create df of gene name, topic, gene id, reweighted_gene_score & original gene score
  de_gene<-reweighted_gene_score %>%
    dplyr::left_join(gene_scores, by = c("gene", "topic")) %>% 
    dplyr::left_join(percent_cells, by = c("gene", "topic"))
  de_gene<-de_gene %>% dplyr::mutate(id=GeneIdMapping(organism=organism)$ids[gene])
  return(de_gene)
}

###############################################################
################## DEs heatmap visualization ##################
###############################################################
# count_mat - t(obj[['RNA']]@counts. The transposed matrix with the expression of genes across cells. Must have - rows = cells, cols = genes.
# l_mat - t(fitted_model$L).The transformed L matrix of the topic modeling. Must have rows = topics, cols = cells.
# genes_of_interest - vector of genes, for example - the concatenation of the DEs of all topics.
# Performs matrix multiplication - l_mat * count_mat = matrix with topics as rows and genes as columns.
# Each cell in that matrix is the estimated expression of that gene in that topic across all cells in the data.
heatmap_topic_gene_expression<-function(count_mat, l_mat, genes_of_interest, filename=NULL,  gaps_row=NULL, ...){
  preprocessed_count_mat<-count_mat[, genes_of_interest]
  genes_topics_mat<-l_mat %*% preprocessed_count_mat
  genes_topics_mat<-sweep(genes_topics_mat, 1, colSums(t(l_mat)), "/")
  colnames(genes_topics_mat)<-colnames(preprocessed_count_mat)
  title<-"Scaled Average Gene Expression Per Topic"
  pheatmap::pheatmap(t(genes_topics_mat), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers = FALSE,
            gaps_row=gaps_row, number_format="%.4f", scale = "row", show_rownames=FALSE, main=title, filename = filename, ...)
  return(genes_topics_mat)
}

###############################################################
#################### DEs ROC visualization ####################
###############################################################
# In case of ground truth, for example - having matching cluster-topic pairs and 
# knowing the cluster DEs, a ROC can be plotted.
# all_genes_names - the universe, the whole pool of genes - rownames(obj[["RNA"]]@count).
# topic_des - a df with all topic DEs (positive & negative, across all topics).
# cluster_des - a df with all cluster DEs (positive & negative, across all clusters).
# cur_topic - the topic we are interested in.
# cur_cluster - the cluster we are interested in.
roc_cluster_topic_des<-function(all_genes_names, topic_des, cluster_des, cur_topic, cur_cluster){
  cluster_des<-cluster_des %>% filter(cluster==cur_cluster & avg_log2FC > 0)
  topic_des<-topic_des %>% filter(topic==cur_topic)
  rownames(topic_des)<-topic_des$gene
  topic_des_positive<-topic_des %>% filter(reweighted_gene_score>0)
  df<-data.frame(gene=all_genes_names, label=FALSE, prediction=-100)
  rownames(df)<-df$gene
  df[df$gene %in% cluster_des$gene, ]$label<-TRUE
  df$prediction<-topic_des[rownames(df), ]$reweighted_gene_score
  second_min_val<-sort(unique(df$prediction))[2]
  df[is.na(df$prediction) == TRUE, ]$prediction<-second_min_val
  if (sum(df$prediction == -Inf) > 0) {
    df[df$prediction == -Inf, ]$prediction<-second_min_val  
  }

  roc_curve<-roc(df$label, df$prediction)
  # points_coord<-return_roc_thresholds(roc_curve, topic_des_positive)
  optimal_point<-coords(roc_curve, "best")
  title<-paste("Cluster=", cur_cluster, "Topic=", cur_topic)
  top_n<-topic_des_positive %>% filter(reweighted_gene_score > as.numeric(optimal_point[1])) %>% nrow()
  sub_title<-paste("THRESHOLD=", optimal_point[1], ", TOP_N=", top_n, ", FPR/1-SPECIFICITY=", (1-optimal_point[2]), ", SENSITIVITY/TPR=", optimal_point[3])
  par(pty="s")
  plot(roc_curve, col=1, lwd=4, legacy.axes=TRUE, print.thres.pch=16, print.thres.col="red", print.thres.cex=0.6, srt=110)
  mtext(title, side=3, line=1)
  mtext(sub_title, side=3, line=0.2)
}


jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <-  length(a) + length(b) - intersection
  return (intersection/union)
}


jaccard_cluster_topics_des <- function (cluster_des, topics_des, max_genes=nrow(topics_des)) {
  cluster_des <- cluster_des %>% filter(avg_log2FC>0) %>% pull(gene)
  topics_des <- topics_des %>%
    filter(reweighted_gene_score>0) %>%
    arrange(desc(reweighted_gene_score)) %>%
    slice_max(order_by = reweighted_gene_score, n=max_genes)

  jaccards <- c()
  for(i in 1:max_genes) {
    jaccards[i] <- jaccard(topics_des %>%
                             slice_max(order_by = reweighted_gene_score, n=i) %>%
                             pull(gene), cluster_des)
  }
  max_jaccard_index <- min(which(jaccards==max(jaccards), arr.ind=T))
  jaccard_graph <- ggplot(as.data.frame(cbind(jaccards=jaccards,
                             topn=seq(1:max_genes),
                             reweighted_gene_score=log(topics_des$reweighted_gene_score),
                             z_score_log=topics_des$z_score_log)),
         aes(x=topn,y=jaccards,color=z_score_log))+
    geom_point() + scale_color_viridis_c(option='turbo')

  score <- log(topics_des$reweighted_gene_score)
  z_score_log <- (score-mean(score))/sd(score)
  message(max(jaccards), " ", max_jaccard_index)
  print(topics_des$z_score_log[max_jaccard_index])
  print(z_score_log[max_jaccard_index])
  histogram_graph <- ggplot(data.frame(z_score=z_score_log), aes(x = z_score)) +
    geom_histogram()  +
    geom_vline(xintercept=c(topics_des$z_score_log[max_jaccard_index],
                            z_score_log[max_jaccard_index]), 
               colour=c('red', 'blue'))
    # geom_vline(mapping = aes(xintercept=,colour='blue'))
  return(ggarrange(jaccard_graph, histogram_graph, ncol = 2))
}
