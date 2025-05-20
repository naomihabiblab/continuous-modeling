#!/usr/bin/env Rscript
library(dplyr)
library(patchwork)
library(ggplot2)
library(rhdf5)
library(future)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(fastTopics)
library(tictoc)
library(stringr)
library(pryr)

source("utils.R")
source("cell_type_analysis/fit_topics.R")


output_dir <- file.path(get_result_dir(), "oligo_topics")
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

k <- Sys.getenv(x = 'SLURM_ARRAY_TASK_ID')

if (is.null(k) | k == '') {
  stop("No k provided ")
}
k <- as.numeric(k)

output_name_prefix <- "all_oligo"
print(output_name_prefix)

input_path <- file.path(SEURAT_OBJECT_DATA_PATH, "oligodendrocytes.h5Seurat")
message("Loading data from ", input_path)
obj <- LoadH5Seurat(input_path,
                    assays = list(SCT=c('counts'),
                                  RNA=c('counts')),
                      reductions = FALSE,
                      graphs = FALSE,
                      neighbors = FALSE,
                      tools = FALSE,
                      misc = FALSE)

print("After loading")
gc()
mem_used()

count.transpose <- t(obj[['SCT']]@counts)
count.transpose <- count.transpose[, colSums(count.transpose) > 0]

print(paste("Number of cells and genes:", dim(count.transpose)))

print("After count")
mem_used()

rm(obj)
gc()

print("After removing")
mem_used()


print(paste("Fiting topic model k=", k))
tic()
topic_fit <- fit_topics(count.transpose, k = k,
                        name=output_name_prefix,
                             control.init = list(nc = 6),
                             control.main = list(numiter = 4, nc = 6),
                             control.refine = list(numiter = 4,extrapolate = TRUE, nc = 6))
toc()

print("Saving")
saveRDS(topic_fit, sprintf("%s_topic_fit.%s.RDS", output_name_prefix, k))


print("Calcualting log likelihood...")
loglik <- loglik_multinom_topic_model(count.transpose, topic_fit)
saveRDS(loglik, sprintf("%s_topic_fit_loglik.%s.RDS", output_name_prefix, k))

rm(loglik)

print("And done!")

print(paste("Loading fit", k))
topic_fit <- readRDS(sprintf("%s_topic_fit.%s.RDS", output_name_prefix, k))

print("Running de_analysis vsnull")
tic()
dfa_vsnull <- de_analysis(fit = topic_fit, lfc.stat = 'vsnull',
                          X = count.transpose, 
                          control = list(nc = 4))
saveRDS(dfa_vsnull, file = sprintf("%s_dfa_vsnull.%s.RDS", output_name_prefix, k))
toc()
print("All done")
