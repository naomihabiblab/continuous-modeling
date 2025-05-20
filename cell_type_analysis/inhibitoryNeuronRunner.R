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

# 264gb for the DE
# 4 cpu

output_dir <- file.path(get_result_dir(), "inhibitory_topics")
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

k <- Sys.getenv(x = 'SLURM_ARRAY_TASK_ID')

if (is.null(k) | k == '') {
  stop("No k provided ")
}
k <- as.numeric(k)

output_name_prefix <- "all_inhibitory"
print(output_name_prefix)

input_path <- file.path(SEURAT_OBJECT_DATA_PATH, "inhibitory.h5Seurat")
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

# TODO: add flag to run only the subtype fit
# inh_classes <- unique(obj@meta.data$annotation)
#
# for (inh_class in inh_classes) {
#   obj_subset <- subset(obj, annotation == inh_class)
#   count.transpose <- t(obj_subset[['SCT']]@counts)
#   count.transpose <- count.transpose[, colSums(count.transpose) > 0]
#
#   print(paste0("Starting on ",  inh_class))
#   output_name_prefix <- sprintf("all_%s_inhibitory", inh_class)
#
#   print(paste("Fiting topic model k=", k))
#   tic()
#   topic_fit <- fit_topic_model(count.transpose, k = k,
#                                control.init = list(nc = 4),
#                                control.main = list(numiter = 4, nc = 4),
#                                control.refine = list(numiter = 4,extrapolate = TRUE, nc = 4))
#   toc()
#
#   print("Saving")
#   saveRDS(topic_fit, sprintf("%s_topic_fit.%s.RDS", output_name_prefix, k))
#   print("And done!")
#   print(paste("Loading fit", k))
#   topic_fit <- readRDS(sprintf("%s_topic_fit.%s.RDS", output_name_prefix, k))
#
#   print("Running de_analysis")
#   dfa <- de_analysis(fit = topic_fit, X = count.transpose, control = list(nc = 2))
#   saveRDS(dfa, file = sprintf("%s_dfa.%s.RDS", output_name_prefix, k))
#
#   print("Running de_analysis vsnull")
#   dfa_vsnull <- de_analysis(fit = topic_fit, lfc.stat = 'vsnull',
#                             X = count.transpose,
#                             control = list(nc = 2))
#   saveRDS(dfa_vsnull, file = sprintf("%s_dfa_vsnull.%s.RDS", output_name_prefix, k))
#
#   print(paste("All done", inh_class))
# }
#
# count.transpose <- t(obj[['SCT']]@counts)
# count.transpose <- count.transpose[, colSums(count.transpose) > 0]

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
                             control.init = list(nc = 8),
                             control.main = list(numiter = 4, nc = 8),
                             control.refine = list(numiter = 4,extrapolate = TRUE, nc = 8))
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


# subset fit
output_name_prefix <- "subset_inhibitory"
print(output_name_prefix)

input_path <- file.path(SEURAT_SUBSET_OBJECT_PATH, "inhibitory.subset.rds")
message("Loading data from ", input_path)
obj <- readRDS(input_path)

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
                             control.init = list(nc = 8),
                             control.main = list(numiter = 4, nc = 8),
                             control.refine = list(numiter = 4,extrapolate = TRUE, nc = 8))
toc()

print("Saving")
saveRDS(topic_fit, sprintf("%s_topic_fit.%s.RDS", output_name_prefix, k))

