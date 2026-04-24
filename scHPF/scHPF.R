library(ComplexHeatmap)

source("R/utils.R")
source("R/data_500_utils.R")
source("R/topic_space_functions.R")
source("R/landscape.analysis.funs.R")
source("R/figures.funs.R")

version_dir <- "v1.1.objects"
cell_type <- "inhibitory"
cell_type_file_name <- "inhibitory"

# cell_type <- "inhibitory_pvalb"
# cell_type_file_name <- "inhibitory_pvalb"

cell_type_500_name <- case_when(cell_type == "inhibitory_pvalb" ~ "inhibitory",
                           cell_type == "excitatory_l2_3" ~ "cux2p",
                           cell_type == "excitatory_l4_6" ~ "cux2m",
                           TRUE ~ cell_type)


K <- case_when(cell_type == "astrocytes" ~ 10,
               cell_type == "inhibitory" ~ 17)

if (cell_type == "astrocytes") {
    USE_FLOAT_32 <- FALSE
} else {
    USE_FLOAT_32 <- TRUE
}
float_suffix <- ifelse(USE_FLOAT_32, "float32", "")

schpf_results_dir <- file.path(get_result_dir(), "scHPF", ifelse(USE_FLOAT_32, paste0(cell_type_file_name, "_float32"), cell_type_file_name))
dir.create(schpf_results_dir, showWarnings = FALSE, recursive = TRUE)

loom_file <- file.path(get_result_dir(), "../scHPF/data", glue::glue("{cell_type_file_name}.loom"))

# Save as loom
CONVERT_TO_LOOM <- FALSE
if (CONVERT_TO_LOOM) {
    input_path <- file.path(RAW_DATA_PATH, version_dir, glue::glue("{cell_type_file_name}.h5Seurat"))
    message("Loading ", input_path)

    obj <- SeuratDisk::LoadH5Seurat(input_path)

    DefaultAssay(obj) <- "RNA"

    SeuratDisk::SaveLoom(obj, loom_file, overwrite=TRUE)
}

# run scHPF
# scHPF prep -i data/astrocytes.loom -o astrocytes -m 10
# schPF train -i astrocytes/filtered.mtx -o astrocytes/ -p astrocytes_schpf_K9 -k 9 -t 4
# scHPF score -m astrocytes/astrocytes_schpf_K9.scHPF_K9_b0_4trials.joblib -o astrocytes -p astrocytes_schpf_K9  -g astrocytes/genes.txt

# python ~/submit_job.py --export PYTHONUNBUFFERED=1  --array 1 --out logs/  -c 2  --memory=320g -J schpf_inhibitory  --conda_env schpf_p37 --no-r-script 'scHPF train -i inhibitory/filtered.mtx -o inhibitory/ -p inhibitory_schpffloat32_k17_t$SLURM_ARRAY_TASK_ID -k 17 -t 1 --float32'
# python bin/scHPF select-best -o tmp/ -p inhibitory_combines_schpffloat32_k17 --save-all -i /Volumes/habib-lab/roi.meir/500/scHPF/inhibitory/inhibitory_schpffloat32_k17_t*.scHPF_K17_b0_1trials.joblib

# Load results




fit_name <- glue::glue("{cell_type}_{float_suffix}K{K}")

schpf_data_dir <- file.path(get_result_dir(), "../scHPF/", cell_type_file_name)
# schpf_data_dir <- file.path("/Users/habib.lab/Documents/Roi/scHPF/tmp")

results_dir_k <- file.path(schpf_results_dir, K)
dir.create(results_dir_k, showWarnings = FALSE)

prefix <- glue::glue("{cell_type_file_name}_schpf{float_suffix}_K{K}")

cell_score_file <- glue::glue("{prefix}.cell_score.txt")
gene_score_file <- glue::glue("{prefix}.gene_score.txt")
maximum_overlaps_file <- glue::glue("{prefix}.maximum_overlaps.txt")
mean_cellscore_fraction_file <- glue::glue("{prefix}.mean_cellscore_fraction.txt")

cell_score <- read.table(file.path(schpf_data_dir, cell_score_file))
colnames(cell_score) <- paste0("P", 1:K)

cell_ids <- h5read(loom_file, "col_attrs/CellID")
projID <- h5read(loom_file, "col_attrs/projid")
rownames(cell_score) <- cell_ids


gene_score <- read.table(file.path(schpf_data_dir, gene_score_file))
colnames(gene_score) <- paste0("P", 1:K)

rownames(gene_score) <- read.table(file.path(get_result_dir(), "../scHPF/", cell_type_file_name, "genes.txt"), header = FALSE)$V1


run_metadata <- jsonlite::read_json(file.path(schpf_data_dir, glue::glue("{prefix}.score_commandline_args.json")))

schpf_group <- paste0(CELL_TYPE_TO_GROUP[cell_type_500_name], "/schpf")

h5createGroup(TOPICS_DATA_H5_PATH, schpf_group)
h5write(gene_score %>% rownames_to_column("gene"), TOPICS_DATA_H5_PATH, paste0(schpf_group, "/gene_score"))
h5write(cell_score %>% rownames_to_column("cell"), TOPICS_DATA_H5_PATH, paste0(schpf_group, "/cell_score"))
h5write(
  obj  = as.character(jsonlite::toJSON(run_metadata, auto_unbox = TRUE)),
  file = TOPICS_DATA_H5_PATH,
  name = paste0(schpf_group, "/score_args")
)


# Load L matrix (cell-topic matrix) from the topic h5 file
L_matrix <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell_type_500_name], "/L")) %>% column_to_rownames("cell")
L_matrix <- L_matrix[rownames(cell_score), ]
# L_matrix: rows = cells, cols = topics

# Compare program scores using correlation (Spearman and Pearson)

for (corr_method in c("spearman", "pearson")) {
    c <- correlate_dataframes(L_matrix, cell_score, method = corr_method)

    add_significance <- TRUE
    if (add_significance) {
        corr.sig <- c$sig
        corr.sig[corr.sig != ''] <- '*'
        cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)
    } else {
        cell_fun <- NULL
    }

    similarity <- c$corr
    similarity[similarity < 0] <- 0
    order <- slanter::slanted_orders(similarity)

    legend <- Legend(title = "adj.pval", pch = c("*"), type = "points", labels = c("<0.05"), title_gp = gpar(fontface='plain'))
    correlation.legend <- Legend(title = glue::glue("correlation {corr_method}"), col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

    withr::with_pdf(file.path(results_dir_k, glue::glue("scHPF_K{K}_correlation_{corr_method}.pdf")), {
        ComplexHeatmap::Heatmap(c$corr,
                                heatmap_legend_param = list(title="scHPF programs correlation", direction="vertical"),
                                col = cyan2orange.correlation,
                                clustering_distance_rows = "spearman",
                                clustering_distance_columns = "spearman",
                                show_row_dend = FALSE,
                                show_column_dend = FALSE,
                                row_title = "fastTopics programs",
                                column_title = "scHPF programs",
                                row_order = order$rows,
                                column_order = order$cols,
                                column_names_rot = 90,
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_title_side = 'bottom',
                                show_heatmap_legend =FALSE,
                                cell_fun = cell_fun) %>% draw(annotation_legend_list=list(legend, correlation.legend),
                                                            column_title=glue::glue("scHPF vs fastTopics {fit_name} cell score correlation"))
    })
}




# Compute correlation between gene loadings (F matrix) and the gene weights from scHPF

# Read F matrix from fastTopics
F_matrix <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell_type_500_name], "/F")) %>% column_to_rownames("gene")
# F_matrix: rows = genes, cols = topics

# Subset to genes present in both F_matrix and scHPF gene weights
common_genes <- intersect(rownames(F_matrix), rownames(gene_score))
F_matrix <- F_matrix[common_genes, , drop=FALSE]


for (corr_method in c("spearman", "pearson")) {
    c <- correlate_dataframes(F_matrix, gene_score, method = corr_method)

    add_significance <- TRUE
    if (add_significance) {
        corr.sig <- c$sig
        corr.sig[corr.sig != ''] <- '*'
        cell_fun <- function(j, i, x, y, w, h, fill) grid.text(corr.sig[i,j], x, y)
    } else {
        cell_fun <- NULL
    }

    similarity <- c$corr
    similarity[similarity < 0] <- 0
    order <- slanter::slanted_orders(similarity)

    legend <- Legend(title = "adj.pval", pch = c("*"), type = "points", labels = c("<0.05"), title_gp = gpar(fontface='plain'))
    correlation.legend <- Legend(title = glue::glue("correlation {corr_method}"), col_fun = cyan2orange.correlation, title_gp = gpar(fontface='plain'))

    withr::with_pdf(file.path(results_dir_k, glue::glue("scHPF_K{K}_gene_correlation_{corr_method}.pdf")), {
        ComplexHeatmap::Heatmap(c$corr,
                                heatmap_legend_param = list(title="scHPF gene correlation", direction="vertical"),
                                col = cyan2orange.correlation,
                                clustering_distance_rows = "spearman",
                                clustering_distance_columns = "spearman",
                                show_row_dend = FALSE,
                                show_column_dend = FALSE,
                                row_title = "fastTopics gene loadings",
                                column_title = "scHPF gene weights",
                                row_order = order$rows,
                                column_order = order$cols,
                                column_names_rot = 90,
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_title_side = 'bottom',
                                show_heatmap_legend =FALSE,
                                cell_fun = cell_fun) %>% draw(annotation_legend_list=list(legend, correlation.legend),
                                                            column_title=glue::glue("scHPF vs fastTopics {fit_name} gene correlation"))
    })
}



overlap_index <- function(set1, set2) {
    length(intersect(set1, set2)) / min(length(set1), length(set2))
}

# Read DE genes from H5 (assume path and group are available)
# You may need to adjust the path/group as appropriate for your data
de_genes <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell_type_500_name], "/de"))


# Compute overlap for top N genes
top_ns <- c(100, 150)
overlap_results <- data.frame(N=integer(), Overlap=numeric())



for (N in top_ns) {
    # Get top N genes for each program 
    gene_score_top_n <- sapply(colnames(gene_score), function(col) {
        gene_score %>%
            slice_max(order_by = !!sym(col), n = N) %>%
            rownames() %>% as_tibble()
    }, simplify = FALSE)

    des_nested <- de_genes %>%
        group_by(topic) %>%
        dplyr::select(gene, topic) %>%
        group_nest()

    mat <- outer(gene_score_top_n,
                 des_nested %>% pull(data),
                 Vectorize(tibble_overlap_coefficient))
    colnames(mat) <- des_nested %>% pull(topic)
    

    withr::with_pdf(file.path(results_dir_k, glue::glue("scHPF_K{K}_gene_overlap_DE_{N}.pdf")), {
        ComplexHeatmap::Heatmap(mat,
                                col = cyan2orange.correlation,
                                show_row_dend = FALSE,
                                show_column_dend = FALSE) %>% draw
    })

}

schpf_abundance <- cbind(cell_score, projid=projID) %>% group_by(projid) %>%
 summarise(across(everything(), mean)) %>% column_to_rownames('projid')


programs_abundance_df <- h5read(TOPICS_DATA_H5_PATH, paste0(CELL_TYPE_TO_GROUP[cell_type_500_name], "/abundance.projid")) 

programs_abundance_wide_df <- programs_abundance_df %>%
    mutate(variable=glue::glue("k{variable}")) %>%
    tidyr::pivot_wider(id_cols = 'projid', names_from = 'variable', values_from = avg) %>%
    arrange(match(projid, rownames(schpf_abundance))) %>%
    filter(projid %in% rownames(schpf_abundance)) %>%
    column_to_rownames('projid')

withr::with_pdf(file.path(results_dir_k, glue::glue("scHPF_K{K}_abundance_vs_mean_score_correlation.pdf")), {
    ComplexHeatmap::Heatmap(
        cor(programs_abundance_wide_df, schpf_abundance),
        name = "correlation",
        col = cyan2orange.correlation,
        clustering_distance_rows = "spearman",
        clustering_distance_columns = "spearman",
        row_title = "schPF programs",
        column_title = "fastTopics programs",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
    ) %>% draw(column_title = glue::glue("{cell_type} - Correlation of abundance vs mean score"))
})



data <- read_h5ad(AD_500_ANDATA_FILE)

traits_columns <- c("sqrt.amyloid", "sqrt.amyloid_mf", "sqrt.tangles", "sqrt.tangles_mf", "cogng_demog_slope")
traits <- data$obsm$meta.data[,c(traits_columns)] %>% mutate(across(all_of(traits_columns), ~as.numeric(as.character(.))))

controls <- data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")])

trait.analysis <- associate.traits(traits,
                                    schpf_abundance[rownames(data),],
                                    controls) %>% dplyr::mutate("topic"="covariate")


withr::with_pdf(file.path(results_dir_k, "TraitAssocationHeatmap.pdf"), {
  plot.trait.associations.heatmap(trait.analysis, params = traits_columns, row.by='covariate')
})

