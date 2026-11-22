# Cell-type specific program dynamics along Alzheimer's disease progression reveal cell-cycle reentry in resilient neurons

Code for the paper. All scripts assume the **project root** as the working directory. Shared R utilities and analysis helpers live in **`R/`** and are sourced as `source("R/...")`.

The code is divided into six parts:

1. **Cell-type analysis**: NMF-based co-expression gene program analysis of the snRNA-seq data for each cell type, program signatures (differentially expressed genes) and pathway enrichment.
   - R files under `cell_type_analysis/`.
   - First run `{cell.type}TopicRunner.R` to fit the programs for the cell type.
   - Then run `{cell.type}_topics.R` to compute the program signatures and pathway enrichment and save the results.
   - `R/topic_loading_permutation.R`: within-cell permutation null for the correlation between two programs' loadings across cells (used for the cell-cycle reentry vs stress/apoptosis programs; usage in the file header).

2. **Program space analysis**: trait associations, BEYOND and multi-cell-type analysis, in `topic_space/` (see [topic_space/README.md](topic_space/README.md)).
   - Main script: `topic_space/topic_space.R`
   - **BEYOND**: applies the BEYOND methodology ([*Green et al*, Nature 2024](https://doi.org/10.1038/s41586-024-07871-6)) to our program atlas.

3. **Validation**: signature scores and trait associations in the validation datasets (PsychAD, CUIMC2, MIT, SEA-AD), in `validation/` (see its README).

4. **Spatial transcriptomics**: see [spatial/README.md](spatial/README.md).

5. **Manuscript code**: the code generating all figures, Extended Data figures and Supplementary Tables in the manuscript.
   - `paper/figures.R`: one block per figure panel.
   - `paper/create.supplementary.tables.R`, with `paper/supplementary_tables_neuronal.R` (the neuronal program sheets of Supplementary Tables 2 and 3).
   - `spatial/spatial_figures.R`: the spatial panels.

6. **Source Data**: the Source Data files (one `Source_Data_<Fig>.xlsx` per figure) are built from the same code that draws the figures, so the data match the plots.
   - `R/source_data.R`: `sd_*` helpers called inside the figure blocks; each call writes the data behind a panel to `{result_dir}/source_data/staging/`. Set `SD_ENABLED=FALSE` to turn them off.
   - `paper/run_figure_blocks.R`: runs selected blocks of `paper/figures.R` (`--list` to see them, `--blocks` to choose). Plots are written to `{result_dir}/source_data/figures_rerun/`, not over the original figures.
   - `paper/run_spatial_source_data.R`: the same for the spatial panels.
   - `paper/source_data_compile.R`: turns the staging tables into the final xlsx files (plus a zip for tables too large for Excel).
   - `paper/source_data_fig3b.R`: the Fig. 3b network tables.

<details>
<summary>Running the code</summary>

To run the code you may need to change the location of the data object directories, mainly:

- `SEURAT_OBJECT_DATA_PATH`
- `SEURAT_SUBSET_OBJECT_PATH`
- `SEA_AD_ROOT_PATH`
- `PROTEOMIC_PATH`
- the result directory in the function `get_result_dir`
</details>

## Citation

If you use this code or the associated results, please cite:

> Meir R, Schwartz G et al. *Cell-type specific program dynamics along Alzheimer's disease progression reveal cell-cycle reentry in resilient neurons*.
> Preprint (earlier title: *Early neuronal reprogramming and cell cycle reentry shape Alzheimer's disease progression*): bioRxiv 2025.06.04.653670; doi: [10.1101/2025.06.04.653670](https://www.biorxiv.org/content/10.1101/2025.06.04.653670v1).

**Interactive app:** explore program dynamics and associations in the [ADynamics](https://habib-lab.shinyapps.io/ADynamics) Shiny app.

## Data and requirements

<details>
<summary>Data accessibility</summary>

All snRNA-seq data used in our study are accessible via [Synapse (syn53366818)](https://www.synapse.org/#!Synapse:syn53366818), including the raw reads, library count matrices and the processed cell atlas as cell-type Seurat objects.
</details>

<details>
<summary>Key intermediate data files: H5 and H5AD</summary>

- **`topics_data.h5`** (path: `{result_dir}/topics_data.h5`, set via `get_result_dir()`)
  - Main HDF5 store for **cell-type-level** program results. Written by the cell-type analysis scripts and by `topic_space.R` and the validation scripts.
  - Layout: one group per cell type (e.g. `neuronal/inhibitory`, `neuronal/excitatory`, `glia/astrocytes`, …). Within each group: `/de` (differential genes), `/pathways`, `/L`, `/F`, `/program_signature` (unique signature genes), abundance and UMAP datasets, and one subgroup per validation project (e.g. `cuimc2`, `psychAD`).

- **`500_topics_all.h5ad`** (path: `{result_dir}/topic_space/500_topics_all.h5ad`)
  - AnnData file for the **integrated program space**: donor-level program abundances, embeddings (e.g. PHATE, UMAP), BEYOND trajectories and dynamics, and metadata. Produced by `topic_space.R` and used by `paper/figures.R` and `paper/create.supplementary.tables.R`.
</details>

<details>
<summary>Software requirements</summary>

## OS requirements
Tested on:
- Gentoo/Linux 2.7
- macOS Monterey (12.4)

## Dependencies
- Main analysis: R (4.0.5), dependencies listed in `renv.lock`
- Validation signature scores (notebooks in `validation/*/`): Python (3.12)
- Spatial: Python (3.12), dependencies listed in `spatial/requirements.txt`

To create or update a snapshot of the installed R dependencies (e.g. after installing new packages), run `renv::snapshot()` in R. This updates `renv.lock` to match the current project library.
</details>

<details>
<summary>Hardware requirements</summary>

## Memory
- fastTopics fitting: some fits require 200–400 GB
- SEA-AD signatures: about 200 GB
- Program analysis: 8 GB should suffice

## CPU
- Most of the code needs only one core.
- fastTopics fitting can take a few days. More cores speed up the fit substantially but may need more memory. The number of cores is set with the `nc` parameter.
</details>
