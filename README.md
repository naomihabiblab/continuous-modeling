# Early neuronal reprogramming and cell cycle reentry shape Alzheimer's disease progression paper's code

All scripts assume the **project root** as the working directory. Shared R utilities and analysis helpers live in **`R/`** and are sourced as `source("R/...")`.

The code is divided into four parts:

1. **Cell-type analysis**: snRNA-seq NMF-based co-expression gene programs analysis for each of the cell types. Program signature (differential expressed genes) and pathway enrichment analyses.
   - R files under `cell_type_analysis/`.
   - First run `{cell.type}TopicRunner.R` to fit the programs for the cell type.
   - Then run `{cell.type}_topics.R` to compute program signature (differential expressed genes), pathway enrichment analyses and save the results.

2. **Program space analysis**: Trait associations and BEYOND and multiple cell type analisys – `topic_space/` (see [topic_space/README.md](topic_space/README.md))
   - Main script: `topic_space/topic_space.R`
   - **BEYOND**: Applying BEYOND methodology ([*Green et al*, Nature 2024](https://doi.org/10.1038/s41586-024-07871-6)) steps over our program atlas.

3. **Validation**
   - Computing signature score in validation datasets (PsychAD, CUIMC2, MIT, SEA-AD) at `validation`. (see dedicated README)

4. **Spatial transcriptomics**
   - See `spatial/README.md`. 

5. **Manuscript code**: The complete code generating all figures, extended data figures and supplementary tables shown in the manuscript.
   - `paper/figures.R`
   - `paper/create.supplementary.tables.R`

<details>
<summary>Running the code</summary>
To run the code you may need to change the location of data object directories. Mainly

- `SEURAT_OBJECT_DATA_PATH`
- `SEURAT_SUBSET_OBJECT_PATH`
- `SEA_AD_ROOT_PATH`
- `PROTEOMIC_PATH`
- result dir in the function `get_result_dir`
</details>

## Citation

If you use this code or the associated results, please cite:

> Meir R, Schwartz G et. al. *Early neuronal reprogramming and cell cycle reentry shape Alzheimer's disease progression*. bioRxiv 2025.06.04.653670; doi: [10.1101/2025.06.04.653670](https://www.biorxiv.org/content/10.1101/2025.06.04.653670v1).

**Interactive app:** Explore program dynamics and associations in an interactive Shiny app: [ADynamics](https://habib-lab.shinyapps.io/ADynamics).

## DATA and requirements

<details>
<summary>Data accessibility</summary>
All snRNA-seq data used in our study is accessible via [Synapse (syn53366818)](https://www.synapse.org/#!Synapse:syn53366818) and contains the raw reads, library count matrices and processed cell atlas in the format of cell-type Seurat objects. 
</details>

<details>
<summary>Key intermidate data files: H5 and H5AD</summary>

- **`topics_data.h5`** (path: `{result_dir}/topics_data.h5`, set via `get_result_dir()`)
  - Main HDF5 store for **cell-type-level** topic results. Written by the cell-type analysis scripts and by `topic_space.R` / validation scripts.
  - Layout: one group per cell type (e.g. `neuronal/inhibitory`, `neuronal/excitatory`, `glia/astrocytes`, …). Within each group: `/de` (differential genes), `/pathways`, `/L`, `/F`, `/program_signature` (unique signature genes), abundance and UMAP datasets, and per-validation-project subproject (e.g. `cuimc2`, `psychAD`).

- **`500_topics_all.h5ad`** (path: `{result_dir}/topic_space/500_topics_all.h5ad`)
  - AnnData file for the **integrated topic space**: donor-level program abundances, embeddings (e.g. PHATE, UMAP), BEYOND trajectories and dynamics, and metadata. Produced by `topic_space.R` and used by `paper/figures.R` and `paper/create.supplementary.tables.R` for cohort description, abundance matrices, and trait associations.
</details>

<details>
<summary>Software requirements</summary>

## OS Requirements
Tested on 
- Gentoo/Linux 2.7
- macOS: Monterey (12.4)

## Dependencies:
- Main analysis: R (4.0.5), dependencies listed in `renv.lock`
- SEA_AD signature score: python (3.12), dependencies listed in `validation/sead/requirements.txt`

To create or update a snapshot of the installed R dependencies (e.g. after installing new packages), run in R: `renv::snapshot()`. This updates `renv.lock` to match the current project library.
</details>

<details>
<summary>Hardware requirements</summary>

## Memory

- fastTopics fitting: some fit requires between 200GB to 400GB
- SEA_AD signatures: around 200GB
- program analysis:  8GB should suffice 

## CPU
- For most code only one core is needed 
- fastTopics fitting: fit can take up few days, increasing the number of core will fasten the fit substantially but may require more memory.
The number of used core can be controlled using the `nc` parameter

## Dependencies:
- Main analysis: R (4.0.5), dependencies listed in `renv.lock`
- Validation signature score: python (3.12), dependencies listed in `validation/requirements.txt` 
- Spatial: python (3.12), dependencies listed in `spatial/requirements.txt` 
</details>
