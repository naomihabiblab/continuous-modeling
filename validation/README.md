# Validation

Validation of program signatures and pathway scores across external snRNA-seq datasets (SeaAD, PsychAD, CUIMC2, MIT). All scripts assume the **project root** as the working directory.

---

## Validation datasets

| Dataset   | Description | Cell types | Notebook | R script |
|----------|-------------|------------|----------|----------|
| **SeaAD**   | SEA-AD (Seattle Alzheimer’s Disease) | inhibitory, excitatory, astrocytes, microglia, oligo, opcs | Yes | `sead/seaAD.R` |
| **PsychAD** | PsychAD / MSSM cohort        | inhibitory, excitatory, astrocytes, microglia, oligo, opcs | Yes | `psychad/psychAD.R` |
| **CUIMC2**  | Columbia single-cell cohort      | inhibitory, excitatory, astrocytes, microglia | No  | `cuimc2.R` |
| **MIT**     | MIT single-cell atlas (resilience) | inhibitory, astrocytes, microglia, oligo, opcs + excitatory (from CSV) | For excitatory neurons (`mit/MIT_gene_signature.ipynb` for excitatory) | `mit.R` |

- **SeaAD** and **PsychAD**: Signature/pathway scores are computed in Python notebooks per cell type; the R scripts consume the notebook outputs (`pathways_and_signatures.csv.gz`, `donor_metadata.csv`) and write to HDF5, then run trait association and community analysis.
- **CUIMC2**: Uses pre-built Seurat objects (`.rds`) per cell type; R computes signatures, pathways, coherence, and trait association.
- **MIT**: Non-excitatory cell types use RDS objects; excitatory uses a precomputed CSV from the MIT notebook (because of memory limitation).

---

## Shared functions

- **`validation.functions.R`** (R): Programs name conversion, `get_genes_by_topic`, `get_pathway_signatures`, coherence helpers, HDF5 read/write (`h5_save_mean_scores`, `h5_save_correlation`, etc.), trait association and program–environment correlation. Used by all validation R scripts.
- **`validation_functions.py`** (Python): Python counterpart for coherence and related utilities; used by notebooks or other Python pipelines as needed.

Both assume project root and paths set in `R/utils.R` and `R/data_500_utils.R` (e.g. `get_result_dir()`, `get_shared_dir()`, `TOPICS_DATA_H5_PATH`).

---

## How to run validation

### 1. Datasets with Jupyter notebooks (SeaAD, PsychAD)

For **SeaAD** and **PsychAD**, run the notebook **first** for each cell type to generate per–cell-type outputs under `Plots/<project>/<cell_type>/` (e.g. `pathways_and_signatures.csv.gz`, `donor_metadata.csv`). Then run the R script.

#### Option A: Papermill (from project root)

Run the notebook once per cell type and write output to a new notebook (or overwrite):

**SeaAD:**

```bash
papermill validation/sead/SEA_AD_gene_signature.ipynb output_seaad_inhibitory.ipynb -p cell_type inhibitory
papermill validation/sead/SEA_AD_gene_signature.ipynb output_seaad_excitatory.ipynb -p cell_type excitatory
# Repeat for: astrocytes, microglia, oligo, opcs
```

**PsychAD:**

```bash
papermill validation/psychad/PSYCH_AD_gene_signature.ipynb output_psychad_inhibitory.ipynb -p cell_type inhibitory
papermill validation/psychad/PSYCH_AD_gene_signature.ipynb output_psychad_excitatory.ipynb -p cell_type excitatory
# Repeat for: astrocytes, microglia, oligo, opcs
```

Valid `cell_type` values: `inhibitory`, `excitatory`, `microglia`, `astrocytes`, `oligo`, `opcs`. You can change `project` or `root_plots_dir` via `-p` if needed (e.g. `-p root_plots_dir ./OtherPlots`).

#### Option B: Manual (run inside the notebook)

1. Open `validation/sead/SEA_AD_gene_signature.ipynb` or `validation/psychad/PSYCH_AD_gene_signature.ipynb`.
2. In the **“Cell type and run options”** section, set `cell_type` in the next code cell to one of: `"inhibitory"`, `"excitatory"`, `"microglia"`, `"astrocytes"`, `"oligo"`, `"opcs"`.
3. Run all cells. Repeat after changing `cell_type` for each cell type you need.

Ensure outputs are written under the same layout the R scripts expect (e.g. `Plots/seaad/<cell_type>/` or `Plots/psychAD/<cell_type>/`).

### 2. R scripts (all validation datasets)

From the **project root**, run the R script for each dataset. Each script loops over the relevant cell types and runs trait association and community analysis after the loop.

```bash
# From project root
Rscript validation/sead/seaAD.R
Rscript validation/psychad/psychAD.R
Rscript validation/cuimc2.R
Rscript validation/mit.R
```

- **SeaAD** and **PsychAD**: Require that the notebooks have already been run for all desired cell types so that `pathways_and_signatures.csv.gz` and `donor_metadata.csv` exist in each `Plots/<project>/<cell_type>/`.
- **CUIMC2**: Requires Seurat objects at the paths defined in the script (e.g. under `get_shared_dir()` / `CUIMC2_OBJECT_DIR`). Skips cell types whose object file is missing.
- **MIT**: Requires RDS objects for inhibitory, astrocytes, microglia, oligo, opcs; excitatory uses a CSV at `Plots/mit/excitatory/pathways_and_signatures.csv` (e.g. from `mit/MIT_gene_signature.ipynb` if used).


## Outputs

- **Per dataset/cell type**: Signature and pathway scores, mean scores by sample/donor, correlation matrices, and (where applicable) coherence results are written into the project HDF5 file (see `TOPICS_DATA_H5_PATH` and `validation.functions.R`).
- **Plots**: PDFs and other figures (trait associations, heatmaps, etc.) are written under `Plots/<project>/<cell_type>/` or the path returned by `get_result_dir()`.

Data paths (e.g. shared data, result dir, HDF5 path) are controlled by `R/utils.R` and `R/data_500_utils.R`; configure them there for your environment.
