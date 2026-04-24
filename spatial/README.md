# Spatial Transcriptomics Analysis

Analysis of spatial transcriptomics data from AD tissue samples using Stereo-seq technology.

## Overview

This directory contains scripts for analyzing spatial co-localization patterns of neuronal stress programs, astrocyte activation, and microglial responses in Alzheimer's disease brain tissue.

## Data Source

- **Technology**: Stereo-seq spatial transcriptomics
- **GEO Accession**: GSE269906
- **Samples**: 12 human brain samples (6 AD, 6 control)
  - **Control**: D02175A4, D02175A6, B01809A3, B01809A4, B01806B5, B01806B6
  - **Moderate AD**: B02008D2, B02009F6
  - **Advanced AD**: B01809C2, C02248B5
  - **Severe AD**: A02092E1, B02008C6

## Workflow

### 1. Data Processing (`data_processing_paper.ipynb`)

**Purpose**: Convert raw Stereo-seq data to annotated AnnData objects

**Steps**:
1. Read raw GEF files from GEO
2. Compute QC metrics
3. Integrate samples using Harmony
4. Cell type annotation using marker genes
5. Save annotated object

**Outputs**:
- Individual sample H5AD files: `{H5AD_OUTPUT_DIR}/{sample_id}.h5ad`
- Integrated object: `{H5AD_OUTPUT_DIR}/all_obj_annotated.h5ad`

### 2. Topic Score Computation (`script_for_paper.ipynb`)

**Purpose**: Compute program and pathway signatures scores for all cells

**Steps**:
1. Load integrated annotated object
2. Extract gene lists from 500 project programs (`topics_data.h5`)
3. Process each cell type separately (scale, variable genes)
4. Compute programs scores using `sc.tl.score_genes()`
5. Binarize scores using quantile thresholds (75%, 80%, 85%)
6. Compute pathway signature scores
7. Identify bimodal distributions using Gaussian Mixture Models
8. Save final object with all scores

**Pathway Signatures**:
- Metallothioneins (MT genes)
- Cell cycle
- DNA repair
- Mitochondrial stress
- Synaptic function
- Oxidative stress response

**Outputs**:
- Cell type-specific objects: `{CELL_TYPES_DIR}/{cell_type}.h5ad`
- Final spatial object: `{FINAL_OBJ_DIR}/spatial_obj_for_paper.h5ad`

### 3. Visualization and Analysis (`spatial_figures.R`)

**Purpose**: Generate publication figures showing spatial patterns

**Analyses**:

#### A. Enrichment Analysis
- Compute co-localization enrichment scores
- Test significance using permutation tests
- Visualize as boxplots with significance stars

#### B. Spatial Visualizations
- Full slice plots showing topic expression patterns
- Neighborhood plots highlighting co-localization
- Cell-level zoom plots with circular neighborhoods

#### C. Cell Type Distribution
- Stacked bar plots of cell type proportions per sample
- Stratified by disease severity

#### D. Score Distributions
- Density plots showing score distributions
- Quantile thresholds visualization
- Bimodal distribution fitting

**Key Functions**:
- `plot_pairs_boxplot()`: Enrichment boxplots with significance
- `plot_topic_xy_all_samples()`: Spatial topic expression
- `plot_cell_neighborhood()`: Zoomed neighborhood view
- `plot_full_slice()`: Full tissue slice with multiple markers

**Outputs**: All plots saved to `{SPATIAL_PLOTS_DIR}/`

## File Descriptions

## Setup and Usage

The configuration automatically detects if running on cluster or local mount.

### Environment Requirements

**Python**:
requirements.txt

**R**:
```r
# Core packages
reticulate  # Python integration
ggplot2, dplyr, tidyverse
ggforce     # geom_circle
ggrastr     # Rasterized plots
withr       # PDF output
```

### Running the Analysis

#### On Cluster

## Integration with main project

### Data Dependencies

- **Topics data**: `topics_data.h5`
  - Used to extract gene lists for each program
  - Cell types: inhibitory, excitatory (cux2+/cux2-), astrocytes, microglia, oligodendrocytes

- **Path utilities**: `utils.R`
  - `get_result_dir()`: Auto-detect cluster vs local for plots
  - `get_shared_dir()`: Auto-detect cluster vs local for shared data

### Outputs

- Plots saved to same structure as main project: `{PLOTS_DIR}/spatial/`

