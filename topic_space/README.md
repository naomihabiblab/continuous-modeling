# Topic space analysis

This directory contains the **program space** (topic space) analysis: donor-level program abundances, trait associations, and BEYOND (trajectories, communities, dynamics).

## Script

- **`topic_space.R`** — Builds the integrated topic space, runs trait associations and BEYOND, and writes outputs to the result directory and to `topics_data.h5`.

## Running

Run from the **project root** (all paths in the script assume project root):

```bash
Rscript topic_space/topic_space.R
```

Or from R:

```r
source("topic_space/topic_space.R")
```

## Inputs

- Cell-type topic scores and metadata (via `read_topic_space()` from `R/topic_space_functions.R`)
- `AD_500_ANDATA_FILE` — AnnData with donor metadata (see `R/data_500_utils.R`)
- Existing per-cell-type results in `topics_data.h5` (e.g. `/de`, `/L`, `/F`) from the cell-type analysis

## Outputs

- **`{result_dir}/topic_space/500_topics_all.h5ad`** — Integrated AnnData (donor × program abundances, embeddings, trait analysis, BEYOND results)
- **`{result_dir}/topic_space/*.pdf`** — Correlation heatmaps and other figures
- **`topics_data.h5`** — Updated with abundance stats, cell-cycle/apoptosis scores, pathway genes JSON paths, etc. (see main [README](../README.md) section “Key data files: H5 and H5AD”)

## BEYOND

The script applies the [BEYOND](https://doi.org/10.1038/s41586-024-07871-6) workflow (Green et al., Nature 2024) to the program atlas: pseudotime, trajectories, community structure, and dynamics.
