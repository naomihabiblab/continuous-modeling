# Early neuronal reprogramming and cell cycle reentry shape Alzheimer’s disease progression paper's code

The code is divided into three parts
1. **Cell-type analysis**: snRNA-seq NMF-based co-expression gene programs analysis for each of the cell types. Program signature (differential expressed genes) and pathway enrichment analyses. 
   
    R files under `cell_type_analysis`. 
   1. First run {cell.type}TopicRunner.R to fit the programs for the cel type
   2. Run {cell.type}_topics.R to compute Program signature (differential expressed genes), pathway enrichment analyses and save the results

2. **Program space analysis**: Trait associations and BEYOND - `topic_space.R` 
   1. **BEYOND**: Applying BEYOND methodology ([*Green et al*, Nature 2024](https://doi.org/10.1038/s41586-024-07871-6)) steps over our program atlas.
3. Validation
   1. Computing signature score in SEA-AD data ([*Gabitto et al. Nat Neurosci 2024*](https://doi.org/10.1038/s41593-024-01774-5)) - jupyter notebook at `validation/sead`
4. **Manuscript code**: The complete code generating all figures, extended data figures and supplementary tables shown in the manuscript.
   - `figures.R`
   - `create.supplementary.tables.R`

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

## DATA and requirements 
<details>
<summary>Data accessibility</summary>
All snRNA-seq data used in our study is accessible via [Synapse (syn53366818)](https://www.synapse.org/#!Synapse:syn53366818) and contains the raw reads, library count matrices and processed cell atlas in the format of cell-type Seurat objects. 
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
- SEA_AD signature score: python (3.12), dependencies listed in `validation/sead/requirements.txt` 
</details>