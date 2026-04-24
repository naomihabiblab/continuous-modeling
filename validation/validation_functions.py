"""
Validation functions for topic analysis.
Python equivalent of validation/validation.functions.R
"""

import os
import re
import json
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from typing import Dict, List, Optional, Union, Tuple
from scipy import stats
from pathlib import Path


# --- Path helpers (mirror R/utils.R get_result_dir, get_shared_dir) ---

CLUSTER_PLOTS_DIR = "/ems/elsc-labs/habib-n/roi.meir/500/Plots/"
LOCAL_PLOTS_DIR = "/Volumes/habib-lab/roi.meir/500/Plots"
CLUSTER_SHARED_DIR = "/ems/elsc-labs/habib-n/Shared/"
LOCAL_SHARED_DIR = "/Volumes/habib-lab/Shared/"


def get_result_dir() -> str:
    """
    Return the path to the results/plots directory.
    Uses cluster path if it exists, otherwise local mount path (same logic as R get_result_dir).
    """
    if os.path.isdir(CLUSTER_PLOTS_DIR):
        return CLUSTER_PLOTS_DIR
    return LOCAL_PLOTS_DIR


def get_shared_dir() -> str:
    """
    Return the path to the shared data directory.
    Uses cluster path if it exists, otherwise local mount path (same logic as R get_shared_dir).
    """
    if os.path.isdir(CLUSTER_SHARED_DIR):
        return CLUSTER_SHARED_DIR
    return LOCAL_SHARED_DIR


# Constants (use path helpers so cluster vs local is detected at import time)
TOPICS_DATA_H5_PATH = os.path.join(get_result_dir(), "topics_data.h5")

CELL_TYPE_TO_GROUP = {
    "inhibitory": "neuronal/inhibitory",
    "cux2p": "neuronal/excitatory/cux2plus",
    "cux2m": "neuronal/excitatory/cux2minus",
    "excitatory": "neuronal/excitatory",
    "microglia": "glia/microglia",
    "astrocytes": "glia/astrocytes",
    "oligo": "glia/oligodendrocytes",
    "opcs": "glia/opcs",
}

NEURONAL_TOPICS_TO_ANNOTATION = {
    "decreasing_early": {"Inh.10.k1": "decreasing_early", "Exc.Cux2-.7.k4": "decreasing_early", "Exc.Cux2+.8.k9": "decreasing_early"},
    "decreasing_late": {"Inh.10.k7": "decreasing_late", "Exc.Cux2-.7.k8": "decreasing_late", "Exc.Cux2+.8.k10": "decreasing_late"},
    "increasing_early": {"Inh.10.k10": "increasing_early", "Exc.Cux2-.7.k3": "increasing_early", "Exc.Cux2+.8.k11": "increasing_early"},
    "increasing_late": {"Inh.10.k8": "increasing_late", "Exc.Cux2-.7.k9": "increasing_late", "Exc.Cux2+.8.k4": "increasing_late"},
}


def convert_topic_name(name: str) -> str:
    """
    Convert topic name to standardized format.
    
    Args:
        name: Topic name to convert
        
    Returns:
        Converted topic name
    """
    # Remove dots, convert to lower, replace ".x." or ".X." with "_"
    name = re.sub(r"\.", "_", name)
    name = name.lower()
    # If pattern like cux2m_x_k1, convert to cux2m_k1
    name = re.sub(r"_\d+_", "_", name)
    name = re.sub(r"exc_", "", name)
    name = re.sub(r"inh_", "", name)
    return name


def get_genes_by_topic(
    cell_type: str,
    h5_path: str = TOPICS_DATA_H5_PATH,
    top_n: int = 50
) -> Dict[str, Dict[str, List[str]]]:
    """
    Get gene lists by topic from H5 file.
    
    Args:
        cell_type: Cell type to get genes for
        h5_path: Path to H5 file
        top_n: Number of top genes to select by z_score_log
        
    Returns:
        Dictionary with keys: genes_by_topic, genes_by_topic_top_n, unique_genes_by_topic
    """
    with h5py.File(h5_path, "r") as h5file:
        if cell_type == "excitatory":
            des_cux2p = pd.DataFrame(h5file[CELL_TYPE_TO_GROUP["cux2p"] + '/de'][()])
            des_cux2m = pd.DataFrame(h5file[CELL_TYPE_TO_GROUP["cux2m"] + '/de'][()])
            des = pd.concat([des_cux2p, des_cux2m], ignore_index=True)
            # Add cell type prefix to topic
            if "cell.type" in des.columns:
                des["topic"] = des["cell.type"].astype(str) + "_" + des["topic"].astype(str)
            elif "cell_type" in des.columns:
                des["topic"] = des["cell_type"].astype(str) + "_" + des["topic"].astype(str)
        else:
            des = pd.DataFrame(h5file[CELL_TYPE_TO_GROUP[cell_type] + '/de'][()])
    
    # Convert bytes to string if necessary
    for col in ['gene', 'topic']:
        if col in des.columns and des[col].dtype == object:
            des[col] = des[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else str(x))
    
    # Create genes_by_topic
    genes_by_topic = des.groupby('topic')['gene'].apply(list).to_dict()
    
    # Add union of genes for neuronal topics
    if cell_type in ["inhibitory", "excitatory"]:
        for union_name, topic_names_raw in NEURONAL_TOPICS_TO_ANNOTATION.items():
            topic_names = [convert_topic_name(name) for name in topic_names_raw.keys()]
            # Only keep topic names that exist in genes_by_topic
            topic_names = [name for name in topic_names if name in genes_by_topic]
            if len(topic_names) > 0:
                union_genes = list(set(
                    gene for topic in topic_names 
                    for gene in genes_by_topic[topic]
                ))
                genes_by_topic[union_name] = union_genes
        
        # Filter to only union topics for neuronal cell types
        genes_by_topic = {k: v for k, v in genes_by_topic.items() 
                         if k in NEURONAL_TOPICS_TO_ANNOTATION.keys()}
    
    # Create genes_by_topic_top_n (top N genes by z_score_log)
    if 'z_score_log' in des.columns:
        genes_by_topic_top_n = (
            des.groupby('topic')
            .apply(lambda x: x.nlargest(top_n, 'z_score_log')['gene'].tolist())
            .to_dict()
        )
    else:
        genes_by_topic_top_n = {k: v[:top_n] for k, v in genes_by_topic.items()}
    
    # Create unique_genes_by_topic (genes appearing in only one topic)
    # Flatten genes_by_topic to get gene-topic pairs
    gene_topic_pairs = [
        (gene, topic) 
        for topic, genes in genes_by_topic.items() 
        for gene in genes
    ]
    gene_topic_df = pd.DataFrame(gene_topic_pairs, columns=['gene', 'topic'])
    
    # Count how many topics each gene appears in
    gene_counts = gene_topic_df.groupby('gene')['topic'].nunique()
    unique_genes = set(gene_counts[gene_counts == 1].index)
    
    # Filter to unique genes per topic
    unique_genes_by_topic = {
        topic: [gene for gene in genes if gene in unique_genes]
        for topic, genes in genes_by_topic.items()
    }
    
    return {
        'genes_by_topic': genes_by_topic,
        'genes_by_topic_top_n': genes_by_topic_top_n,
        'unique_genes_by_topic': unique_genes_by_topic
    }


def get_pathway_signatures(
    cell_type: str,
    results_dir: Optional[str] = None
) -> Dict[str, List[str]]:
    """
    Get pathway gene signatures.
    
    Args:
        cell_type: Cell type to get pathways for
        results_dir: Base results directory (default: get_result_dir())
        
    Returns:
        Dictionary mapping pathway names to gene lists
    """
    if results_dir is None:
        results_dir = get_result_dir()
    # Try to load from JSON file
    pathways_to_genes_path = os.path.join(results_dir, "topic_space", f"Pathways_genes_{cell_type}.json")
    
    pathway_to_genes = {}
    if os.path.exists(pathways_to_genes_path):
        with open(pathways_to_genes_path, 'r') as f:
            pathway_to_genes = json.load(f)
    
    # For neuronal cell types, also load from natasha directory
    if cell_type in ["inhibitory", "excitatory"]:
        natasha_path = os.path.join(results_dir, "natasha", "cuimc2", cell_type, "PathwayToGenes.json")
        if os.path.exists(natasha_path):
            with open(natasha_path, 'r') as f:
                neuronal_pathways = json.load(f)
                pathway_to_genes.update(neuronal_pathways)
    
    return pathway_to_genes



def compute_gene_set_coherence(
    expr_matrix: np.ndarray,
    genes: Union[List[str], List[int]],
    gene_names: Optional[List[str]] = None,
    n_permutations: int = 50,
    method: str = 'pearson',
    seed: int = 1,
    show_progress: bool = True,
    n_cores: Optional[int] = None,
    temp_folder: Optional[str] = None
) -> Dict:
    """
    Compute mean pairwise correlation for a gene set and compare to random background.
    
    This is the Python equivalent of compute_gene_set_coherence() in validation.functions.R
    
    Args:
        expr_matrix: Expression matrix (genes x cells), e.g., from adata.X or adata.raw.X
        genes: Either a list of gene names (strings) or gene indices (integers)
        gene_names: List of all gene names in expr_matrix (required if genes are strings).
                   Typically adata.var_names or adata.var['gene_name']
        n_permutations: Number of random background permutations (default: 50)
        method: Correlation method: 'spearman' or 'pearson' (default: 'pearson')
        seed: Random seed for reproducibility (default: 1)
        show_progress: Show progress bar (default: True, ignored when n_cores > 1)
        n_cores: Number of cores for parallel processing (default: all available cores)
                Set to 1 for sequential processing with progress bar
        temp_folder: Directory for joblib memory mapping temp files (default: None uses system temp).
                    Use this if /tmp or /dev/shm is full. Example: '/scratch/tmp'
        
    Returns:
        Dictionary with observed_mean, bg_means, p_value, z_score, effect sizes
        
    Note:
        Caching pairwise correlations is NOT recommended because:
        1. Random permutations use different gene sets each time (minimal overlap)
        2. Pre-computing all pairs (20k x 20k) requires ~3GB RAM + huge compute time
        3. The bottleneck is correlating across cells, not looking up gene pairs
    """
    import multiprocessing
    
    if seed is not None:
        np.random.seed(seed)
    
    # Determine number of cores
    if n_cores is None:
        n_cores = multiprocessing.cpu_count()
    
    # Handle gene names vs indices
    if gene_names is not None and len(genes) > 0 and isinstance(genes[0], str):
        # Convert gene names to indices
        gene_name_to_idx = {name: i for i, name in enumerate(gene_names)}
        valid_genes = [g for g in genes if g in gene_name_to_idx]
        gene_indices = [gene_name_to_idx[g] for g in valid_genes]
        all_gene_indices = list(range(len(gene_names)))
    else:
        # Assume genes are already indices
        gene_indices = list(genes)
        valid_genes = gene_indices
        all_gene_indices = list(range(expr_matrix.shape[0]))
    
    n_genes = len(gene_indices)
    
    if n_genes < 2:
        print("Warning: Less than 2 valid genes found. Cannot compute correlation.")
        return {
            'observed_mean': np.nan,
            'bg_means': np.nan,
            'bg_mean': np.nan,
            'bg_sd': np.nan,
            'p_value': np.nan,
            'z_score': np.nan,
            'cohens_d': np.nan,
            'percentile_rank': np.nan,
            'fold_enrichment': np.nan,
            'cliffs_delta': np.nan,
            'n_genes': n_genes,
            'n_permutations': n_permutations
        }
    
    def compute_mean_cor(indices, expr_mat, corr_method):
        """
        Compute mean pairwise correlation between genes (rows).
        expr_sub has shape (n_genes, n_cells).
        We want gene-gene correlation, NOT cell-cell correlation.
        
        Sparse-aware: uses sparse operations to avoid converting to dense.
        Only the small (n_genes x n_genes) correlation matrix is dense.
        """
        from scipy import sparse
        
        expr_sub = expr_mat[indices, :]
        n_genes_sub = expr_sub.shape[0]
        n_cells = expr_sub.shape[1]
        
        if n_genes_sub < 2:
            return np.nan
        
        is_sparse = sparse.issparse(expr_sub)
        
        if corr_method == 'spearman':
            # Spearman requires ranking, which needs dense data
            # Convert only the subset (small: n_genes x n_cells)
            if is_sparse:
                expr_dense = expr_sub.toarray()
            else:
                expr_dense = np.asarray(expr_sub)
            result = stats.spearmanr(expr_dense, axis=1)
            cor_mat = result.correlation if hasattr(result, 'correlation') else result[0]
        else:
            # Pearson correlation - sparse-aware computation
            # Formula: cor(i,j) = (E[Xi*Xj] - μi*μj) / (σi * σj)
            
            if is_sparse:
                # Sparse-efficient computation
                # Row means: mean of each gene across cells
                row_means = np.asarray(expr_sub.mean(axis=1)).ravel()  # (n_genes,)
                
                # Row variances using E[X²] - E[X]²
                row_sq_means = np.asarray(expr_sub.power(2).mean(axis=1)).ravel()
                row_vars = row_sq_means - row_means ** 2
                row_stds = np.sqrt(np.maximum(row_vars, 0))  # Avoid negative from numerical issues
                
                # Dot products: E[Xi * Xj] * n_cells = Xi @ Xj.T
                # (n_genes, n_cells) @ (n_cells, n_genes) = (n_genes, n_genes) - small!
                dot_products = (expr_sub @ expr_sub.T).toarray()  # Only this is dense
                
                # Covariance: E[Xi*Xj] - μi*μj
                cov = dot_products / n_cells - np.outer(row_means, row_means)
                
                # Correlation: cov / (σi * σj)
                std_outer = np.outer(row_stds, row_stds)
                std_outer = np.where(std_outer == 0, np.nan, std_outer)
                cor_mat = cov / std_outer
                
                # Clip to [-1, 1] due to numerical precision
                cor_mat = np.clip(cor_mat, -1.0, 1.0)
            else:
                # Dense computation
                expr_sub = np.asarray(expr_sub)
                mean = expr_sub.mean(axis=1, keepdims=True)
                Xc = expr_sub - mean
                cov = (Xc @ Xc.T) / (n_cells - 1)
                std = np.sqrt(np.diag(cov))
                std = np.where(std == 0, np.nan, std)
                cor_mat = cov / np.outer(std, std)
        
        # Handle case where cor_mat is scalar (2 genes)
        if np.isscalar(cor_mat):
            return float(cor_mat)
        
        # Ensure 2D
        cor_mat = np.atleast_2d(cor_mat)
        
        # Extract upper triangle (excluding diagonal) and compute mean
        upper_tri = cor_mat[np.triu_indices_from(cor_mat, k=1)]
        return np.nanmean(upper_tri)
    
    # Compute observed mean correlation for the gene set
    print(f"Computing observed mean correlation for gene set of size {n_genes} using {n_cores} cores...")
    import time
    start_time = time.time()
    observed_mean = compute_mean_cor(gene_indices, expr_matrix, method)
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.2f} seconds")
    
    # Compute background distribution using random gene sets of the same size
    print(f"Computing background distribution ({n_permutations} permutations, using {n_cores} cores)...")
    
    # Pre-generate all random gene sets for reproducibility
    if seed is not None:
        np.random.seed(seed + 1)  # Different seed for background
    random_gene_sets = [
        np.random.choice(all_gene_indices, size=n_genes, replace=False).tolist()
        for _ in range(n_permutations)
    ]
    
    if n_cores > 1:
        # Memory-efficient parallel computation using joblib
        # joblib uses memory mapping to share large arrays between processes
        # This avoids the pickle overhead of multiprocessing.Pool
        try:
            from joblib import Parallel, delayed
            import os
            
            # Set temp folder for joblib memory mapping if specified
            # This is useful when /tmp or /dev/shm is full
            if temp_folder is not None:
                os.makedirs(temp_folder, exist_ok=True)
                os.environ['JOBLIB_TEMP_FOLDER'] = temp_folder
                print(f"  Using temp folder: {temp_folder}")
            
            # Try to use tqdm for a nice progress bar
            try:
                from tqdm import tqdm
                tasks = tqdm(random_gene_sets, desc="  Permutations", disable=not show_progress)
            except ImportError:
                tasks = random_gene_sets
                if show_progress:
                    print("  (Install tqdm for progress bar: pip install tqdm)")
            
            # max_nbytes='100M' triggers memory mapping for arrays > 100MB
            # This means the expr_matrix is memory-mapped, not copied to each worker
            # temp_folder redirects mmap files away from /tmp if specified
            bg_means = np.array(
                Parallel(n_jobs=n_cores, backend='loky', max_nbytes='100M', 
                         temp_folder=temp_folder)(
                    delayed(compute_mean_cor)(indices, expr_matrix, method) 
                    for indices in tasks
                )
            )
        except ImportError:
            # Fallback to sequential if joblib not available
            print("  Warning: joblib not installed, falling back to sequential processing")
            print("  Install with: pip install joblib")
            n_cores = 1
    
    if n_cores == 1:
        # Sequential computation with progress bar
        bg_means = np.zeros(n_permutations)
        start_time = time.time()
        
        for i in range(n_permutations):
            bg_means[i] = compute_mean_cor(random_gene_sets[i], expr_matrix, method)
            
            if show_progress:
                elapsed = time.time() - start_time
                avg_per_iter = elapsed / (i + 1)
                remaining = avg_per_iter * (n_permutations - i - 1)
                pct = 100 * (i + 1) / n_permutations
                
                # Format time
                def format_time(secs):
                    secs = int(secs)
                    if secs < 3600:
                        return f"{secs // 60:02d}:{secs % 60:02d}"
                    else:
                        return f"{secs // 3600:02d}:{(secs % 3600) // 60:02d}:{secs % 60:02d}"
                
                bar_width = 30
                filled = int(bar_width * (i + 1) / n_permutations)
                bar = '=' * filled + ' ' * (bar_width - filled)
                
                if (i + 1) == 1 or (i + 1) == n_permutations or (i + 1) % max(1, n_permutations // 20) == 0:
                    print(f"\r[{bar}] {pct:3.0f}% | {i+1}/{n_permutations} | "
                          f"Elapsed: {format_time(elapsed)} | ETA: {format_time(remaining)}   ", 
                          end='', flush=True)
        
        if show_progress:
            total_time = time.time() - start_time
            print(f"\n  Completed in {total_time:.1f} seconds" if total_time < 60 
                  else f"\n  Completed in {total_time/60:.1f} minutes")
    
    # Compute empirical p-value (one-tailed: is observed higher than random?)
    # Add 1 to numerator and denominator for conservative estimate (Phipson & Smyth, 2010)
    p_value = (np.sum(bg_means >= observed_mean) + 1) / (n_permutations + 1)
    
    # Compute effect size metrics
    bg_mean = np.nanmean(bg_means)
    bg_sd = np.nanstd(bg_means, ddof=1)  # Use ddof=1 to match R's sd()
    
    # Z-score (standardized difference)
    z_score = (observed_mean - bg_mean) / bg_sd if bg_sd > 0 else np.nan
    
    # Cohen's d (same as z-score in this context, but conventionally reported)
    cohens_d = z_score
    
    # Percentile rank: proportion of background values below observed
    # Higher = observed is more extreme (better coherence)
    percentile_rank = np.mean(bg_means < observed_mean) * 100
    
    # Fold enrichment: ratio of observed to background mean
    # >1 means gene set is more coherent than random
    fold_enrichment = observed_mean / bg_mean if bg_mean != 0 else np.nan
    
    # Cliff's delta: non-parametric effect size
    # Ranges from -1 to 1; >0 means observed tends to be higher than background
    # Formula: (# times obs > bg - # times obs < bg) / total comparisons
    cliffs_delta = (np.sum(observed_mean > bg_means) - np.sum(observed_mean < bg_means)) / n_permutations
    
    return {
        'observed_mean': observed_mean,
        'bg_means': bg_means,
        'bg_mean': bg_mean,
        'bg_sd': bg_sd,
        'p_value': p_value,
        'z_score': z_score,
        'cohens_d': cohens_d,
        'percentile_rank': percentile_rank,
        'fold_enrichment': fold_enrichment,
        'cliffs_delta': cliffs_delta,
        'n_genes': n_genes,
        'n_permutations': n_permutations
    }


def compute_coherence_multiple(
    expr_matrix: np.ndarray,
    gene_sets: Dict[str, List[str]],
    gene_names: Optional[List[str]] = None,
    n_permutations: int = 50,
    method: str = 'pearson',
    seed: int = 1,
    show_progress: bool = True,
    n_cores: Optional[int] = None,
    temp_folder: Optional[str] = None
) -> pd.DataFrame:
    """
    Compute coherence for multiple gene sets with FDR correction.
    
    This is the Python equivalent of compute_coherence_multiple() in validation.functions.R
    
    Args:
        expr_matrix: Expression matrix (genes x cells)
        gene_sets: Named dictionary of gene lists (gene names per set)
        gene_names: List of all gene names in expr_matrix (required if genes are strings)
        n_permutations: Number of random background permutations (default: 50)
        method: Correlation method: 'spearman' or 'pearson' (default: 'pearson')
        seed: Random seed (default: 1)
        show_progress: Show progress bar for each gene set (ignored when n_cores > 1)
        n_cores: Number of cores for parallel processing (default: all available cores)
        temp_folder: Directory for joblib memory mapping temp files (default: None uses system temp).
                    Use this if /tmp or /dev/shm is full. Example: '/scratch/tmp'
        
    Returns:
        DataFrame with coherence statistics for each gene set, including FDR-adjusted p-values
    """
    import multiprocessing
    
    if n_cores is None:
        n_cores = multiprocessing.cpu_count()
    
    print(f"Using {n_cores} cores for parallel processing")
    
    if seed is not None:
        np.random.seed(seed)
    
    results = []
    for name, genes in gene_sets.items():
        print(f"\nComputing coherence for {name} of size {len(genes)}...")
        res = compute_gene_set_coherence(
            expr_matrix=expr_matrix,
            genes=genes,
            gene_names=gene_names,
            n_permutations=n_permutations,
            method=method,
            seed=None,  # Don't reset seed for each gene set
            show_progress=show_progress,
            n_cores=n_cores,
            temp_folder=temp_folder
        )
        results.append({
            'gene_set': name,
            'n_genes': res['n_genes'],
            'observed_mean_cor': res['observed_mean'],
            'bg_means_cor': res['bg_means'],
            'bg_mean_cor': res['bg_mean'],
            'bg_sd_cor': res['bg_sd'],
            'z_score': res['z_score'],
            'cohens_d': res['cohens_d'],
            'percentile_rank': res['percentile_rank'],
            'fold_enrichment': res['fold_enrichment'],
            'cliffs_delta': res['cliffs_delta'],
            'p_value': res['p_value'],
        })
    
    result_df = pd.DataFrame(results)
    
    # Adjust p-values for multiple testing (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    valid_pvals = result_df['p_value'].notna()
    if valid_pvals.sum() > 0:
        _, result_df.loc[valid_pvals, 'p_adj'], _, _ = multipletests(
            result_df.loc[valid_pvals, 'p_value'], method='fdr_bh'
        )
    else:
        result_df['p_adj'] = np.nan
    
    return result_df
