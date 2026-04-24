# =============================================================================
# PATH CONFIGURATION FOR SPATIAL ANALYSIS
# =============================================================================
"""
Centralized path configuration for spatial transcriptomics analysis.
Automatically detects cluster vs local environment.


Required data:
- Raw GEF files from GSE269906 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269906)
- 500 project topics data (topics_data.h5)

Outputs:
- H5AD files: {SPATIAL_DATA_DIR}/h5ad_files/
- Processed objects: {SPATIAL_DATA_DIR}/for_paper/
- Plots: {PLOTS_DIR}/spatial/

"""

import os
import glob

# =============================================================================
# PATH DETECTION FUNCTIONS
# =============================================================================

def get_plots_dir():
    """
    Get plots output directory.
    
    Returns:
        str: Path to plots directory
        
    Raises:
        FileNotFoundError: If neither cluster nor local path exists
    """
    cluster_path = "/ems/elsc-labs/habib-n/roi.meir/500/Plots"
    local_path = "/Volumes/habib-lab/roi.meir/500/Plots"
    
    if os.path.exists(cluster_path):
        return cluster_path
    elif os.path.exists(local_path):
        return local_path
    else:
        raise FileNotFoundError(
            "Cannot find plots directory. "
            f"Tried:\n  - {cluster_path}\n  - {local_path}"
        )

def get_spatial_base_dir():
    """
    Detect if running on cluster or local mount for Gali's spatial data.
    
    Returns:
        str: Path to spatial data directory
        
    Raises:
        FileNotFoundError: If neither cluster nor local path exists
    """
    plots_dir = get_plots_dir()
    return os.path.join(plots_dir, "spatial")


def get_shared_dir():
    """
    Get shared data directory (500 project data).
    
    Returns:
        str: Path to shared directory
        
    Raises:
        FileNotFoundError: If neither cluster nor local path exists
    """
    cluster_path = "/ems/elsc-labs/habib-n/Shared"
    local_path = "/Volumes/habib-lab/Shared"
    
    if os.path.exists(cluster_path):
        return cluster_path
    elif os.path.exists(local_path):
        return local_path
    else:
        raise FileNotFoundError(
            "Cannot find shared directory. "
            f"Tried:\n  - {cluster_path}\n  - {local_path}"
        )


def get_gef_path(sample_id: str, raw_data_dir: str) -> str:
    """
    Return absolute path to the .cellbin.gef for a given sample_id.
    
    Uses glob pattern GSM*_{sample_id}.cellbin.gef to find the file.
    
    Args:
        sample_id: Sample identifier (e.g., "A02092E1")
        raw_data_dir: Directory containing raw GEF files
    
    Returns:
        str: Path to GEF file
        
    Raises:
        FileNotFoundError: If no matching GEF file found
    """
    pattern = os.path.join(raw_data_dir, f"GSM*_{sample_id}.cellbin.gef")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No .cellbin.gef file found for sample '{sample_id}' "
            f"under {raw_data_dir} matching pattern {os.path.basename(pattern)}"
        )
    return matches[0]

# =============================================================================
# DIRECTORY PATHS
# =============================================================================

# Base directories (auto-detected)
SPATIAL_BASE_DIR = get_spatial_base_dir()
SHARED_DIR = get_shared_dir()
PLOTS_DIR = get_plots_dir()

# Raw data from GEO
RAW_DATA_DIR = os.path.join(SHARED_DIR, "GSE269906_RAW")

# Processed data output directories
H5AD_OUTPUT_DIR = os.path.join(SPATIAL_BASE_DIR, "h5ad_files/sc_26_10_new")
PAPER_OUTPUT_DIR = os.path.join(SPATIAL_BASE_DIR, "for_paper")
PAPER_H5AD_DIR = os.path.join(PAPER_OUTPUT_DIR, "h5ad_files")
CELL_TYPES_DIR = os.path.join(PAPER_OUTPUT_DIR, "cell_types")
FINAL_OBJ_DIR = os.path.join(PAPER_OUTPUT_DIR, "final_obj")

# Create output directories if they don't exist
os.makedirs(H5AD_OUTPUT_DIR, exist_ok=True)
os.makedirs(PAPER_H5AD_DIR, exist_ok=True)
os.makedirs(CELL_TYPES_DIR, exist_ok=True)
os.makedirs(FINAL_OBJ_DIR, exist_ok=True)

# 500 project data (topics, DE results, etc.)
TOPICS_DATA_PATH = os.path.join(PLOTS_DIR, "500_topics/topics_data.h5")

# Output for plots/figures
SPATIAL_PLOTS_DIR = SPATIAL_BASE_DIR
os.makedirs(SPATIAL_PLOTS_DIR, exist_ok=True)

# Enrichment data paths (for script_for_paper.ipynb)
ENRICHMENT_BASE_DIR = os.path.join(PAPER_OUTPUT_DIR, "enrichment_data")
ENRICHMENT_PROGRAMS_DIR = os.path.join(ENRICHMENT_BASE_DIR, "programs")
ENRICHMENT_PATHWAYS_DIR = os.path.join(ENRICHMENT_BASE_DIR, "pathways")
ENRICHMENT_PROGRAMS_PATH = os.path.join(ENRICHMENT_PROGRAMS_DIR, "all_pvals_programs.csv")
ENRICHMENT_PATHWAYS_PATH = os.path.join(ENRICHMENT_PATHWAYS_DIR, "all_pvals_pathways.csv")

# Create enrichment directories if they don't exist
os.makedirs(ENRICHMENT_PROGRAMS_DIR, exist_ok=True)
os.makedirs(ENRICHMENT_PATHWAYS_DIR, exist_ok=True)

# =============================================================================
# SPECIFIC DATA FILE PATHS
# =============================================================================

# Final processed spatial object
SPATIAL_OBJECT_PATH = os.path.join(FINAL_OBJ_DIR, "spatial_obj_for_paper.h5ad")

# Annotated merged object
ALL_ANNOTATED_PATH = os.path.join(H5AD_OUTPUT_DIR, "all_obj_annotated.h5ad")

# =============================================================================
# CONFIGURATION SUMMARY
# =============================================================================

# Determine environment
IS_CLUSTER = "ems" in SPATIAL_BASE_DIR
ENVIRONMENT = "CLUSTER" if IS_CLUSTER else "LOCAL"

# Print configuration
print("=" * 80)
print("SPATIAL ANALYSIS CONFIGURATION")
print("=" * 80)
print(f"Environment: {ENVIRONMENT}")
print()
print("Data Directories:")
print(f"  Spatial base: {SPATIAL_BASE_DIR}")
print(f"  Raw GEF files: {RAW_DATA_DIR}")
print(f"  Shared (500 project): {SHARED_DIR}")
print()
print("Output Directories:")
print(f"  H5AD files: {H5AD_OUTPUT_DIR}")
print(f"  Paper objects: {FINAL_OBJ_DIR}")
print(f"  Plots: {SPATIAL_PLOTS_DIR}")
print()
print("Key Data Files:")
print(f"  Topics data: {TOPICS_DATA_PATH}")
print(f"  Spatial object: {SPATIAL_OBJECT_PATH}")
print()
print("=" * 80)

# =============================================================================
# VALIDATION
# =============================================================================

# Verify critical paths exist
if not os.path.exists(TOPICS_DATA_PATH):
    print(f"⚠️  WARNING: Topics data not found: {TOPICS_DATA_PATH}")
    print("   This file is required for computing topic scores.")
    
if not os.path.exists(RAW_DATA_DIR):
    print(f"⚠️  WARNING: Raw data directory not found: {RAW_DATA_DIR}")
    print("   This directory is required for data_processing_paper.ipynb")

# =============================================================================
# SAMPLE DEFINITIONS
# =============================================================================

# Sample lists for analysis
CASE_SAMPLES = ["A02092E1", "B02009F6", "C02248B5", "B02008C6", "B02008D2", "B01809C2"]
CONTROL_SAMPLES = ["D02175A4", "D02175A6", "B01809A3", "B01809A4", "B01806B5", "B01806B6"]

# Sample metadata (for grouping by disease severity)
SAMPLE_METADATA = {
    # Control
    "B01806B5": {"diagnosis": "control", "level": "control"},
    "B01806B6": {"diagnosis": "control", "level": "control"},
    "B01809A3": {"diagnosis": "control", "level": "control"},
    "B01809A4": {"diagnosis": "control", "level": "control"},
    "D02175A4": {"diagnosis": "control", "level": "control"},
    "D02175A6": {"diagnosis": "control", "level": "control"},
    # Moderate
    "B02008D2": {"diagnosis": "case", "level": "moderate"},
    "B02009F6": {"diagnosis": "case", "level": "moderate"},
    # Advanced
    "B01809C2": {"diagnosis": "case", "level": "advanced"},
    "C02248B5": {"diagnosis": "case", "level": "advanced"},
    # Severe
    "A02092E1": {"diagnosis": "case", "level": "severe"},
    "B02008C6": {"diagnosis": "case", "level": "severe"},
}

moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]

# =============================================================================
# END CONFIGURATION
# =============================================================================
