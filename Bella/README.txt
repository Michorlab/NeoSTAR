================================================================================
MANUSCRIPT REVISION README: SUMMARY OF CODE MODIFICATIONS & ADDITIONS
================================================================================
Author: Isabella Pabón
Date: July 2026
Project: Triple-Negative Breast Cancer (TNBC) single-cell RNA-seq Manuscript Revision
================================================================================

This README documents the precise engineering modifications, bug fixes, statistical
safeguards, and publication-grade figures implemented across the analysis pipeline 
for Figures 1, 2, and 3.

--------------------------------------------------------------------------------
FILE 1: Figure 1 Pipeline (Core QC, Annotation UMAPs, & TNBC Subtype Calling)
--------------------------------------------------------------------------------
* Deprecation Fixes:
  - Upgraded older Seurat syntax to prevent fatal errors under modern Seurat versions[cite: 1].
  - Converted the obsolete `slot = "data"` parameter to `layer = "data"` inside the 
    `GetAssayData()` call for dot plots and feature plots[cite: 1].
* Environment & Dependency Management:
  - Added targeted global options (`options(sccomp_models_dir = ...)`) and forced matching 
    environment variables (`Sys.setenv(SCCOMP_MODELS_DIR = ...)`) to dynamically configure 
    local package compilation paths for `sccomp` without breaking on local directory permissions[cite: 1].
* Robust Data Conversions:
  - Fixed a silent failure in the downstream data aggregation pipeline by explicitly wrapping 
    the raw list metrics in `dplyr::as_tibble()` before performing `dplyr::left_join` and 
    `geom_segment` visual plotting[cite: 1]. This ensures clean table structures for FDR text alignment.
* Cosmological & Working Directory Management:
  - Standardized all script evaluations to point cleanly to the local data directory structure[cite: 1]. 

--------------------------------------------------------------------------------
FILE 2: Figure 2 Pipeline (Epithelial Architecture, Heterogeneity, & Metaprograms)
--------------------------------------------------------------------------------
* Deprecation Fixes:
  - Swapped out all instances of the deprecated `slot` parameter in favor of the active 
    `layer = "data"` and `layer = "scale.data"` syntax across all epithelial sub-clustering, 
    module scoring, and `DoHeatmap` evaluations[cite: 2].
* Statistical Robustness & Patient-Level Aggregations (Pseudoreplication Fixes):
  - Engineered a custom Poisson-based simulation pipeline (`n_sim = 10000`) to rigorously 
    evaluate TROP2 (`TACSTD2`) intratumoral heterogeneity (ITH) across treatment cohorts, 
    accounting for structural zero-inflation by standardizing the drop-out effect size (Z-score)[cite: 2].
  - Fixed pseudoreplication issues in down-stream analysis sections by collapsing single-cell 
    metaprogram (MP) data down to true patient-level averages before conducting comparative testing 
    (e.g., Wilcoxon tests across clinical response cohorts now evaluate n = 24 patients rather 
    than n = hundreds of thousands of cells, protecting against falsely deflated p-values)[cite: 2].
* Visual Sophistication & Aesthetics:
  - Refined the multi-rank non-negative matrix factorization (cNMF) plotting script by 
    introducing customized magma gradient palettes coupled with detailed row split divisions[cite: 2].
  - Optimized sample annotation sidebars for the final Jaccard similarity matrices, using 
    clean programmatic rendering blocks to maximize figure reproducibility[cite: 2].

--------------------------------------------------------------------------------
FILE 3: Figure 3 & Supplementary Figure S6 Pipeline (TME & Macrophage Overhauls)
--------------------------------------------------------------------------------
* Primary Boxplot Replacement (The Reviewer Fix):
  - Completely replaced the old cell-level violin expressions in favor of patient-level 
    signature boxplots ("Fig3C_CD8T_boxplots.pdf")[cite: 3]. 
  - Engineered the pipeline to calculate the mean module score per patient sample first, 
    map out directional nudging and labeled jitter dots via geom_text_repel, and apply a 
    robust Wilcoxon test directly to sample means to address reviewer pseudoreplication concerns[cite: 3].
  - Curated and wrapped long signature titles programmatically for a crisp publication layout[cite: 3].

================================================================================
END OF README
================================================================================