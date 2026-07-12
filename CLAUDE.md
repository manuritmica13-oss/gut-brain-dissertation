# Gut-Brain Axis Dissertation

MSc Bioengineering, Imperial College London. Supervisor: Dr. Hashemi.

## Project
Computational modelling of gut microbiome functional profiles across MDD and PD cohorts. Testing whether depression-based gut signatures can predict preclinical Parkinson's (AUC > 0.75).

## Structure
- data/raw/ — downloaded datasets (NCBI SRA / ENA, never modify)
- data/processed/ — QC outputs, batch-corrected feature tables
- scripts/qc/ — Phase 1: QIIME2, HUMAnN3, MMUPHin
- scripts/clustering/ — Phase 2: PERMANOVA, DESeq2, ANCOMBC
- scripts/network/ — Phase 3: WGCNA
- scripts/ml/ — Phase 4: Random Forest, XGBoost, SHAP, DCA
- scripts/biomarkers/ — Phase 5: annotation, Metabolomics Workbench
- results/ — all outputs and figures
- envs/ — conda environment files

## Datasets
PD: PRJNA834801, PRJEB55464, PRJNA808166, PRJNA588035
MDD: PRJNA510556, PRJNA762199, PRJNA1083304

## Languages
Python and R. Conda for environment management. HPC: Imperial College cluster.
