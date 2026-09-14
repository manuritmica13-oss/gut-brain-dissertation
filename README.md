# Gut-Brain Axis Metabolites as Early Biomarkers of Parkinson's Disease

**MSc Dissertation — Engineering for Biomedicine**
**Imperial College London | Professor Parastoo Hashemi's Research Group**
**Author: Manuela Coronado**

---

## Overview

This repository contains the full bioinformatics and machine learning pipeline developed for my MSc dissertation, investigating whether gut microbiome functional profiles from Major Depressive Disorder (MDD) patients can predict preclinical Parkinson's Disease (PD) via tryptophan-kynurenine pathway dysregulation.

**Research Question:** Can gut microbiome functional profiles from MDD cohorts predict preclinical Parkinson's Disease through kynurenine pathway dysregulation?

**Biological rationale:** MDD patients show elevated IDO1 activity and disrupted tryptophan-kynurenine metabolism — overlapping with pathways implicated in PD neurodegeneration. Whether this gut-brain axis dysregulation precedes PD onset remains unclear.

---

## Pipeline Overview

Raw sequencing data → QIIME2 (amplicon processing) → HUMAnN3 (functional pathway analysis) → MMUPHin (batch correction) → SMOTE (class imbalance handling) → RF / XGBoost classifiers → SHAP (feature attribution) → LOSO cross-validation

---

## Repository Structure

data/ — Metadata and processed feature tables
envs/ — Conda environment files
scripts/ — Analysis scripts (Python and R)
results/ — Output figures, model results, and reports
CLAUDE.md — AI usage declaration

---

## Methods

### Datasets
- 7 independent metagenomic cohorts spanning MDD and PD patient populations
- ~1,884 samples across studies from multiple countries
- Cohorts were sourced from publicly available databases and processed independently before integration

### Tools & Libraries

QIIME2 — Amplicon sequence processing and taxonomic classification
HUMAnN3 — Functional pathway reconstruction from metagenomic data
MMUPHin — Batch effect correction across heterogeneous cohorts
scikit-learn — Random Forest classifier
XGBoost — Gradient boosting classifier
SMOTE — Synthetic oversampling for class imbalance correction
SHAP — Model interpretability and feature attribution
Python / R — Pipeline development and statistical analysis
Imperial HPC (SLURM) — High-performance computing cluster

### Evaluation
- Primary metric: AUPRC (Area Under the Precision-Recall Curve) — chosen over AUROC due to class imbalance
- Validation strategy: Leave-One-Study-Out (LOSO) cross-validation to assess generalisation across cohorts
- Permutation testing: 1,000 permutation iterations to confirm results are not due to chance

---

## Key Findings

- Kynurenine pathway dysregulation — specifically IDO1-driven tryptophan catabolism — was identified as a key differentiating signal between MDD and preclinical PD microbiome profiles
- SMOTE-balanced RF and XGBoost classifiers demonstrated above-chance AUPRC under LOSO validation
- SHAP feature attribution identified specific kynurenine pathway metabolite features driving classification

---

## AI Usage Declaration

This project was developed with the assistance of AI tools, in accordance with Imperial College London's guidance on responsible AI use in academic work.

Tools used:

Claude (Anthropic) — VS Code Extension: Used throughout the development of this pipeline for code debugging, methodology discussion, structuring analytical approaches, and drafting documentation. Claude was used as an interactive coding assistant — all code was reviewed, validated, and understood by the author before implementation. Claude did not generate any scientific conclusions independently; all interpretation and research decisions were made by the author.

GitHub Copilot: Used for code autocompletion during scripting. All suggestions were reviewed and verified by the author.

Scope of AI assistance: AI tools assisted with the process of writing and debugging code, not with scientific reasoning or interpretation. All analytical decisions — choice of metrics, validation strategy, model selection, and biological interpretation — were made independently by the author.

The use of these tools is declared in full transparency. The intellectual content, scientific judgement, and research contributions of this dissertation are entirely the author's own.

---


---

This repository is made public for transparency and reproducibility purposes. Data used in this study are publicly available through their respective original publications.
