# Multi-Omics Integration Project

## Project Overview

This project integrates multi-omics data from TCGA and METABRIC breast cancer cohorts to build a prognostic analysis pipeline. The main goal is to align eigengene features, clinical survival data, and learned embeddings to support survival modeling and cohort-based interpretation.

## Files Included

- `phase1.ipynb`
  - Phase 1 notebook covering data loading, eigengene preprocessing, TCGA and METABRIC cleaning, normalization, and dataset integration.
  - It prepares the core molecular features used for downstream survival and graph-based analysis.

- `phase_2.ipynb`
  - Phase 2 continuation notebook that builds on Phase 1.
  - Performs survival alignment, embedding similarity analysis, graph-based validation, and LLM-style interpretation.
  - Includes survival modeling, validation of patient similarity, and risk profiling.

## Key Concepts

- **Eigengene integration**: Transpose and clean module eigengene matrices from TCGA and METABRIC, then align shared features.
- **Survival alignment**: Match clinical survival records to patient embeddings and eigengenes for prognostic modeling.
- **Graph embeddings**: Use learned embeddings to identify similar patients and validate survival consistency.
- **Interpretation layer**: Build structured clinical interpretation prompts and outputs based on model-derived patient profiles.

## Why these files

The notebook pair (`phase1.ipynb` and `phase_2.ipynb`) represent the complete project workflow. They are meant to be used together, with Phase 2 extending and validating the analysis from Phase 1.

## How to use

1. Open `phase1.ipynb` first to run data preparation and integration.
2. Continue with `phase_2.ipynb` to perform modeling, validation, and interpretation.

## Notes

- No report file is included in this repository by design as requested.
- This info file provides context for reviewers and employers about the project scope and flow.
