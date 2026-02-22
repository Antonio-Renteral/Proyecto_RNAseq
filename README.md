# RNA-seq Differential Expression Analysis

## Overview

This project performs a differential gene expression analysis using RNA-seq data from Act_CD4 cells. The analysis compares genotypes to identify significantly upregulated and downregulated genes.

The workflow includes quality control, gene filtering, statistical modeling with DESeq2, and visualization of results.

## Dataset

Processed gene-level count data were obtained using the `recount3` package.

The analysis focuses on:
- Cell type: Act_CD4
- Comparison: Bap1 mutant vs wild-type (Bap1wt/wt_CD4Cre)

## Methods

- Quality control based on gene assignment proportion
- Filtering of low-count genes
- Differential expression analysis using **DESeq2**
- Variance stabilizing transformation (VST)
- Visualization:
  - PCA plot
  - MA plot
  - Volcano plot
  - Heatmap of top differentially expressed genes

## Repository Structure

Proyecto_RNAseq/
│
├── code/ # Analysis scripts
├── processed-data/ # RDS objects
├── raw-data/ # Raw input files
├── plots/ # Generated figures
└── README.md

## How to Reproduce

1. Run `DownloadData.R` to download and preprocess the dataset.
2. Run `Model.R` to perform differential expression analysis.
3. Run `Visualization.R` to generate figures.

## Author

Addiel Antonio Platas Renteral
LCG – RNA-seq Project (2026)
