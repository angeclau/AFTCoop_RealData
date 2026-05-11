---
title: "AFTCoop Real Case Studies"
output: rmarkdown::github_document
---

# AFTCoop Applied to TCGA Cancer Data

This repository contains the scripts required to reproduce the real-data analyses presented in the manuscript: C. Angelini, D. De Canditiis, I. De Feis, and A. Iuliano. *Cooperative AFT models for multi-omics data integration* submitted (2025).

The results refer to Section 6 (Tables 4-7) and to the Supplementary information file (Section S5 and Figures S5-S8).

The examples illustrate the use of the **AFTCoop** R package for fitting cooperative Accelerated Failure Time (AFT) survival regression models integrating two omics data views:

- RNA-seq gene expression (that constitutes matrix **U**)
- DNA methylation (that constitutes matrix **Z**)

## Description

To execute the analysis of the case studies:

-  Install the R package **AFTCoop** from GitHub (results refer to version 0.2.4).
-  Download the R script *main_TCGA.R* in this repository.
-  Open the  *main_TCGA.R* in RStudio and set the working directory to the source file location
-  Download the R objects containing the real data (i.e., folder **data/**) and locate the folder in your working directory.
-  Run the  *main_TCGA.R* script with the given parameter configuration.
-  Output will be saved in the folder **Results**
-  After the aalysis is complete, run the *main_boxplot_TCGA.R* script to produce the figures.

## 🧪 AFTCoop Installation

For now, install the latest version of **AFTCoop** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("angeclau/AFTCoop")
```

For more information see: https://github.com/angeclau/AFTCoop.

## 🧪 Running AFTCoop simulation study

The R function *main_TCGA.R* serves as a wrapper for the entire real data analysis that consist of 4 dataset.

Processed data are available in folder **data/**  (retrieved and assembled from linkedomics (https://www.linkedomics.org/data_download/).

In the folder **data/** there are 4 datasets:

- TCGA_BRCA_2Data.RData  
- TCGA_COLORECTAL_2Data.RData  
- TCGA_GLIOMA_2Data.RData  
- TCGA_UCEC_2Data.RData

For each dataset:

- **U:** normalized RNA-seq expression matrix RSEM counts transformed as log2(value + 1) generated from the Illumina HiSeq platform
- **Z**: DNA methylation matrix obtained from the Illumina HumanMethylation450K (HM450K) platform represented as centered Beta-values (value - 0.5) aggregated at the gene level
- **Y**: survival outcome information (overall survival time and censoring indicator).

The pre-computed output is saved in the *Results** folder and consists of several boxplots (in PDF format) and several tables (in txt format).
For example the *weibull_table_manuscript_mgpval_BRCA.txt* contains the results for the BRCA.
The *Boxplot_nsplit_Cindex_UnoCindex_comb_50_weibull_BRCA.pdf* the boxplots.


## 📚 Citation
please cite:
1. C. Angelini, D. De Canditiis, I. De Feis, A. Iuliano. *Cooperative AFT models for multi-omics data integration* in preparation (2025)

### 🏛 Funding
This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.
