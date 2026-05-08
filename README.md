---
title: "AFTCoop real data examples"
output: rmarkdown::github_document
---

# AFTCoop applied to cancer data from TGCA

These scripts enable the execution of the real data example described in the manuscript by C. Angelini, D. De Canditiis, I. De Feis, and A. Iuliano. *Cooperative AFT models for multi-omics data integration* submitted (2025).
 
Such an example aims to illustrate the use of the R package **AFTCoop** in fitting cooperative AFT survival regression models with two omics views, corresponding to RNA-seq and DNA Methylation.

## Description

To execute the simulation:

-  Install the R package AFTCoop from GitHub.
-  Download the R script main_TCGA.R in this repository.
-  Download the R object containing the real data (i.e., folder data/).
-  Open the  main_TCGA.R in RStudio and set the working directory to the source file location.
-  Run the  main_TCGA.R script with the given parameter configuration.
-  Output will be saved in the folder Results
-  Run the main_boxplot_TCGA.R script to produce the figures.

## 🧪 AFTCoop Installation

For now, install the latest version of **AFTCoop** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("angeclau/AFTCoop")
```

## 🧪 Running AFTCoop simulation study

The R function main_TCGA.R serves as a wrapper for the entire real data analysis.

Processed data are available in fodet **data/**  (retrieved from linkedomics (\url{https://www.linkedomics.org/data_download/}).

- Matrix **U** corresponds to normalized RNA-seq count data (RSEM, $log_2(\text{Val} + 1)$) from the Illumina HiSeq platform; 

- Matrix **Z** corresponds to DNA methylation obtained from the Illumina HumanMethylation450K (HM450K) platform and reported as Beta-values (centered as value-$0.5$) and aggregated at the gene level. 

- **Y** corresponds to the overall survival time and censoring indicator.

## 📚 Citation
please cite:
1. C. Angelini, D. De Canditiis, I. De Feis, A. Iuliano. *Cooperative AFT models for multi-omics data integration* in preparation (2025)

### 🏛 Funding
This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.
