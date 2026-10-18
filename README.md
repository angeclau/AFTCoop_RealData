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

1. Install the R package **AFTCoop** from GitHub (results refer to version 0.2.4).
2. Download the R script *main_TCGA.R* in this repository.
3.  Open the  *main_TCGA.R* in RStudio and set the working directory to the source file location
4.  Download the R objects containing the real data (i.e., folder **data/**) and locate the folder in your working directory.
5. Run the  *main_TCGA.R* script with the given parameter configuration.
6.  Output will be saved in the folder **Results**. Here you can find results for Tables 4-7 as described below.
7.  After the analysis is complete, run the *main_boxplot_TCGA.R* script to produce the supplementary figures.

Folder **Results** contains
- weibull_latex_table_mgpval_BRCA.txt
- weibull_latex_table_mgpval_COLORECTAL.txt
- weibull_latex_table_mgpval_GLIOMA.txt
- weibull_latex_table_mgpval_UCEC.txt

that correspond to Tables 4-7 in Section 6 of the maniscript.

## 🧪 AFTCoop Installation

For now, install the latest version of **AFTCoop** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("angeclau/AFTCoop")
```

Recall that for installing a package from GitHub you need the "devtools" library.
For more information about AFTCoop installation, see: https://github.com/angeclau/AFTCoop


## Other dependencies

This analysis workflow requires the following R packages:

``` r
library(devtools)
library(survival) 
library(Hmisc) 
library(flexsurv)
library(cvTools)
library(tidyverse) 
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(survAUC)
library(systemfonts) 
library(extrafont)
```

You can install all required packaged using

``` r
install.packages(c("survival", "Hmisc", "flexsurv", "cvTools",
                   "tidyverse", "hrbrthemes", "viridis" ,"ggpubr", "gridExtra", "survAUC",
                   "systemfonts", "extrafont", "devtools"))
```

see Session info below for more details.

## 🧪 Running AFTCoop Real case studies

The R function *main_TCGA.R* serves as a wrapper for the entire real data analysis that consist of 4 datasets.

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

The processed datasets are also available online at https://github.com/angeclau/AFTCoop_RealData/tree/main/data

The intermediate outputs are saved in the **Results** folder and consists of several boxplots (in PDF format) and several tables (in txt format).
For example the *weibull_table_manuscript_mgpval_BRCA.txt* contains the results for the BRCA.
The *Boxplot_nsplit_Cindex_UnoCindex_comb_50_weibull_BRCA.pdf* the boxplots.

## 🧪 Other useful Info about the real data analysis

The case study results shown in the manuscript were obtained using R
version 4.5.3 on a MacBook pro Apple M4 Max using only parallelization
over rho values and default settings.

Overall, the execution of **main_TCGA.R** requires about 10 hours.
With 
- TCGA_BRCA requiring about 4.5 h (for all 50 independent spits)
- TCGA_COLORECTAL requiring about 2 h (for all 50 independent spits).
- TCGA_GLIOMA requiring about 3 h (for all 50 independent spits).
- TCGA_UCEC requiring about 45 minutes (for all 50 independent spits).

## Session info


The results in the manuscript were obtained with the environment
summarized below. Click on the section to expand the full list of
package versions.

| Setting  | Value                                  |
|----------|----------------------------------------|
| R        | 4.5.3 (2026-03-11)                     |
| OS       | macOS 27.0.1 (aarch64, darwin20)       |
| IDE      | RStudio 2026.09.0+174 (desktop)        |
| Pandoc   | 3.10                                   |
| Quarto   | 1.10.18                                |
| Locale   | en_US.UTF-8                            |
| Timezone | Europe/Rome                            |
| Date     | 2026-10-16                             |

<details>
<summary><b>Package versions</b></summary>

```text
 package      * version    date (UTC) lib source
 abind          1.4-8      2024-09-12 [1] CRAN (R 4.5.0)
 AFTCoop      * 0.2.4      2026-04-28 [1] local (AFTCoop_0.2.4.tar.gz)
 backports      1.5.1      2026-04-03 [1] CRAN (R 4.5.2)
 base64enc      0.1-6      2026-02-02 [1] CRAN (R 4.5.2)
 bitops         1.1-0      2026-07-30 [1] CRAN (R 4.5.2)
 broom          1.0.13     2026-05-14 [1] CRAN (R 4.5.2)
 car            3.1-5      2026-02-03 [1] CRAN (R 4.5.2)
 carData        3.0-6      2026-01-30 [1] CRAN (R 4.5.2)
 caTools      * 1.18.4     2026-07-20 [1] CRAN (R 4.5.2)
 checkmate      2.3.4      2026-02-03 [1] CRAN (R 4.5.2)
 cli            3.6.6      2026-04-09 [1] CRAN (R 4.5.2)
 cluster        2.1.8.3    2026-07-30 [1] CRAN (R 4.5.2)
 colorspace     2.1-4      2026-10-08 [1] CRAN (R 4.5.2)
 cvTools        0.3.3      2024-03-13 [1] CRAN (R 4.5.0)
 data.table     1.18.6.1   2026-08-24 [1] CRAN (R 4.5.2)
 DEoptimR       1.2-2      2026-09-25 [1] CRAN (R 4.5.2)
 deSolve        1.42       2026-03-20 [1] CRAN (R 4.5.2)
 digest         0.6.39     2025-11-19 [1] CRAN (R 4.5.2)
 dplyr        * 1.2.1      2026-04-03 [1] CRAN (R 4.5.2)
 evaluate       1.0.5      2025-08-27 [1] CRAN (R 4.5.0)
 farver         2.1.2      2024-05-13 [1] CRAN (R 4.5.0)
 fastmap        1.2.0      2024-05-15 [1] CRAN (R 4.5.0)
 flexsurv       2.3.2      2024-08-17 [1] CRAN (R 4.5.0)
 forcats      * 1.0.1      2025-09-25 [1] CRAN (R 4.5.0)
 foreign        0.8-91     2026-01-29 [1] CRAN (R 4.5.3)
 Formula        1.2-6      2026-08-03 [1] CRAN (R 4.5.2)
 generics       0.1.4      2025-05-09 [1] CRAN (R 4.5.0)
 ggplot2      * 4.0.3      2026-04-22 [1] CRAN (R 4.5.2)
 ggpubr       * 1.0.0      2026-07-06 [1] CRAN (R 4.5.2)
 ggsignif       0.6.4      2022-10-13 [1] CRAN (R 4.5.0)
 glue           1.8.1      2026-04-17 [1] CRAN (R 4.5.2)
 gridExtra      2.3.1      2026-06-25 [1] CRAN (R 4.5.2)
 gtable         0.3.6      2024-10-25 [1] CRAN (R 4.5.0)
 Hmisc        * 5.3-0      2026-09-06 [1] CRAN (R 4.5.2)
 hms            1.1.4      2025-10-17 [1] CRAN (R 4.5.0)
 htmlTable      2.5.0      2026-04-22 [1] CRAN (R 4.5.2)
 htmltools      0.5.9      2025-12-04 [1] CRAN (R 4.5.2)
 htmlwidgets    1.6.4      2023-12-06 [1] CRAN (R 4.5.0)
 knitr          1.52       2026-09-06 [1] CRAN (R 4.5.2)
 lattice        0.23-1     2026-08-12 [1] CRAN (R 4.5.2)
 lifecycle      1.0.5      2026-01-08 [1] CRAN (R 4.5.2)
 lubridate    * 1.9.5      2026-02-04 [1] CRAN (R 4.5.2)
 magrittr       2.0.5      2026-04-04 [1] CRAN (R 4.5.2)
 Matrix         1.7-6      2026-07-25 [1] CRAN (R 4.5.2)
 mstate         0.3.3      2024-07-11 [1] CRAN (R 4.5.0)
 muhaz          1.2.6.5    2026-07-29 [1] CRAN (R 4.5.2)
 mvtnorm        1.4-2      2026-07-12 [1] CRAN (R 4.5.2)
 nnet           7.3-21     2026-08-03 [1] CRAN (R 4.5.2)
 numDeriv       2016.8-1.1 2019-06-06 [1] CRAN (R 4.5.0)
 otel           0.2.0      2025-08-29 [1] CRAN (R 4.5.0)
 pillar         1.11.1     2025-09-17 [1] CRAN (R 4.5.0)
 pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.5.0)
 purrr        * 1.2.2      2026-04-10 [1] CRAN (R 4.5.2)
 quadprog       1.5-8      2019-11-20 [1] CRAN (R 4.5.0)
 R6             2.6.1      2025-02-15 [1] CRAN (R 4.5.0)
 RColorBrewer   1.1-3      2022-04-03 [1] CRAN (R 4.5.0)
 Rcpp           1.1.2      2026-07-05 [1] CRAN (R 4.5.2)
 readr        * 2.2.0      2026-02-19 [1] CRAN (R 4.5.2)
 rlang          1.3.0      2026-07-05 [1] CRAN (R 4.5.2)
 rmarkdown      2.32       2026-09-01 [1] CRAN (R 4.5.2)
 robustbase     0.99-7     2026-02-05 [1] CRAN (R 4.5.2)
 rpart          4.1.27     2026-03-27 [1] CRAN (R 4.5.2)
 RSpectra       0.16-2     2024-07-18 [1] CRAN (R 4.5.0)
 rstatix        1.1.0      2026-07-23 [1] CRAN (R 4.5.2)
 rstudioapi     0.19.0     2026-06-11 [1] CRAN (R 4.5.2)
 S7             0.2.2      2026-04-22 [1] CRAN (R 4.5.2)
 scales         1.4.0      2025-04-24 [1] CRAN (R 4.5.0)
 sessioninfo    1.2.4      2026-06-04 [1] CRAN (R 4.5.2)
 statmod        1.5.2      2026-05-17 [1] CRAN (R 4.5.2)
 stringi        1.8.9      2026-08-04 [1] CRAN (R 4.5.2)
 stringr      * 1.6.0      2025-11-04 [1] CRAN (R 4.5.0)
 survAUC      * 1.4-0      2025-09-26 [1] CRAN (R 4.5.0)
 survival     * 3.8-12     2026-09-09 [1] CRAN (R 4.5.2)
 tibble       * 3.3.1      2026-01-11 [1] CRAN (R 4.5.2)
 tidyr        * 1.3.2      2025-12-19 [1] CRAN (R 4.5.2)
 tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.5.0)
 tidyverse    * 2.0.0      2023-02-22 [1] CRAN (R 4.5.0)
 timechange     0.4.0      2026-01-29 [1] CRAN (R 4.5.2)
 tzdb           0.5.0      2025-03-15 [1] CRAN (R 4.5.0)
 vctrs          0.7.3      2026-04-11 [1] CRAN (R 4.5.2)
 viridis      * 0.6.5      2024-01-29 [1] CRAN (R 4.5.0)
 viridisLite  * 0.4.3      2026-02-04 [1] CRAN (R 4.5.2)
 withr          3.0.3      2026-06-19 [1] CRAN (R 4.5.2)
 xfun           0.61       2026-09-16 [1] CRAN (R 4.5.2)

 [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
 * ── Packages attached to the search path.
```

</details>
## 📚 Citation
please cite:
1. C. Angelini, D. De Canditiis, I. De Feis, A. Iuliano. *Cooperative AFT models for multi-omics data integration* in preparation (2025)

### 🏛 Funding
This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.
