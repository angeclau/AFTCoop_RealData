
# This function create the C-index e Uno C-index boxplot for the real studies in the manuscript

# Load packages ---
library(tidyverse)
library(ggpubr)
library(gridExtra)

# Set the output directory
dir_out <- paste0(getwd(),"/Results/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out)
}

source("boxplot_measure.R") ## for generating manuscript figures (boxplots)

# Set parameters for plot 
tumor  <- c("BRCA", "COLORECTAL", "GLIOMA", "UCEC")
model  <- "weibull"
rank   <- "mgpval"
nsplit <- 50

for(i in 1:length(tumor)){
  
  file_name <- paste0(dir_out, model,"_Cindex_nsplits_", rank, "_", tumor[i], ".txt")
  Cindex <- read.table(file = file_name, sep = "\t",header = TRUE,check.names=FALSE)
  bxp_Cindex <- boxplot_measure(nsplit, Cindex, "Harrell's C-index", c(0.4,1))
  
  file_name <- paste0(dir_out, model,"_CUnoindex_nsplits_", rank, "_", tumor[i], ".txt")
  UnoCindex <- read.table(file = file_name, sep = "\t",header = TRUE,check.names=FALSE)
  bxp_UnoCindex <- boxplot_measure(nsplit, UnoCindex, "Uno C-index", c(0.4,1))
  
  # Reduce spacing between panels
  bxp_Cindex <- bxp_Cindex + theme(plot.margin = margin(5,2,5,2))
  bxp_UnoCindex <- bxp_UnoCindex + theme(plot.margin = margin(5,2,5,2))
  
  # Save combined boxplot
  file_out <- paste0(dir_out, "Boxplot_nsplit_Cindex_UnoCindex_comb_", nsplit, "_", model, "_", tumor[i], ".pdf")
  combined_plot2 <- grid.arrange(bxp_Cindex, bxp_UnoCindex, ncol = 2, nrow = 1)
  ggsave(file_out, plot = combined_plot2, width = 12, height = 4, device = cairo_pdf)
  
}
