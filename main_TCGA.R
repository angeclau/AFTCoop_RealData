# This function run the real data example for the manuscript

#----------------------------------------------------------
# Real case - TCGA data
#----------------------------------------------------------

# Set the output directory
dir_out <- paste0(getwd(),"/Results/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out)
}

#----------------------------------------------------------
# Load required packages
#----------------------------------------------------------

myseed <- 123
set.seed(myseed)

library(AFTCoop)
library(caTools)
library(Hmisc)
library(survival)
library(survAUC)
library(tidyverse)
library(viridis)
library(ggpubr)

#----------------------------------------------------------
# Set input parameters
#----------------------------------------------------------

model <- "weibull"                 # "weibull", "lognormal", "loglogistic" 
topn_U <-  500                     # top-ranked variables for U
topn_Z <-  500                     # top-ranked variables for Z
rank <- "mgpval"                   # screening type - "absmg", "mg", "mgpadj", "mgpval"
nfolds <- 5                        # number of folds
lam_min <- T                       # T: CV chooses minimum, F: CV chooses 1sd
rho_values <- c(1,0.25,0.5,0.75)   # agreement parameter
nlambda <- 100                     # Number of lambda values to consider in the regularization path
lambda.ratio.min <- 0.01           # Minimum ratio of lambda max to lambda min
parallel.rho <- TRUE               # Whether to use parallel computing for rho values
parallel.cv <- TRUE                # Whether to use parallel computing for cross-validation
ncore_max_rho <- 4                 # Maximum number of cores to use across different rho values
ncore_max_cv <- 5                  # Maximum number of cores to use during cross-validation
iplot <- FALSE                     # If TRUE, generate CV plot

# Note on Mac Os with Apple Silicon processor and R 4.5.2 it is better to avoid nested parallelism)
if (getRversion()>= "4.5.2" && R.version$arch =="aarch64"){ 
  print("Removed the nested parallelism") 
  parallel_cv <- FALSE
  ncore_max_cv <- 1
}

#----------------------------------------------------------
# Load data: U:= RNASeq;Z:= Methylation; Y:= (T,delta) Surv 
#----------------------------------------------------------

dataset_files <- c(
  BRCA = "data/TCGA_BRCA_2Data.RData",
  COLORECTAL = "data/TCGA_COLORECTAL_2Data.RData",
  GLIOMA = "data/TCGA_GLIOMA_2Data.RData",
  UCEC = "data/TCGA_UCEC_2Data.RData"
)

tcga_data <- lapply(dataset_files, function(f) {
  env <- new.env()
  load(f, envir = env)
  as.list(env)
})

#----------------------------------------------------------
# Run this for each TCGA dataset
#----------------------------------------------------------

for (tumor in names(tcga_data)){
  
  cat("\n====================\n")
  cat("Dataset:", tumor, "\n")
  cat("====================\n")
  
  U <- tcga_data[[tumor]]$U
  Z <- tcga_data[[tumor]]$Z
  Y <- tcga_data[[tumor]]$Y
  
  # Set argument "time" in survAUC::UnoC(), i.e. tau constant
  unotime <- min(quantile(Y$time,0.90),max(Y$time[Y$status!=0]))
  
  cat("Patient numbers: ", nrow(U), "\n")
  cat("RNASeq variables: ", ncol(U), "\n")
  cat("Methy variables: ", ncol(Z), "\n")
  cat("Perc CR%: ", length(which(Y$status==0))/nrow(Y), "\n")
  
  #---------------------------------------------------------
  # Variable screening
  #---------------------------------------------------------
  
  message("----------------------------------")
  
  # Indices of the top-ranked features in the initial set of the two views
  index_u <- variable_screening(U, Y$time, Y$status, family = "AFT", model = model, rank = rank, topn = topn_U)
  cat("Screened variables in U: ", length(index_u), "\n")
  
  index_z <- variable_screening(Z, Y$time, Y$status, family = "AFT", model = model, rank = rank, topn = topn_Z)
  cat("Screened variables in Z: ", length(index_z), "\n")
  
  # --- Save ---
  save(index_u,index_z, file = paste0(dir_out, model,"_screening_", rank, "_", tumor, ".RData"))
  
  #---------------------------------------------------------
  # Run data analysis on several splitted realizations of training/test
  #----------------------------------------------------------
  
  nsplits <- 50                                   # number of splits to run
  coeff_U <- matrix(0, length(index_u), nsplits)  # matrix to store coefficients for U
  coeff_Z <- matrix(0, length(index_z), nsplits)  # matrix to store coefficients for Z
  coeff_coop <- list()                            # list to store coefficients for coop
  
  # initialize vectors to store performance metrics across solitting runs
  C_index_U <- NULL 
  C_index_uno_U <- NULL
  C_index_Z <- NULL
  C_index_uno_Z <- NULL
  C_index_coop <- matrix(0, nsplits, length(rho_values))
  colnames(C_index_coop) <- paste0("rho_", rho_values)
  C_index_uno_coop <- matrix(0, nsplits, length(rho_values))
  colnames(C_index_uno_coop) <- paste0("rho_", rho_values)
  
  for(isplit in c(1:nsplits)){
    
    set.seed(myseed + isplit)
    
    message("----------------------------------")
    cat("Split run: ", isplit, "\n")
    
    # Split the data into 80% training and 20% testing each time
    train_ind <- sample(seq_len(nrow(Y)), size = floor(0.8 * nrow(Y)))
    
    # U, Z, and Y have patients as rows (listed in the same order)
    ## Define the training set 
    U_train <- U[train_ind, ] 
    Z_train <- Z[train_ind, ] 
    Y_train <- Y[train_ind, ] 
    
    ## Define the test set 
    U_test <- U[-train_ind, ]
    Z_test <- Z[-train_ind, ]
    Y_test <- Y[-train_ind, ]
    
    # Percentage of the censored data on training and testing set
    # cat("Perc cens train: ", length(which(Y_train$status==0))/nrow(Y_train), "\n")
    # cat("Perc cens test: ", length(which(Y_test$status==0))/nrow(Y_test), "\n")
    # Percentage of the events on training and testing set
    # cat("Perc events train: ",length(which(Y_train$status==1))/nrow(Y_train), "\n")
    # cat("Perc events test: ",length(which(Y_test$status==1))/nrow(Y_test), "\n")
    
    #----------------------------------------------------------   
    # Training phase 
    #----------------------------------------------------------
    
    # sigma estimation using survreg {survival}
    fit_survreg <- survreg(Surv(Y_train$time,Y_train$status) ~ 1, dist = model, scale = 0)
    sigma.est <- exp(fit_survreg$icoef[2]) # scale parameter
    #cat("est sigma only intercept",sigma.est, "\n")
    
    #----------------------------------------------------------
    
    case <- "onlyU" # i,e., only gene expression
    beta_Utrain_screen <- aft_coop(U = U_train[,index_u], Z = Z_train[,index_z],
                                   Y = log(Y_train$time), delta = Y_train$status,
                                   sigma = sigma.est, nfolds = nfolds, model = model, case = case, 
                                   rho_values = rho_values[1], lam_min = lam_min,
                                   nlambda = nlambda, lambda.ratio.min = lambda.ratio.min, 
                                   parallel.rho = parallel.rho, parallel.cv = parallel.cv, 
                                   ncore_max_rho = ncore_max_rho, ncore_max_cv = ncore_max_cv, 
                                   seed = myseed, iplot = iplot)
    
    coeff_U[,isplit] <- beta_Utrain_screen
    rownames(coeff_U) <- colnames(U_train[,index_u])
    
    #----------------------------------------------------------
    
    # predicted linear scores
    pred_Utest_screen <- predict_aft_coop(U = U_test[,index_u], Z = Z_test[,index_z], 
                                          mU = colMeans(U_train[,index_u]), 
                                          beta = beta_Utrain_screen[,1], 
                                          case = case)
    
    # Harrell's C-index 
    c.index_U <- rcorr.cens(pred_Utest_screen, Surv(Y_test$time, Y_test$status))
    cat('C-index Only U = ', as.numeric(c.index_U[1]), "\n")
    
    # C-statistic by Uno et al.
    c.index_uno_U <- UnoC(Surv(Y_train$time, Y_train$status), Surv(Y_test$time, Y_test$status), -pred_Utest_screen, unotime)
    cat('Uno C-index Only U = ', as.numeric(c.index_uno_U), "\n")
    
    C_index_U  <- c(C_index_U, as.numeric(c.index_U[1]))
    C_index_uno_U  <- c(C_index_uno_U, as.numeric(c.index_uno_U))
    
    #----------------------------------------------------------
    
    case <- "onlyZ"  # i,e., only methylation
    beta_Ztrain_screen <- aft_coop(U = U_train[,index_u], Z = Z_train[,index_z],
                                   Y = log(Y_train$time), delta = Y_train$status,
                                   sigma = sigma.est, nfolds = nfolds, model = model, case = case, 
                                   rho_values = rho_values[1], lam_min = lam_min,
                                   nlambda = nlambda, lambda.ratio.min = lambda.ratio.min, 
                                   parallel.rho = parallel.rho, parallel.cv = parallel.cv, 
                                   ncore_max_rho = ncore_max_rho, ncore_max_cv = ncore_max_cv, 
                                   seed = myseed, iplot = iplot)
    
    coeff_Z[,isplit] <- beta_Ztrain_screen
    rownames(coeff_Z) <- colnames(Z_train[,index_z])
    
    #----------------------------------------------------------
    
    # predicted linear scores
    pred_Ztest_screen <- predict_aft_coop(U_test[,index_u], Z_test[,index_z], 
                                          mZ = colMeans(Z_train[,index_z]), 
                                          beta = beta_Ztrain_screen[,1], 
                                          case = case)
    
    # Harrell's C-index
    c.index_Z <- rcorr.cens(pred_Ztest_screen, Surv(Y_test$time, Y_test$status))
    cat('C-index Only Z = ', as.numeric(c.index_Z[1]), "\n")
    
    # C-statistic by Uno et al. 
    c.index_uno_Z <- UnoC(Surv(Y_train$time, Y_train$status), Surv(Y_test$time, Y_test$status), -pred_Ztest_screen, unotime)
    cat('Uno C-index Only Z = ', as.numeric(c.index_uno_Z), "\n")
    
    C_index_Z  <- c(C_index_Z, as.numeric(c.index_Z[1]))
    C_index_uno_Z  <- c(C_index_uno_Z, as.numeric(c.index_uno_Z))
    
    #----------------------------------------------------------
    
    case <- "coop"  # i,e.,  Gene expression + Methylation
    beta_cooptrain_screen <- aft_coop(U = U_train[,index_u], Z = Z_train[,index_z],
                                      Y = log(Y_train$time), delta = Y_train$status,
                                      sigma = sigma.est, nfolds = nfolds, model = model, case = case, 
                                      rho_values = rho_values, lam_min = lam_min,
                                      nlambda = nlambda, lambda.ratio.min = lambda.ratio.min, 
                                      parallel.rho = parallel.rho, parallel.cv = parallel.cv, 
                                      ncore_max_rho = ncore_max_rho, ncore_max_cv = ncore_max_cv, 
                                      seed = myseed, iplot = iplot)
    
    coeff_coop[[paste0("Split_", isplit)]] <- beta_cooptrain_screen
    
    #----------------------------------------------------------   
    
    # predicted linear scores
    for(r in 1:length(rho_values)){
      pred_cooptest_screen <- predict_aft_coop(U_test[,index_u], Z_test[,index_z], 
                                               mU = colMeans(U_train[,index_u]), 
                                               mZ = colMeans(Z_train[,index_z]),
                                               beta = beta_cooptrain_screen[,r], 
                                               case = case)
      
      # Harrell's C-index
      c.index_coop <- rcorr.cens(pred_cooptest_screen, Surv(Y_test$time, Y_test$status))
      cat('C-index Coop = ', as.numeric(c.index_coop[1]), "\n")
      
      # C-statistic by Uno et al. 
      c.index_uno_coop <- UnoC(Surv(Y_train$time, Y_train$status), Surv(Y_test$time, Y_test$status), -pred_cooptest_screen, unotime)
      cat('Uno C-index Coop = ', as.numeric(c.index_uno_coop), "\n")
      
      C_index_coop[isplit, r] <- as.numeric(c.index_coop[1])
      C_index_uno_coop[isplit, r] <- as.numeric(c.index_uno_coop)
      
    }
    
    #----------------------------------------------------------   
    
    message("----------------------------------")
    
    # Save regression coefficient
    df_U <- data.frame(Symbol = colnames(U_train[,index_u]), coefficients = beta_Utrain_screen)
    cat('Number of selected coefficient in Only U = ', sum(abs(df_U$coefficients)>0), "\n")
    
    df_Z <- data.frame(Symbol = colnames(Z_train[,index_z]), coefficients = beta_Ztrain_screen)
    cat('Number of selected coefficient in Only Z = ', sum(abs(df_Z$coefficients)>0), "\n")
    
    colnames(beta_cooptrain_screen) <- rho_values
    df_coop <- data.frame(Symbol = c(colnames(U_train[,index_u]),colnames(Z_train[,index_z])), coefficients = beta_cooptrain_screen)
    cat('Number of selected coefficient in Coop rho_1 = ', sum(abs(df_coop[,2])>0), "\n") # Early fusion
    cat('Number of selected coefficient in Coop rho_0.25 = ', sum(abs(df_coop[,3])>0), "\n") # rho=0.25
    cat('Number of selected coefficient in Coop rho_0.50 = ', sum(abs(df_coop[,4])>0), "\n") # rho=0.50
    cat('Number of selected coefficient in Coop rho_0.75 = ', sum(abs(df_coop[,5])>0), "\n") # rho=0.75
    
  }
  
  message("----------------------------------")
  
  # --- C-index ---
  df_C <- data.frame(C_index_U, C_index_Z, C_index_coop)
  colnames(df_C) <- c("only U","only Z","early fusion","rho=0.25","rho=0.5","rho=0.75")
  mean_C_index <- apply(df_C, 2, mean) 
  df_C_index <- data.frame(
    case = c("only U","only Z","early fusion","rho=0.25","rho=0.5","rho=0.75"),
    mean = mean_C_index
  )
  rownames(df_C_index) <- 1:6
  cat("Average C-index: \n",
      paste0(df_C_index$case, " = ", round(df_C_index$mean, 4), collapse = " | "),
      "\n")
  
  # --- Save ---
  file_name <- paste0(dir_out, model,"_Cindex_nsplits_", rank, "_", tumor, ".txt")
  write.table(df_C, file = file_name, sep = "\t", row.names = F, col.names = T)
  
  #---------------------------------------------------------
  
  # --- Uno-C-index ---
  df_uno_C <- data.frame(C_index_uno_U, C_index_uno_Z, C_index_uno_coop)
  colnames(df_uno_C) <- c("only U","only Z","early fusion","rho=0.25","rho=0.5","rho=0.75")

    mean_uno_C_index <- apply(df_uno_C, 2, mean) 
  df_uno_C_index <- data.frame(
    case =  c("only U","only Z","early fusion","rho=0.25","rho=0.5","rho=0.75"),
    mean = mean_uno_C_index
  )
  rownames(df_uno_C_index) <- 1:6
  cat("Average Uno C-index: \n",
      paste0(df_uno_C_index$case, " = ", round(df_uno_C_index$mean, 4), collapse = " | "),
      "\n")
  
  # --- Save ---
  file_name <- paste0(dir_out, model,"_CUnoindex_nsplits_", rank, "_", tumor, ".txt")
  write.table(df_uno_C, file = file_name, sep = "\t", row.names = F, col.names = T)
  
  colnames(coeff_U) <- paste0("Split_", 1:nsplits)
  colnames(coeff_Z) <- paste0("Split_", 1:nsplits)
  coeff_coop <- lapply(coeff_coop, function(df) {
    colnames(df) <- paste0("rho_", rho_values)
    rownames(df) <- c(rownames(coeff_U), rownames(coeff_Z))
    df
  })
  
  # --- Save beta coefficients ---
  file_out <- paste0(dir_out, model,"_beta_", rank, "_", tumor, ".RData")
  save(coeff_U, coeff_Z, coeff_coop, file = file_out)
  
  message("----------------------------------")
  
  #-------------------------------------------
  # Shared genes - Voting (set the desired thresholed e.g., 70%)
  #-------------------------------------------
  voting_thresh <- 0.70
  
  # --- Only U ---
  perc_U <- rowSums(coeff_U != 0) / ncol(coeff_U)
  selected_genes_U <- rownames(coeff_U)[perc_U >=  voting_thresh]
  NSR_U <- length(selected_genes_U)/topn_U
  cat("Selected genes Only U (",voting_thresh,"%):", NSR_U, "\n")
  
  # --- Only Z ---
  perc_Z <- rowSums(coeff_Z != 0) / ncol(coeff_Z)
  selected_genes_Z <- rownames(coeff_Z)[perc_Z >=  voting_thresh]
  NSR_Z <- length(selected_genes_Z)/topn_Z
  cat("Selected genes Only Z (",voting_thresh,"%):", NSR_Z, "\n")
  
  # --- Coop ---
  selected_genes_coop <- list()
  NSR_coop <- NULL
  nselected_genes_coop<-vector("integer", length = length(rho_values)) 
  for (i in 1:length(rho_values)) {
    combined_col_matrix <- sapply(coeff_coop, function(sim) sim[, i])
    
    freq_selection <- rowSums(combined_col_matrix != 0) / length(coeff_coop)
    genes_over_threshold <- rownames(coeff_coop[[1]])[freq_selection >= voting_thresh]
    rho_label <- paste0("rho_", rho_values[i])
    selected_genes_coop[[rho_label]] <- genes_over_threshold
    nselected_genes_coop[i]<-length(genes_over_threshold)
    NSR <- length(genes_over_threshold)/(topn_U+topn_Z)
    NSR_coop <- rbind(NSR_coop, NSR)
    cat("rho =", rho_values[i], "| Selected genes (",voting_thresh,"%) :", NSR, "\n")
  }
  # Voting genes (%)
  NSR <- c(NSR_U, NSR_Z, NSR_coop)
  
  # --- Table in the manuscript for each dataset ---
  table_manuscript <- data.frame(Case=c("only U","only Z","early fusion","rho=0.25","rho=0.5","rho=0.75"), 
                                 "C_index"= df_C_index$mean, "Uno_C-_ndex" = df_uno_C_index$mean, NSR)
  
  # --- Save table data ---
  file_table <- paste0(dir_out, model,"_table_manuscript_", rank, "_", tumor, ".txt")
  write.table(table_manuscript, file = file_table, sep = "\t", row.names = F, col.names = T)
  
  # --- Save selected genes name ---
  file_out <- paste0(dir_out, model,"_selected_genes_", rank, "_", tumor, "_voting",voting_thresh,"%", ".RData")
  save(selected_genes_U, selected_genes_Z, selected_genes_coop, file = file_out)
  
  # --- latex table ---
  table_latex <- paste0(
    "\\begin{tabular}{lr}\n",
    "\\hline\n",
    "Case & C-index & Uno C-index & NSR (",voting_thresh,"%) \\\\\n",
    "\\hline\n",
    paste0(
      df_uno_C_index$case, " & ", round(df_C_index$mean, 4), " & ", round(df_uno_C_index$mean, 4), "& ", NSR, " \\\\",
      collapse = "\n"
    ),
    "\n\\hline\n",
    "\\end{tabular}"
  )
  
  cat(table_latex)
  
  # --- Save latex table ---
  file_res <- paste0(dir_out, model,"_latex_table_", rank, "_", tumor, ".txt")
  cat(table_latex, file = file_res, append = FALSE)
  
}
