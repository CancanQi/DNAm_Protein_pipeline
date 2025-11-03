rm(list=ls())
# load package

suppressPackageStartupMessages(library(optparse))

library(foreign)
library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huber??s estimation of the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(R.utils)
library(matrixStats)
library(plyr)

# defination of args
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input metadata file path", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output results file path", metavar = "FILE"),
  make_option(c("-c", "--cols"), type = "character", help = "Columns to extract (comma-separated, X or Y should be the first)", metavar = "COLS"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print progress messages [default %default]")
)

# get args
opt <- parse_args(OptionParser(option_list = option_list))

# args checking
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$cols)) {
  cat("❌ Error: You must specify input (-i), output (-o) and columns (-c).\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

# confirm args
if (opt$verbose) {
  cat("✅ Parameters received:\n")
  cat("Input file: ", opt$input, "\n")
  cat("Output file:", opt$output, "\n")
  cat("Columns:    ", opt$cols, "\n\n")
}

# load metadata
df.meta <- read.csv(opt$input, stringsAsFactors = FALSE,row.names = 1)
# load DNAm change path to right one
Mvals<-readRDS("/share/users/qicancan/DNAmethylation/demo/DNAm.beta.demo.RDS")

# get cols to extract from args
cols_to_extract <- unlist(strsplit(opt$cols, split = ","))

# check cols existance
missing_cols <- setdiff(cols_to_extract, colnames(df.meta))
if (length(missing_cols) > 0) {
  stop(paste("❌ The following columns do not exist in the data:", paste(missing_cols, collapse = ", ")))
}

### match samples
keep<-intersect(df.meta$DNAmID,colnames(Mvals))
print(paste0("N of matched samples: ",length(keep)))

### prepare input data
M_matrix<-t(Mvals[,match(keep,colnames(Mvals))])
PHENO<-df.meta[,cols_to_extract]

### add GLM function
LMtest <- function(meth_matrix, methcol, Y, covariates = NULL) {
  # extract specific cpg site
  probe_vals <- meth_matrix[, methcol]
  
  # reformat
  df <- data.frame(meth = probe_vals, Y = Y)
  
  if (!is.null(covariates)) {
    df <- cbind(df, covariates)
  }
  
  # remove NA
  df_complete <- df[complete.cases(df), ]
  n_samples <- nrow(df_complete)
  
  # generate formula
  if (!is.null(covariates)) {
    formula_str <- paste("meth ~ Y + ", paste(colnames(covariates), collapse = "+"))
  } else {
    formula_str <- "meth ~ Y"
  }
  
  # lm fit or glm fit
  mod <- lm(as.formula(formula_str), data = df_complete)
  cf <- summary(mod)$coefficients
  
  # returen probeID, Estimate, Std.Error, P-value, sample size
  data.frame(
    probeID = colnames(meth_matrix)[methcol],
    BETA = cf["Y", "Estimate"],
    SE = cf["Y", "Std. Error"],
    P_VAL = cf["Y", "Pr(>|t|)"],
    N = n_samples
  )
}

# get data frame of all covariates
covariates <- PHENO[, 2:ncol(PHENO)]   # pay attention to the columns!

# run the analysis 
system.time(
  ind.res <- mclapply(
    seq_len(ncol(M_matrix)),
    LMtest,
    meth_matrix = M_matrix,
    Y = PHENO[,1],
    covariates = covariates
  )
)

# merge results
all.results <- ldply(ind.res, rbind)

# calculate lambda value 
lambda<-median(qchisq(as.numeric(as.character(all.results$P_VAL)),df=1,lower.tail = F),na.rm = T)/qchisq(0.5,1)
print(paste0("lambda value is ",round(lambda,digits = 3)))

write.table(all.results, opt$output,na="NA",row.names = FALSE)
gzip(opt$output)
