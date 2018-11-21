# Collect the stability of Function Perturbation and EFSIS from different datasets and calculate the p-value

library(foreign)
library(reshape2)
library(ggplot2)
library(scmamp)

cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

num.folds <- 10

# Load the data
args <- commandArgs(TRUE)
data.set <- args[1]
script.version <- args[2]
# script.version <- 'efsosParalStab'
path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
data.file <- list.files(path = path.data, pattern = '.arff')
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
percent.sel.fea <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5) / 100

# Get the general information about this dataset

num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features

df.merge <- data.frame()  # to record the final merged dataframe
names.methods <- c('SAM', 'GeoDE', 'ReliefF', 'Information Gain', 'Function Perturbation','EFSIS')
stabs <- matrix(nrow = length(nums.sel.fea), ncol = length(names.methods))
j <- 0
for (i in nums.sel.fea) {
  j = j + 1
  file.name <- paste(path.data, 'num', i, '-stab-', script.version, '.txt', sep = '')
  stab <- read.table(file.name)[c(1:5, 7), ]  # read the file, skip the line for 'EFSOS'
  stabs[j, ] <- unlist(stab)
}

df.stabs <- as.data.frame.matrix(stabs)
colnames(df.stabs) <- names.methods

stabs.fp.efsis <- df.stabs[, c(5,6)]

# wilcox_test <- wilcox.test(stabs.fp.efsis$`Function Perturbation`, stabs.fp.efsis$EFSIS, paired = TRUE,
#                       alternative = 'less')
wilcox_test <- wilcoxonSignedTest(unlist(stabs.fp.efsis$`Function Perturbation`), unlist(stabs.fp.efsis$EFSIS))
print(stabs.fp.efsis$`Function Perturbation`)
print(stabs.fp.efsis$EFSIS)
print(paste('p-value', wilcox_test$p.value))
