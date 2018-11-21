# Use wilcoxon Signed Test to check if there are significant differences between the best method and the others

library(foreign)
library(scmamp)
# Function to calculate the standard variance of each row
rowSd <- function(x) {
  return(apply(x, 1, sd))
}

# The function to the best method(s): max mean & min std. May be multiple
Bests <- function(df){  # df is the dataframe of the auc matrix
  ave <- rowMeans(df, na.rm = F, dims = 1)
  pos.max <- which(ave == max(ave))  # the index of the method with highest average
  ave.best <- df[pos.max, ]  # save the values of that method as another variable, may be more than one
  if (length(pos.max) > 1){
    var <- rowSd(ave.best)
    pos.min <- which(var == min(var))  # the index of the method with min standard deviation, still, may be more than one
    var.best <- ave.best[pos.min, ]
    best <- var.best
  } else {
    best <- ave.best
  }
  return(best)  # return the best ones
}

# The function to the best method(s): max mean & min std. If there are multiple, pick the first one
Best <- function(df){  # df is the dataframe of the auc matrix
  ave <- rowMeans(df, na.rm = F, dims = 1)
  pos.max <- which(ave == max(ave))  # the index of the method with highest average
  ave.best <- df[pos.max, ]  # save the values of that method as another variable, may be more than one
  if (length(pos.max) > 1){
    var <- rowSd(ave.best)
    pos.min <- which(var == min(var))  # the index of the method with min standard deviation
    var.best <- ave.best[pos.min, ]
    best <- var.best
    if (length(pos.min) > 1){
      best <- var.best[1, ]
    } else {
      best <- var.best
    }
  } else {
    best <- ave.best
  }
  return(best)  # return the best one
}

# Load the data
args <- commandArgs(TRUE)
data.set <- args[1]
script.version <- args[2]
path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
data.file <- list.files(path = path.data, pattern = '.arff')
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
percent.sel.fea <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5) / 100  # percentages of selected features

# Get the general information about this dataset

num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
# nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features

num.folds <- 10
names.methods <- c('SAM', 'GeoDE', 'RelifF', 'Information Gain', 'Function Perturbation','EFSIS')

auc.mean.list <- data.frame(row.names = names.methods)
auc.sd.list <- data.frame(row.names = names.methods)

for (i in percent.sel.fea) {
  method.worse <- c()
  num.sel.fea <- ceiling(num.fea * i)  # the numbers of selected features
  file.name <- paste('../data/', data.set, '/num', num.sel.fea, '-auc-', script.version, '.txt', sep = '')
  df <- read.table(file.name)[c(1:5, 7), ]  # read the AUC file
  num.methods <- nrow(df)
  
  # ======= if there are multiple best ones, pick the first one ====
  bests <- Bests(df)
  best <- Best(df)
  for (j in c(1:num.methods)) {  # compare the auc of each method with the best one
    auc.competitor <- df[j, ]
    if (!identical(unlist(best), unlist(auc.competitor))){
      # p.value <- wilcox.test(unlist(best), unlist(auc.competitor), paired = TRUE,
      #                        alternative = 'greater')$p.value
      p.value <- wilcoxonSignedTest(unlist(best), unlist(auc.competitor))$p.value
      if (p.value < 0.05){
          method.worse <- c(method.worse, names.methods[j])
      }
    }
  }
  
  print(paste('Percent: ', i*100, '%', sep = ''))
  print(paste('The best method(s):', row.names(bests)))
  if (length(method.worse) > 0) {
    print(paste('Worse method than the best one:', paste(method.worse, collapse = ', '))) 
  }

  auc.mean <- rowMeans(df)
  auc.sd <- rowSd(df)
  auc.mean.list[as.character(i)] <- auc.mean
  auc.sd.list[as.character(i)] <- auc.sd
}

num.row <- nrow(auc.mean.list)
num.col <- ncol(auc.mean.list)
table <- matrix(nrow = num.row, ncol = num.col)
for (i in c(1:num.row)){
  for (j in c(1:num.col)){
    table[i, j] <- paste(format(round(auc.mean.list[i, j], 2), nsmall = 2), '\u00b1', format(round(auc.sd.list[i, j], 2), nsmall = 2), '\t', collapse = ' ')
  }
}

file.name <- paste('mean-sd', data.set, script.version, '.txt', sep = '_')
write.table(table, file.name, quote = F, row.names = F, col.names = F)

