# Use paired t-test to check if there are significant differences between ensemble (function perturbation and EFSIS) and best individual (with the highest average AUC among all individuals)

# Function to calculate the variance of each row
rowVar <- function(x) {
  return(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

# The function to calculate the p-adj for each method compared with the individual method with the highest average AUC
BestIndv <- function(df, vec){  # df is the dataframe of the auc matrix, vec gives the positions of the individual method (row numbers)
  indv <- df[vec, ]  # extract the individual methods
  ave <- rowMeans(indv, na.rm = F, dims = 1)
  pos.max <- which(ave == max(ave))  # the index of the individual method with highest average
  indv.ave.best <- indv[pos.max, ]  # save the values of that method as another variable, may be more than one
  if (length(pos.max) > 1){
    var <- rowVar(indv.ave.best)
    pos.min <- which(var == min(var))  # the index of the individual method with min variance
    indv.var.best <- indv.ave.best[pos.min, ]
    if (length(pos.min) > 1){
      indv.best <- indv.var.best[1, ]
    } else {
      indv.best <- indv.var.best
    }
  } else {
    indv.best <- indv.ave.best
  }
  return(indv.best)  # return the best individual AUC list
}

# Load the data
data.set <- 'ProstateSingh'
path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
data.file <- list.files(path = path.data, pattern = '.arff')
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
percent.sel.fea <- c(0.1, 0.2, 0.4, 0.7, 1, 1.5, 2, 3, 4, 5) / 100  # percentages of selected features

# Get the general information about this dataset

num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
# nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features

num.folds <- 10
names.methods <- c('SAM', 'GeoDE', 'RelifF', 'Information Gain', 'Function Perturbation', 'EFSIS', 'CSS')
vec.indv <- c(1:4)  # the indexes of all the individual methods

auc.mean.list <- data.frame(row.names = names.methods)
auc.var.list <- data.frame(row.names = names.methods)

for (i in percent.sel.fea) {
  method.sig <- c()
  num.sel.fea <- ceiling(num.fea * i)  # the numbers of selected features
  file.name <- paste('../data/', data.set, '/num', num.sel.fea, '-auc-main-btsp-css.txt', sep = '')
  df <- read.table(file.name)  # read the AUC file
  num.methods <- nrow(df)
  indv.best <- BestIndv(df, vec.indv)
  for (j in c(1:num.methods)) {  # compare the auc of each method with the best individual one
    auc.competitor <- df[j, ]
    if (!identical(unlist(indv.best), unlist(auc.competitor))){
      p.value <- t.test(unlist(indv.best), unlist(auc.competitor), paired = T)$p.value
      if (p.value < 0.05){
        method.sig <- c(method.sig, names.methods[j])
      }
    }
  }
  if (length(method.sig) > 0){
    print(paste('Percent:', i, 'method:', method.sig))
  }
  auc.mean <- rowMeans(df)
  auc.var <- rowVar(df)
  auc.mean.list[as.character(i)] <- auc.mean
  auc.var.list[as.character(i)] <- auc.var
}

