# Use paired t-test to check if EFSIS is significantly better or worse than function perturbation

# Function to calculate the standard variance of each row
rowSd <- function(x) {
  return(apply(x, 1, sd))
}

# Load the data
data.set <- 'DLBCL'
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
auc.sd.list <- data.frame(row.names = names.methods)

for (i in percent.sel.fea) {
  method.better <- c()
  method.worse <- c()
  num.sel.fea <- ceiling(num.fea * i)  # the numbers of selected features
  file.name <- paste('../data/', data.set, '/num', num.sel.fea, '-auc-main-btsp-css.txt', sep = '')
  df <- read.table(file.name)  # read the AUC file
  num.methods <- nrow(df)
  auc.func <- df[5, ]
  auc.efsis <- df[6, ]
  if (!identical(unlist(auc.func), unlist(auc.efsis))) {
    p.value <- t.test(unlist(auc.func), unlist(auc.efsis), paired = T)$p.value
    if (p.value < 0.05) {
      print(paste('Percent:', i))
      if (mean(unlist(auc.func)) < mean(unlist(auc.efsis))) {
        print('EFSIS is significantly BETTER than Function Perturbation')
      } else {
        print('EFSIS is significantly WORSE than Function Perturbation')
      }
    }
  }
  auc.mean <- rowMeans(df)
  auc.sd <- rowSd(df)
  auc.mean.list[as.character(i)] <- auc.mean
  auc.sd.list[as.character(i)] <- auc.sd
}


