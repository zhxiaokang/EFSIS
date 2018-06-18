# plot the boxplot for auc from efsis_diff_10fold.R

library(reshape2)
library(ggplot2)

cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

num.folds <- 10
# Load the data
data.set <- 'AML'
path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
data.file <- list.files(path = path.data, pattern = '.arff')
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
percent.sel.fea <- c(0.1, 0.2, 0.4, 0.7, 1, 1.5, 2, 3, 4, 5) / 100

# Get the general information about this dataset

num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features

df.merge <- data.frame()  # to record the final merged dataframe

for (i in nums.sel.fea) {
  file.name <- paste('../data/', data.set, '/num', i, '-auc-main-btsp-css.txt', sep = '')
  df <- read.table(file.name)  # read the file
  # set the row names as one column
  df$method <- row.names(df)
  # add the number of selected features as one column
  # df$num <- i
  df.melt <- melt(df, id.vars = 'method')
  # merge the current dataframe to the current merged one
  df.merge <- cbind.all(df.merge, df.melt[3])
}
methods <- row.names(df)
colnames(df.merge) <- nums.sel.fea
df.merge <- as.data.frame(df.merge)
df.merge$method <- c(rep(methods, num.folds))
df.melt <- melt(df.merge, id.vars = 'method')
df.melt$method <- factor(df.melt$method, levels = methods)

# ggplot(data = df.melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=method))

p1 <- ggplot(data = df.melt, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=method))
p1 <- p1 + facet_wrap( ~ variable, scales="free")

# get the average rankings (averaging on 10 folds)
rank.method <- data.frame(matrix(NA, nrow = length(methods), ncol = num.folds))
rownames(rank.method) <- methods
count <- 0
for (i in nums.sel.fea) {
  count <- count + 1
  file.name <- paste('../data/', data.set, '/num', i, '-auc-main-btsp-css.txt', sep = '')
  df <- read.table(file.name)  # read the file
  df.rank <- df
  for (j in ncol(df)){
    df.rank[, j] <- rank(-df[, j])
  }
  rank.method[, count] <- rowMeans(df.rank)
  names(rank.method)[count] <- paste(i, '_features')
}
rank.method.df <- as.data.frame(t(rank.method))
rank.method.df$num_sel_fea <- nums.sel.fea
rank.method.melt <- melt(rank.method.df, id.vars = 'num_sel_fea')
p2 <- ggplot(data = rank.method.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line()
# print(p2)


