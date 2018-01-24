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
nums.sel.fea <- seq(4, 36, 2)  # the numbers of selected features
df.merge <- data.frame()  # to record the final merged dataframe

for (i in nums.sel.fea) {
  file.name <- paste('../data/Leukemia/num', i, '-auc.txt', sep = '')
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

# ggplot(data = df.melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=method))

p <- ggplot(data = df.melt, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=method))
p + facet_wrap( ~ variable, scales="free")

