# plot the boxplot for auc from efsis_diff_10fold.R

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
data.set <- 'AML'
path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
data.file <- list.files(path = path.data, pattern = '.arff')
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
percent.sel.fea <- c(0.1, 0.2, 0.4, 0.7, 1, 1.5, 2, 3, 4, 5) / 100

# Get the general information about this dataset

num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features

df.merge <- data.frame()  # to record the final merged dataframe
names.methods <- c('sam', 'geode', 'ref', 'infog', 'func', 'efsis', 'css')
stabs <- matrix(nrow = length(nums.sel.fea), ncol = length(names.methods))
j <- 0
for (i in nums.sel.fea) {
  j = j + 1
  file.name <- paste(path.data, 'num', i, '-stab-main-btsp-css.txt', sep = '')
  stab <- read.table(file.name)  # read the file
  stabs[j, ] <- unlist(stab)
}

df.stabs <- as.data.frame.matrix(stabs)
colnames(df.stabs) <- names.methods
df.stabs$num.sel.fea <- nums.sel.fea
stabs.long <- melt(df.stabs, id = 'num.sel.fea')
pic <- ggplot(data = stabs.long, aes(x = num.sel.fea, y = value, linetype = variable, colour = variable)) + geom_line() 

# # get the average performance
# ranks.stabs <- df.stabs
# for (i in c(1:nrow(df.stabs))){
#   ranks.stabs[i, ] <- rank(-df.stabs[i, ])
# }
# ranks.stabs.ave <- colMeans(ranks.stabs)
# 
# # Friedman test
# friedmanPost(ranks.stabs)
# 
# # TurkeyHSD test
# df.stabs <- as.data.frame.matrix(stabs)
# colnames(df.stabs) <- names.methods
# df.stabs <- as.data.frame(t(df.stabs))
# df.stabs$method <- rownames(df.stabs)
# df.stabs.melt <- melt(df.stabs, id.vars = 'method')
# aov.stabs <- aov(df.stabs.melt$value ~ df.stabs.melt$method + df.stabs.melt$variable)
# posthoc.stabs <- TukeyHSD(x=aov.stabs, 'df.stabs.melt$method', conf.level=0.95)
