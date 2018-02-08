# plot the boxplot for auc from efsis_diff_10fold.R

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
nums.sel.fea <- seq(4, 36, 2)  # the numbers of selected features
df.merge <- data.frame()  # to record the final merged dataframe
names.methods <- c('sam', 'geode', 'ref', 'chs', 'form_sam-geode', 'form_ref-chs', 'form_sam-geode-ref-chs', 
                  'consensus_sam-geode', 'consensus_ref_chs', 'consensus_sam-geode-ref-chs')
stabs <- matrix(nrow = length(nums.sel.fea), ncol = length(names.methods))
j <- 0
for (i in nums.sel.fea) {
  j = j + 1
  file.name <- paste('../data/DLBCL/num', i, '-stab.txt', sep = '')
  stab <- read.table(file.name)  # read the file
  stabs[j, ] <- unlist(stab)
}

df.stabs <- as.data.frame.matrix(stabs)
colnames(df.stabs) <- names.methods
df.stabs$num.sel.fea <- nums.sel.fea
stabs.long <- melt(df.stabs, id = 'num.sel.fea')
ggplot(data = stabs.long, aes(x = num.sel.fea, y = value, linetype = variable, colour = variable)) + geom_line() 

# get the average performance
ranks.stabs <- df.stabs
for (i in c(1:nrow(df.stabs))){
  ranks.stabs[i, ] <- rank(-df.stabs[i, ])
}
ranks.stabs.ave <- colMeans(ranks.stabs)

friedmanPost(ranks.stabs)
