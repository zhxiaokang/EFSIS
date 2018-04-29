# To calculate the p-adj for each method compared with the individual method with the highest average AUC

# Function to calculate the variance of each row
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# The function to calculate the p-adj for each method compared with the individual method with the highest average AUC
pAdj <- function(df, vec){  # df is the dataframe of the auc matrix, vec gives the positions of the individual method (row numbers)
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
  df.t <- data.frame(t(df))
  df.t$control <- unlist(indv.best)
  p.raw <- friedmanPost(df.t, control = ncol(df.t))[-ncol(df.t)]
  p.adj <- p.adjust(p.raw, method = 'BH')
  output.list <- list('p.raw' = p.raw, 'p.adj' = p.adj)
}

num.folds <- 10
nums.sel.fea <- seq(4, 50, 2)  # the numbers of selected features
names.methods <- c('sam', 'geode', 'ref', 'chs', 'form_sam-geode', 'form_ref-chs', 'form_sam-geode-ref-chs', 
                   'consensus_sam-geode', 'consensus_ref_chs', 'consensus_sam-geode-ref-chs')
dataset <- 'CNS'
vec <- c(1:10)  # the indexes of all the individual methods

p.raw.list <- data.frame(row.names = names.methods)  # save the p for all #sel
p.adj.list <- data.frame(row.names = names.methods)  # save the p-adj for all #sel
auc.mean.list <- data.frame(row.names = names.methods)

for (i in nums.sel.fea) {
  file.name <- paste('../data/', dataset, '/num', i, '-auc-stab10-rank.txt', sep = '')
  df <- read.table(file.name)  # read the file
  output.padj <- pAdj(df, vec)
  p.raw <- output.padj$p.raw
  p.adj <- output.padj$p.adj
  p.raw.list[as.character(i)] <- p.raw
  p.adj.list[as.character(i)] <- p.adj
  method.sig <- names.methods[which(p.raw < 0.05)]
  if (length(method.sig) > 0){
    print(paste('num. fold:', i, 'method:', method.sig))
  }
  auc.mean <- colMeans(df)
  auc.mean.list[as.character(i)] <- auc.mean
}

write.table(p.raw.list, paste('../data/', dataset, '/p_value-auc-stab10-rank.txt', sep = ''), quote = FALSE, sep = "\t", row.names = T, col.names = T)
write.table(auc.mean.list, paste('../data/', dataset, '/auc-mean-stab10-rank.txt', sep = ''), quote = FALSE, sep = "\t", row.names = T, col.names = T)
