# plot the boxplot for auc from efsis_diff_10fold.R

library(reshape2)
library(ggplot2)
library(cowplot)

cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

num.folds <- 10
nums.sel.fea <- seq(4, 50, 2)  # the numbers of selected features
df.merge <- data.frame()  # to record the final merged dataframe
dataset <- 'CNS'

# get the average accuracy (AUC) (averaging on 10 folds)
avr.method.noStab <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(avr.method.noStab) <- methods
var.method.noStab <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(var.method.noStab) <- methods
count <- 0
for (i in nums.sel.fea) {
  count <- count + 1
  file.name <- paste('../data/', dataset, '/num', i, '-auc-stab10-noStab.txt', sep = '')
  df <- read.table(file.name)  # read the file
  avr.method.noStab[, count] <- rowMeans(df)
  var.method.noStab[, count] <- apply(df, 1, var)
  names(avr.method.noStab)[count] <- paste(i, '_features')
  names(var.method.noStab)[count] <- paste(i, '_features')
}
avr.method.noStab.df <- as.data.frame(t(avr.method.noStab))
avr.method.noStab.df$num_sel_fea <- nums.sel.fea

var.method.noStab.df <- as.data.frame(t(var.method.noStab))
var.method.noStab.df$num_sel_fea <- nums.sel.fea

#=====Stability value
avr.method.Stab <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(avr.method.Stab) <- methods
var.method.Stab <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(var.method.Stab) <- methods
count <- 0
for (i in nums.sel.fea) {
  count <- count + 1
  file.name <- paste('../data/', dataset, '/num', i, '-auc-stab10.txt', sep = '')
  df <- read.table(file.name)  # read the file
  avr.method.Stab[, count] <- rowMeans(df)
  var.method.Stab[, count] <- apply(df, 1, var)
  names(avr.method.Stab)[count] <- paste(i, '_features')
  names(var.method.Stab)[count] <- paste(i, '_features')
}
avr.method.Stab.df <- as.data.frame(t(avr.method.Stab))
avr.method.Stab.df$num_sel_fea <- nums.sel.fea

var.method.Stab.df <- as.data.frame(t(var.method.Stab))
var.method.Stab.df$num_sel_fea <- nums.sel.fea

#=====Stability ranking
avr.method.stabRank <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(avr.method.stabRank) <- methods
var.method.stabRank <- data.frame(matrix(NA, nrow = length(methods), ncol = 1))
rownames(var.method.stabRank) <- methods
count <- 0
for (i in nums.sel.fea) {
  count <- count + 1
  file.name <- paste('../data/', dataset, '/num', i, '-auc-stab10-rank.txt', sep = '')
  df <- read.table(file.name)  # read the file
  avr.method.stabRank[, count] <- rowMeans(df)
  var.method.stabRank[, count] <- apply(df, 1, var)
  names(avr.method.stabRank)[count] <- paste(i, '_features')
  names(var.method.stabRank)[count] <- paste(i, '_features')
}
avr.method.stabRank.df <- as.data.frame(t(avr.method.stabRank))
avr.method.stabRank.df$num_sel_fea <- nums.sel.fea

var.method.stabRank.df <- as.data.frame(t(var.method.stabRank))
var.method.stabRank.df$num_sel_fea <- nums.sel.fea

# ========form_sam-geode
form.sam.geode <- data.frame(noStab = avr.method.noStab.df$`form_sam-geode`, stab = avr.method.Stab.df$`form_sam-geode`, stabRank = avr.method.stabRank.df$`form_sam-geode`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
form.sam.geode.melt <- melt(form.sam.geode, id.vars = 'num_sel_fea')
p1 <- ggplot(data = form.sam.geode.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of sam & geode \nwith formula')

# ========form_ref-chs
form.ref.chs <- data.frame(noStab = avr.method.noStab.df$`form_ref-chs`, stab = avr.method.Stab.df$`form_ref-chs`, stabRank = avr.method.stabRank.df$`form_ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
form.ref.chs.melt <- melt(form.ref.chs, id.vars = 'num_sel_fea')
p2 <- ggplot(data = form.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of ref & chs \nwith formula')

# ========form_sam-geode-ref-chs
form.sam.geode.ref.chs <- data.frame(noStab = avr.method.noStab.df$`form_sam-geode-ref-chs`, stab = avr.method.Stab.df$`form_sam-geode-ref-chs`, stabRank = avr.method.stabRank.df$`form_sam-geode-ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
form.sam.geode.ref.chs.melt <- melt(form.sam.geode.ref.chs, id.vars = 'num_sel_fea')
p3 <- ggplot(data = form.sam.geode.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of sam & geode ref & chs \nwith formula')

# ========consensus_sam-geode
consensus.sam.geode <- data.frame(noStab = avr.method.noStab.df$`consensus_sam-geode`, stab = avr.method.Stab.df$`consensus_sam-geode`, stabRank = avr.method.stabRank.df$`consensus_sam-geode`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
consensus.sam.geode.melt <- melt(consensus.sam.geode, id.vars = 'num_sel_fea')
p4 <- ggplot(data = consensus.sam.geode.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of sam & geode \nwith consensus ranking')

# ========consensus_ref-chs
consensus.ref.chs <- data.frame(noStab = avr.method.noStab.df$`consensus_ref-chs`, stab = avr.method.Stab.df$`consensus_ref-chs`, stabRank = avr.method.stabRank.df$`consensus_ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
consensus.ref.chs.melt <- melt(consensus.ref.chs, id.vars = 'num_sel_fea')
p5 <- ggplot(data = consensus.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of ref & chs \nwith consensus ranking')

# ========consensus_sam-geode-ref-chs
consensus.sam.geode.ref.chs <- data.frame(noStab = avr.method.noStab.df$`consensus_sam-geode-ref-chs`, stab = avr.method.Stab.df$`consensus_sam-geode-ref-chs`, stabRank = avr.method.stabRank.df$`consensus_sam-geode-ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(avr.method.Stab.df))
consensus.sam.geode.ref.chs.melt <- melt(consensus.sam.geode.ref.chs, id.vars = 'num_sel_fea')
p6 <- ggplot(data = consensus.sam.geode.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Average AUC across 10 folds') + ggtitle('Ensemble of sam & geode ref & chs \nwith consensus ranking')

fig1 <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

# ================== Variance ================

# ========form_sam-geode
form.sam.geode <- data.frame(noStab = var.method.noStab.df$`form_sam-geode`, stab = var.method.Stab.df$`form_sam-geode`, stabRank = var.method.stabRank.df$`form_sam-geode`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
form.sam.geode.melt <- melt(form.sam.geode, id.vars = 'num_sel_fea')
p11 <- ggplot(data = form.sam.geode.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of sam & geode \nwith formula')

# ========form_ref-chs
form.ref.chs <- data.frame(noStab = var.method.noStab.df$`form_ref-chs`, stab = var.method.Stab.df$`form_ref-chs`, stabRank = var.method.stabRank.df$`form_ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
form.ref.chs.melt <- melt(form.ref.chs, id.vars = 'num_sel_fea')
p12 <- ggplot(data = form.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of ref & chs \nwith formula')

# ========form_sam-geode-ref-chs
form.sam.geode.ref.chs <- data.frame(noStab = var.method.noStab.df$`form_sam-geode-ref-chs`, stab = var.method.Stab.df$`form_sam-geode-ref-chs`, stabRank = var.method.stabRank.df$`form_sam-geode-ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
form.sam.geode.ref.chs.melt <- melt(form.sam.geode.ref.chs, id.vars = 'num_sel_fea')
p13 <- ggplot(data = form.sam.geode.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of sam & geode ref & chs \nwith formula')

# ========consensus_sam-geode
consensus.sam.geode <- data.frame(noStab = var.method.noStab.df$`consensus_sam-geode`, stab = var.method.Stab.df$`consensus_sam-geode`, stabRank = var.method.stabRank.df$`consensus_sam-geode`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.sam.geode.melt <- melt(consensus.sam.geode, id.vars = 'num_sel_fea')
p14 <- ggplot(data = consensus.sam.geode.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of sam & geode \nwith consensus ranking')

# ========consensus_ref-chs
consensus.ref.chs <- data.frame(noStab = var.method.noStab.df$`consensus_ref-chs`, stab = var.method.Stab.df$`consensus_ref-chs`, stabRank = var.method.stabRank.df$`consensus_ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.ref.chs.melt <- melt(consensus.ref.chs, id.vars = 'num_sel_fea')
p15 <- ggplot(data = consensus.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of ref & chs \nwith consensus ranking')

# ========consensus_sam-geode-ref-chs
consensus.sam.geode.ref.chs <- data.frame(noStab = var.method.noStab.df$`consensus_sam-geode-ref-chs`, stab = var.method.Stab.df$`consensus_sam-geode-ref-chs`, stabRank = var.method.stabRank.df$`consensus_sam-geode-ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.sam.geode.ref.chs.melt <- melt(consensus.sam.geode.ref.chs, id.vars = 'num_sel_fea')
p16 <- ggplot(data = consensus.sam.geode.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Variance AUC across 10 folds') + ggtitle('Ensemble of sam & geode ref & chs \nwith consensus ranking')

fig2 <- plot_grid(p11, p12, p13, p14, p15, p16, ncol = 3)

# ==================== Difference of variance between consensus ranking and formula ========

# ========consensus_sam-geode - form_sam-geode
consensus.sam.geode <- data.frame(noStab = var.method.noStab.df$`consensus_sam-geode` - var.method.noStab.df$`form_sam-geode`, stab = var.method.Stab.df$`consensus_sam-geode` - var.method.Stab.df$`form_sam-geode`, stabRank = var.method.stabRank.df$`consensus_sam-geode` - var.method.stabRank.df$`form_sam-geode`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.sam.geode.melt <- melt(consensus.sam.geode, id.vars = 'num_sel_fea')
p21 <- ggplot(data = consensus.sam.geode.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Difference of Variance of AUC across 10 folds') + ggtitle('"Consensus Ranking" - "formula" \nEnsemble of sam & geode')

# ========consensus_ref-chs - form_ref-chs
consensus.ref.chs <- data.frame(noStab = var.method.noStab.df$`consensus_ref-chs` - var.method.noStab.df$`form_ref-chs`, stab = var.method.Stab.df$`consensus_ref-chs` - var.method.Stab.df$`form_ref-chs`, stabRank = var.method.stabRank.df$`consensus_ref-chs` - var.method.stabRank.df$`form_ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.ref.chs.melt <- melt(consensus.ref.chs, id.vars = 'num_sel_fea')
p22 <- ggplot(data = consensus.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Difference of Variance of AUC across 10 folds') + ggtitle('"Consensus Ranking" - "formula" \nEnsemble of ref & chs')

# ========consensus_sam-geode-ref-chs - form_sam-geode-ref-chs
consensus.sam.geode.ref.chs <- data.frame(noStab = var.method.noStab.df$`consensus_sam-geode-ref-chs` - var.method.noStab.df$`form_sam-geode-ref-chs`, stab = var.method.Stab.df$`consensus_sam-geode-ref-chs` - var.method.Stab.df$`form_sam-geode-ref-chs`, stabRank = var.method.stabRank.df$`consensus_sam-geode-ref-chs` - var.method.stabRank.df$`form_sam-geode-ref-chs`, num_sel_fea = nums.sel.fea, row.names = row.names(var.method.Stab.df))
consensus.sam.geode.ref.chs.melt <- melt(consensus.sam.geode.ref.chs, id.vars = 'num_sel_fea')
p23 <- ggplot(data = consensus.sam.geode.ref.chs.melt, aes(num_sel_fea, y=value, linetype = variable, colour = variable)) + geom_line() + xlab('Number of selected features') + ylab('Difference of Variance of AUC across 10 folds') + ggtitle('"Consensus Ranking" - "formula" \nEnsemble of sam & geode ref & chs')

fig3 <- plot_grid(p21, p22, p23, ncol = 3)
