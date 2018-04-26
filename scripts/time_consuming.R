remove(list = ls())

library(foreign)
library(caret)
library(samr)
library(GeoDE)
library(FSelector)  # to use Chi-Square
library(FSelectorRcpp)  # to use Information Gain
library(CORElearn)  # to use ReliefF                                                                                                                                                                                       
library(RankAggreg)
library(e1071)
library(pROC)
library(reshape)
library(ggplot2)

# ============ Parameters definition ===========
path.script <- setwd('./')  # the path to the script
path.data <- '../data/CNS/'  # path to the data
data.file <- 'cns.arff'
nums.sel.fea <- c(4, 30, 50)
k.folds <- 10  # k-fold cross validation
num.round <- 50  # number of rounds of resampling for EFSIS
seed.10fold <- 12345

# =============== Data Preparation ===============

# Load the data
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature

# Get the general information about this dataset
num.sample <- nrow(data.raw)  # number of samples
num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
pos.label <- num.fea + 1  # the position of the label marked in the dataset
fea.name <- colnames(data.raw[1:num.fea])

# Split into training and test sets
labels <- data.raw[, pos.label]
label.control <- levels(labels)[1]
label.treat <- levels(labels)[2]  # only proper for 2-class classification problem
index.control <- which(labels == label.control)
index.treat <- which(labels == label.treat)
set.seed(seed = seed.10fold)
pos.control.train.list <- createFolds(index.control, k.folds, T, T)  # the function gives the position of samples based on the 1st parameter
set.seed(seed = seed.10fold)
pos.treat.train.list <- createFolds(index.treat, k.folds, T, T)

time.sam <- data.frame('fold' = c(1:k.folds))
time.geode <- data.frame('fold' = c(1:k.folds))
time.ref <- data.frame('fold' = c(1:k.folds))
time.chs <- data.frame('fold' = c(1:k.folds))
time.ig <- data.frame('fold' = c(1:k.folds))

for (num.sel.fea in c(nums.sel.fea)){
  cat('Start the loop of', num.sel.fea, 'selected features:\n')
  
  # =============== k-folds CV scheme ===============
  
  sel.fea.sam.folds <- data.frame()
  sel.fea.geode.folds <- data.frame()
  sel.fea.ref.folds <- data.frame()
  sel.fea.chs.folds <- data.frame()
  sel.fea.efsis.form.sam.geode.folds <- data.frame()
  sel.fea.efsis.form.ref.chs.folds <- data.frame()
  sel.fea.efsis.form.sam.geode.ref.chs.folds <- data.frame()
  sel.fea.efsis.consensus.sam.geode.folds <- data.frame()
  sel.fea.efsis.consensus.ref.chs.folds <- data.frame()
  sel.fea.efsis.consensus.sam.geode.ref.chs.folds <- data.frame()
  
  auc.all.folds <- data.frame(row.names = c('sam', 'geode', 'ref', 'chs', 'form_sam-geode', 'form_ref-chs', 'form_sam-geode-ref-chs', 
                                            'consensus_sam-geode', 'consensus_ref-chs', 'consensus_sam-geode-ref-chs'))
  time.sam.sel <- c()
  time.geode.sel <- c()
  time.ref.sel <- c()
  time.chs.sel <- c()
  time.ig.sel <- c()
  
  for (i in c(1:k.folds)){
    
    # =============== Data Preparation ===============
    cat('Fold', i, ':\n')
    index.train <- c(index.control[pos.control.train.list[[i]]], index.treat[pos.treat.train.list[[i]]])
    index.train.control <- index.control[pos.control.train.list[[i]]]
    index.train.treat <- index.treat[pos.treat.train.list[[i]]]
    index.test <- setdiff(c(index.control, index.treat), index.train)
    data.train <- data.raw[index.train, ]  # row -> sample, column -> feature
    data.test <- data.raw[index.test, ]
    
    #     # save the data
    #     data.train.file.name <- paste(i, 'fold', '-', 'train', '-', data.file, sep = '')
    #     data.test.file.name <- paste(i, 'fold', '-', 'test', '-', data.file, sep = '')
    #     write.arff(data.train, paste(path.data, data.train.file.name, sep = ''))
    #     write.arff(data.test, paste(path.data, data.test.file.name, sep = ''))
    
    fea.name <- colnames(data.train[1:num.fea])
    label.train <- data.train[, pos.label]
    
    # Prepare for the training and test data
    
    x.train.temp <- t(data.train[, 1:num.fea])  # usually row -> features, transform if not
    x.test.temp <- t(data.test[, 1:num.fea])
    
    # Normalization (Avoid using info from test set)
    
    x.mean <- apply(x.train.temp, 1, mean)
    x.sd <- apply(x.train.temp, 1, sd)
    x.train <- (x.train.temp - x.mean) / x.sd  # row --> feature, column --> sample
    x.test <- (x.test.temp - x.mean) / x.sd
    x.train[is.na(x.train)] <- 0
    x.test[is.na(x.test)] <- 0
    
    y.train <- data.train[, pos.label]
    y.test <- data.test[, pos.label]
    
    # =============== Feature Selection using SAM ===============
    
    # print('Start selecting features using SAM')
    start.time <- proc.time()
    label.geode <- factor(y.train)
    levels(label.geode) <- c('1', '2')
    data.sam <- list(x=x.train, y=factor(label.geode), genenames=paste("gene",as.character(fea.name),sep=""), geneid=as.character(fea.name), logged2=TRUE)  
    invisible(capture.output(samr.obj <- samr(data.sam, resp.type="Two class unpaired", nperms=100)))
    invisible(capture.output(delta.table <- samr.compute.delta.table(samr.obj)))
    del <- -3
    siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data.sam, delta.table)
    siggenes.table$genes.lo[,4] <- -as.numeric((siggenes.table$genes.lo[,4]))  # reverse the Score(d) of genes.lo
    siggenes <- rbind(siggenes.table$genes.up, siggenes.table$genes.lo)  # combine the up and low genes
    siggenes <- siggenes[order(siggenes[,4],decreasing = TRUE),]  # sort the genes according to the Score(d)
    fea.rank.sam.pos <- siggenes[, 1]  # all genes sorted according to Score(d) (significance of differentially expressed)
    fea.rank.sam.pos <- as.numeric(fea.rank.sam.pos) - 1  # but this is the position, need to be transferred to ID (row name)
    fea.rank.name.sam <- rep(0, num.fea)
    for (j in c(1:num.fea)){
      fea.rank.name.sam[j] <- fea.name[fea.rank.sam.pos[j]]
    }
    sel.fea.sam <- fea.rank.name.sam[1:num.sel.fea]
    end.time <- proc.time()
    time.taken <- end.time - start.time
#     cat('Time taken by SAM:', '\n')
#     cat(time.taken[1], '\n')
    time.sam.sel <- c(time.sam.sel, time.taken[[1]])
    
    # =============== Feature Selection using GeoDE ===============
    
    # print('Start selecting features using GeoDE')
    start.time <- proc.time()
    gammas <- 1
    data.geode <- data.frame(fea.name, x.train)
    label.geode <- factor(y.train)
    levels(label.geode) <- c('1', '2')  # GeoDE::chdirAnalysis only accepts '1' and '2' as label factors
    chdir.analysis <- chdirAnalysis(data.geode, factor(label.geode), gammas, CalculateSig=FALSE, nnull=10)
    fea.rank.name.geode <- names(chdir.analysis$results[[1]])  # the ranked feature list
    sel.fea.geode <- fea.rank.name.geode[1:num.sel.fea]
    end.time <- proc.time()
    time.taken <- end.time - start.time
#     cat('Time taken by GeoDE:', '\n')
#     cat(time.taken[1], '\n')
    time.geode.sel <- c(time.geode.sel, time.taken[1])
    
    # =============== Feature Selection using ReliefF ===============
    
    # print('Start selecting features using ReliefF')
    start.time <- proc.time()
    data.ref <- data.frame(t(x.train), y.train, check.names = F)  # add the param to avoid changing '-' to '.'
    estReliefF <- attrEval('y.train', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
    names(estReliefF) <- fea.name  # It's very annoying that 'attrEval' will change the '-' in the names to '.'
    fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
    sel.fea.ref <- names(fea.rank.ref[1:num.sel.fea])
    end.time <- proc.time()
    time.taken <- end.time - start.time
#     cat('Time taken by ReliefF:', '\n')
#     cat(time.taken[1], '\n')
    time.ref.sel <- c(time.ref.sel, time.taken[1])

    # =============== Feature Selection using Chi-Squared ===============
    
    # print('Start selecting features using Chi-Squared')
    start.time <- proc.time()
    data.chs <- data.frame(t(x.train), y.train, check.names = F)
    weights <- chi.squared(y.train~., data.chs)
    fea.rank.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
    sel.fea.chs <- row.names(fea.rank.chs[1:num.sel.fea, , drop = F])
    end.time <- proc.time()
    time.taken <- end.time - start.time
#     cat('Time taken by Chi-Squared:', '\n')
#     cat(time.taken[1], '\n')
    time.chs.sel <- c(time.chs.sel, time.taken[1])
    
    # =============== Feature Selection using Information Gain ===============
    
    # print('Start selecting features using Information Gain')
    start.time <- proc.time()
    data.ig <- data.frame(t(x.train), y.train, check.names = F)
    weights <- information_gain(y.train~., data.ig)
    fea.rank.ig <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
    sel.fea.ig <- row.names(fea.rank.ig[1:num.sel.fea, , drop = F])
    end.time <- proc.time()
    time.taken <- end.time - start.time
    #     cat('Time taken by Chi-Squared:', '\n')
    #     cat(time.taken[1], '\n')
    time.ig.sel <- c(time.ig.sel, time.taken[1])
  }
  time.sam[paste(num.sel.fea, '_features', ssep = '')] <- time.sam.sel
  time.geode[paste(num.sel.fea, '_features', ssep = '')] <- time.geode.sel
  time.ref[paste(num.sel.fea, '_features', ssep = '')] <- time.ref.sel
  time.chs[paste(num.sel.fea, '_features', ssep = '')] <- time.chs.sel
  time.ig[paste(num.sel.fea, '_features', ssep = '')] <- time.ig.sel
}
    
#     # =============== Feature Selection using Chi-Squared ===============
#     
#     # print('Start selecting features using Chi-Squared')
#     start.time <- Sys.time()
#     data.chs <- data.frame(t(x.train), y.train, check.names = F)
#     weights <- chi.squared(y.train~., data.chs)
#     fea.rank.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
#     sel.fea.chs <- row.names(fea.rank.chs[1:num.sel.fea, , drop = F])
#     end.time <- Sys.time()
#     time.taken <- end.time - start.time
#     print('Time taken by Chi-Squared:')
#     time.taken
#     print(paste('Time taken by Chi-Squared:', time.taken))

