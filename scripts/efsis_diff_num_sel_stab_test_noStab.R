# EFSIS (Ensemble of Feature Selection Integrating Stability
# Ensemble of SAM, GeoDE, Chi-Squared and ReliefF, using both formula and consensus ranking to integrate them. Compare with SAM, GeoDE, ReliefF, Chi-Squared
# Use 10-fold cross-validation to get 10 prediction accuracies
# Inherited from 'efsis_diff_num_sel_stab_test.R', but not using stability when integrating functions
remove(list = ls())

library(foreign)
library(caret)
library(samr)
library(GeoDE)
library(FSelector)  # to use Chi-Square
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
nums.sel.fea <- seq(4, 50, 2)
k.folds <- 10  # k-fold cross validation
num.round <- 5  # number of rounds of resampling for EFSIS
seed.10fold <- 12345

# ============ function definition ===========

# define the function of calculating stability
stab <- function(dataFrame){
  num.round <- num.round  # number of resampling rounds
  union <- Reduce(union, dataFrame)  # F in the formula: the list of all features, which have been selected in at least one of n sampling steps
  whole <- unlist(list(dataFrame))  # all the features in List (with replicate)
  sum <- 0
  for(fea in union){
    freq <- length(which(whole == fea))
    sum <- sum + freq
  }
  stab.na <- sum / (num.round * length(union))
  return(stab.na)
}

# define the function to add column to an empty dataframe / matrix
cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# define the function of EFSIS
efsis <- function(){
  
  # Feature selection using resampling
  
  # the dataframe collecting all the scores for each feature
  
  fea.rank.merge.sam <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.sam) <- fea.name
  
  fea.rank.merge.geode <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.geode) <- fea.name
  
  fea.rank.merge.ref <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.ref) <- fea.name
  
  fea.rank.merge.chs <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.chs) <- fea.name
  
  # the dataframe collecting the top features for all rounds
  
  fea.top.sam <- data.frame()
  fea.top.geode <- data.frame()
  fea.top.ref <- data.frame()
  fea.top.chs <- data.frame()
  
  ## score the features based on each round of resampling
  
  seed.resample <- 1234
  for (i in c(1:num.round)){
    seed.resample <- seed.resample + 1
    set.seed(seed = seed.resample)
    id.train.resample <- c(sample(which(label.train == levels(label.train)[1]), num.resample.control), sample(which(label.train == levels(label.train)[2]), num.resample.treat))
    x.train.resample <- x.train[, id.train.resample]
    y.train.resample <- y.train[id.train.resample]
    
    # use SAM to rank the features
    label.geode <- factor(y.train.resample)
    levels(label.geode) <- c('1', '2')
    data.sam <- list(x=x.train.resample, y=factor(label.geode), genenames=paste("gene",as.character(fea.name),sep=""), geneid=as.character(fea.name), logged2=TRUE)  
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
    fea.top.sam <- cbind.all(fea.top.sam, fea.rank.name.sam[1:(num.sel.fea*10)])  # add the top features from this round to the whole record
    fea.rank.sam <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.sam) <- fea.rank.name.sam
    fea.rank.merge.sam <- merge(fea.rank.merge.sam, fea.rank.sam, by = 'row.names')
    row.names(fea.rank.merge.sam) <- fea.rank.merge.sam$Row.names
    fea.rank.merge.sam <- fea.rank.merge.sam[, -1]
    
    
    # use GeoDE to rank the features
    gammas <- 1
    data.geode <- data.frame(fea.name, x.train.resample)
    label.geode <- factor(y.train.resample)
    levels(label.geode) <- c('1', '2')  # GeoDE::chdirAnalysis only accepts '1' and '2' as label factors
    chdir.analysis <- chdirAnalysis(data.geode, factor(label.geode), gammas, CalculateSig=FALSE, nnull=10)
    fea.rank.name.geode <- names(chdir.analysis$results[[1]])  # the ranked feature list for this round
    fea.top.geode <- cbind.all(fea.top.geode, fea.rank.name.geode[1:(num.sel.fea*10)])  # all the top features for this round to the whole record
    fea.rank.geode <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.geode) <- fea.rank.name.geode
    fea.rank.merge.geode <- merge(fea.rank.merge.geode, fea.rank.geode, by = 'row.names')
    row.names(fea.rank.merge.geode) <- fea.rank.merge.geode$Row.names
    fea.rank.merge.geode <- fea.rank.merge.geode[, -1]
    
    # use ReliefF to rank the features
    data.ref <- data.frame(t(x.train.resample), y.train.resample, check.names = F)  # add the param to avoid changing '-' to '.'
    estReliefF <- attrEval('y.train.resample', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
    names(estReliefF) <- fea.name  # This command needs to be added because it's very annoying that 'attrEval' will change the '-' in the names to '.'
    fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
    fea.rank.ref <- data.frame(importance = fea.rank.ref)
    fea.rank.name.ref <- rownames(fea.rank.ref)  # the ranked feature list for this round
    fea.top.ref <- cbind.all(fea.top.ref, fea.rank.name.ref[1:(num.sel.fea*10)])  # add the top features for this round to the whole record
    colnames(fea.rank.ref) <- 'rank'
    fea.rank.ref['rank'] <- seq(1, num.fea)
    fea.rank.merge.ref <- merge(fea.rank.merge.ref, fea.rank.ref, by = 'row.names')
    row.names(fea.rank.merge.ref) <- fea.rank.merge.ref$Row.names
    fea.rank.merge.ref <- fea.rank.merge.ref[, -1]
    
    # use Chi-Squared to rank the features
    data.chs <- data.frame(t(x.train.resample), y.train.resample, check.names = F)
    weights <- chi.squared(y.train.resample~., data.chs)
    fea.rank.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
    fea.rank.name.chs <- rownames(fea.rank.chs)  # the ranked feature list for this round
    fea.top.chs <- cbind.all(fea.top.chs, fea.rank.name.chs[1:(num.sel.fea*10)])  # add the top features for this round to the whole record
    colnames(fea.rank.chs) <- 'rank'
    fea.rank.chs['rank'] <- seq(1, num.fea)
    fea.rank.merge.chs <- merge(fea.rank.merge.chs, fea.rank.chs, by = 'row.names')
    row.names(fea.rank.merge.chs) <- fea.rank.merge.chs$Row.names
    fea.rank.merge.chs <- fea.rank.merge.chs[, -1]
  }
  
  # =============== Stability Calculation ===============
  
  stab.sam <- stab(fea.top.sam)
  stab.geode <- stab(fea.top.geode)
  stab.ref <- stab(fea.top.ref)
  stab.chs <- stab(fea.top.chs)
  
  # =============== Sub-Final Score Calculation ===============
  ## calculate the sub-final score of each feature
  
  ### sub-final score for SAM
  for (line in c(1:num.fea)){
    fea.rank.merge.sam$final.rank[line] <- prod(fea.rank.merge.sam[line, ][-1])
  }
  fea.rank.merge.sam.order <- fea.rank.merge.sam[order(fea.rank.merge.sam$final.rank), , drop = F]
  fea.rank.merge.sam.order$final.rank <- seq(1, num.fea)
  
  ### sub-final score for GeoDE
  for (line in c(1:num.fea)){
    fea.rank.merge.geode$final.rank[line] <- prod(fea.rank.merge.geode[line, ][-1])
  }
  fea.rank.merge.geode.order <- fea.rank.merge.geode[order(fea.rank.merge.geode$final.rank), , drop = F]
  fea.rank.merge.geode.order$final.rank <- seq(1, num.fea)
  
  ### sub-final score for ref
  for (line in c(1:num.fea)){
    fea.rank.merge.ref$final.rank[line] <- prod(fea.rank.merge.ref[line, ][-1])
  }
  fea.rank.merge.ref.order <- fea.rank.merge.ref[order(fea.rank.merge.ref$final.rank), , drop = F]
  fea.rank.merge.ref.order$final.rank <- seq(1, num.fea)
  
  ### sub-final score for chs
  for (line in c(1:num.fea)){
    fea.rank.merge.chs$final.rank[line] <- prod(fea.rank.merge.chs[line, ][-1])
  }
  fea.rank.merge.chs.order <- fea.rank.merge.chs[order(fea.rank.merge.chs$final.rank), , drop = F]
  fea.rank.merge.chs.order$final.rank <- seq(1, num.fea)
  
  # =============== Final Score Calculation ===============
  
  # use formula to get the EFSIS ranking list
  fea.rank.merge.efsis <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.efsis) <- fea.name
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.sam.order[, 'final.rank', drop = F], by = 'row.names', all = T)
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.geode.order[, 'final.rank', drop = F], by = 'row.names')
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.ref.order[, 'final.rank', drop = F], by = 'row.names', all = T)
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.chs.order[, 'final.rank', drop = F], by = 'row.names')
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  colnames(fea.rank.merge.efsis) <- c('final.rank', 'sam.rank', 'geode.rank', 'ref.rank', 'chs.rank')
  
  for (line in c(1:num.fea)){
    # integrage SAM, GeoDE
    fea.rank.merge.efsis$final.rank.sam.geode[line] <- (fea.rank.merge.sam$final.rank[line]) * 
      (fea.rank.merge.geode$final.rank[line])
    # integrage Ref, Chs
    fea.rank.merge.efsis$final.rank.ref.chs[line] <- (fea.rank.merge.ref$final.rank[line]) * 
      (fea.rank.merge.chs$final.rank[line])
    # integrage SAM, GeoDE, Ref, Chs
    fea.rank.merge.efsis$final.rank.sam.geode.ref.chs[line] <- (fea.rank.merge.sam$final.rank[line]) * 
      (fea.rank.merge.geode$final.rank[line]) * (fea.rank.merge.ref$final.rank[line]) * 
      (fea.rank.merge.chs$final.rank[line])
  }
  ## Pick top features as feature subset
  # integrage SAM, GeoDE
  fea.order.efsis.sam.geode <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank.sam.geode), ]
  sel.fea.efsis.form.sam.geode <- row.names(fea.order.efsis.sam.geode[1:num.sel.fea, ])
  # integrage Ref, Chs
  fea.order.efsis.ref.chs <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank.ref.chs), ]
  sel.fea.efsis.form.ref.chs <- row.names(fea.order.efsis.ref.chs[1:num.sel.fea, ])
  # integrage SAM, GeoDE, Ref, Chs
  fea.order.efsis.sam.geode.ref.chs <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank.sam.geode.ref.chs), ]
  sel.fea.efsis.form.sam.geode.ref.chs <- row.names(fea.order.efsis.sam.geode.ref.chs[1:num.sel.fea, ])
  
  # use RankAggreg to get the efsis ranking list
  rank.list.sam <- row.names(fea.rank.merge.sam[order(fea.rank.merge.sam[, 'final.rank']), ])
  rank.list.geode <- row.names(fea.rank.merge.geode[order(fea.rank.merge.geode[, 'final.rank']), ])
  rank.list.ref <- row.names(fea.rank.merge.ref[order(fea.rank.merge.ref[, 'final.rank']), ])
  rank.list.chs <- row.names(fea.rank.merge.chs[order(fea.rank.merge.chs[, 'final.rank']), ])
  # integrage SAM, GeoDE
  rank.matrix.sam.geode <- rbind(rank.list.sam, rank.list.geode)
  invisible(capture.output(rank.list.efsis.sam.geode <- RankAggreg(rank.matrix.sam.geode, num.sel.fea, NULL, method = 'GA')))
  sel.fea.efsis.consensus.sam.geode <- rank.list.efsis.sam.geode$top.list
  # integrage Ref, Chs
  rank.matrix.ref.chs <- rbind(rank.list.ref, rank.list.chs)
  invisible(capture.output(rank.list.efsis.ref.chs <- RankAggreg(rank.matrix.ref.chs, num.sel.fea, NULL, method = 'GA')))
  sel.fea.efsis.consensus.ref.chs <- rank.list.efsis.ref.chs$top.list
  # integrage SAM, GeoDE, Ref, Chs
  rank.matrix.sam.geode.ref.chs <- rbind(rank.list.sam, rank.list.geode, rank.list.ref, rank.list.chs)
  invisible(capture.output(rank.list.efsis.sam.geode.ref.chs <- RankAggreg(rank.matrix.sam.geode.ref.chs, num.sel.fea, NULL, method = 'GA')))
  sel.fea.efsis.consensus.sam.geode.ref.chs <- rank.list.efsis.sam.geode.ref.chs$top.list
  
  # =============== END of Selecting Features using efsis ===============
  output.list <- list('sel.fea.efsis.form.sam.geode' = sel.fea.efsis.form.sam.geode, 'sel.fea.efsis.form.ref.chs' = sel.fea.efsis.form.ref.chs, 'sel.fea.efsis.form.sam.geode.ref.chs' = sel.fea.efsis.form.sam.geode.ref.chs, 
                      'sel.fea.efsis.consensus.sam.geode' = sel.fea.efsis.consensus.sam.geode, 'sel.fea.efsis.consensus.ref.chs' = sel.fea.efsis.consensus.ref.chs, 'sel.fea.efsis.consensus.sam.geode.ref.chs' = sel.fea.efsis.consensus.sam.geode.ref.chs,
                      'stab.sam' = stab.sam, 'stab.geode' = stab.geode, 'stab.ref' = stab.ref, 'stab.chs' = stab.chs)
}


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

# Loop of different numbers of of selected features
for (num.sel.fea in c(nums.sel.fea)){
  print(paste('Start the loop of', num.sel.fea, 'selected features'))
  
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
  
  for (i in c(1:k.folds)){
    
    # =============== Data Preparation ===============
    
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
    
    # =============== Feature Selection using GeoDE ===============
    
    # print('Start selecting features using GeoDE')
    gammas <- 1
    data.geode <- data.frame(fea.name, x.train)
    label.geode <- factor(y.train)
    levels(label.geode) <- c('1', '2')  # GeoDE::chdirAnalysis only accepts '1' and '2' as label factors
    chdir.analysis <- chdirAnalysis(data.geode, factor(label.geode), gammas, CalculateSig=FALSE, nnull=10)
    fea.rank.name.geode <- names(chdir.analysis$results[[1]])  # the ranked feature list
    sel.fea.geode <- fea.rank.name.geode[1:num.sel.fea]
    
    # =============== Feature Selection using ReliefF ===============
    
    # print('Start selecting features using ReliefF')
    data.ref <- data.frame(t(x.train), y.train, check.names = F)  # add the param to avoid changing '-' to '.'
    estReliefF <- attrEval('y.train', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
    names(estReliefF) <- fea.name  # It's very annoying that 'attrEval' will change the '-' in the names to '.'
    fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
    sel.fea.ref <- names(fea.rank.ref[1:num.sel.fea])
    
    # =============== Feature Selection using Chi-Squared ===============
    
    # print('Start selecting features using Chi-Squared')
    data.chs <- data.frame(t(x.train), y.train, check.names = F)
    weights <- chi.squared(y.train~., data.chs)
    fea.rank.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
    sel.fea.chs <- row.names(fea.rank.chs[1:num.sel.fea, , drop = F])
    
    # =============== Feature Selection using EFSIS ===============
    
    # print('Start selecting features using EFSIS')
    num.resample.control <- ceiling(length(index.train.control) / num.round)  # number of sampled samples in each round
    num.resample.treat <- ceiling(length(index.train.treat) / num.round)
    output.list.efsis <- efsis()
    sel.fea.efsis.form.sam.geode <- output.list.efsis$sel.fea.efsis.form.sam.geode
    sel.fea.efsis.form.ref.chs <- output.list.efsis$sel.fea.efsis.form.ref.chs
    sel.fea.efsis.form.sam.geode.ref.chs <- output.list.efsis$sel.fea.efsis.form.sam.geode.ref.chs
    sel.fea.efsis.consensus.sam.geode <- output.list.efsis$sel.fea.efsis.consensus.sam.geode
    sel.fea.efsis.consensus.ref.chs <- output.list.efsis$sel.fea.efsis.consensus.ref.chs
    sel.fea.efsis.consensus.sam.geode.ref.chs <- output.list.efsis$sel.fea.efsis.consensus.sam.geode.ref.chs
    stab.sam <- output.list.efsis$stab.sam
    stab.geode <- output.list.efsis$stab.geode
    stab.ref <- output.list.efsis$stab.ref
    stab.chs <- output.list.efsis$stab.chs
    
    # save the selected features
    
    sel.fea.sam.folds <- cbind.all(sel.fea.sam.folds, sel.fea.sam)
    sel.fea.geode.folds <- cbind.all(sel.fea.geode.folds, sel.fea.geode)
    sel.fea.ref.folds <- cbind.all(sel.fea.ref.folds, sel.fea.ref)
    sel.fea.chs.folds <- cbind.all(sel.fea.chs.folds, sel.fea.chs)
    sel.fea.efsis.form.sam.geode.folds <- cbind.all(sel.fea.efsis.form.sam.geode.folds, sel.fea.efsis.form.sam.geode)
    sel.fea.efsis.form.ref.chs.folds <- cbind.all(sel.fea.efsis.form.ref.chs.folds, sel.fea.efsis.form.ref.chs)
    sel.fea.efsis.form.sam.geode.ref.chs.folds <- cbind.all(sel.fea.efsis.form.sam.geode.ref.chs.folds, sel.fea.efsis.form.sam.geode.ref.chs)
    sel.fea.efsis.consensus.sam.geode.folds <- cbind.all(sel.fea.efsis.consensus.sam.geode.folds, sel.fea.efsis.consensus.sam.geode)
    sel.fea.efsis.consensus.ref.chs.folds <- cbind.all(sel.fea.efsis.consensus.ref.chs.folds, sel.fea.efsis.consensus.ref.chs)
    sel.fea.efsis.consensus.sam.geode.ref.chs.folds <- cbind.all(sel.fea.efsis.consensus.sam.geode.ref.chs.folds, sel.fea.efsis.consensus.sam.geode.ref.chs)
    
    # print('Finished selecting features')
    
    # =============== END of Selecting Features ===============
    
    # =============== Prediction ===============
    
    # Prepare the data sets for classificaion, usually row --> sample, col --> feature
    
    train.sam <- t(x.train[sel.fea.sam, ])
    train.geode <- t(x.train[sel.fea.geode, ])
    train.ref <- t(x.train[sel.fea.ref, ])
    train.chs <- t(x.train[sel.fea.chs, ])
    train.efsis.form.sam.geode <- t(x.train[sel.fea.efsis.form.sam.geode, ])
    train.efsis.form.ref.chs <- t(x.train[sel.fea.efsis.form.ref.chs, ])
    train.efsis.form.sam.geode.ref.chs <- t(x.train[sel.fea.efsis.form.sam.geode.ref.chs, ])
    train.efsis.consensus.sam.geode <- t(x.train[sel.fea.efsis.consensus.sam.geode, ])
    train.efsis.consensus.ref.chs <- t(x.train[sel.fea.efsis.consensus.ref.chs, ])
    train.efsis.consensus.sam.geode.ref.chs <- t(x.train[sel.fea.efsis.consensus.sam.geode.ref.chs, ])
    
    test.sam <- t(x.test[sel.fea.sam, ])
    test.geode <- t(x.test[sel.fea.geode, ])
    test.ref <- t(x.test[sel.fea.ref, ])
    test.chs <- t(x.test[sel.fea.chs, ])
    test.efsis.form.sam.geode <- t(x.test[sel.fea.efsis.form.sam.geode, ])
    test.efsis.form.ref.chs <- t(x.test[sel.fea.efsis.form.ref.chs, ])
    test.efsis.form.sam.geode.ref.chs <- t(x.test[sel.fea.efsis.form.sam.geode.ref.chs, ])
    test.efsis.consensus.sam.geode <- t(x.test[sel.fea.efsis.consensus.sam.geode, ])
    test.efsis.consensus.ref.chs <- t(x.test[sel.fea.efsis.consensus.ref.chs, ])
    test.efsis.consensus.sam.geode.ref.chs <- t(x.test[sel.fea.efsis.consensus.sam.geode.ref.chs, ])
    
    # Use SVM to test the performance
    
    ## Prepare the data sets for SVM
    
    train.sam.svm <- data.frame(cbind(train.sam, y.train))
    test.sam.svm <- data.frame(test.sam)
    
    train.geode.svm <- data.frame(cbind(train.geode, y.train))
    test.geode.svm <- data.frame(test.geode)
    
    train.ref.svm <- data.frame(cbind(train.ref, y.train))
    test.ref.svm <- data.frame(test.ref)
    
    train.chs.svm <- data.frame(cbind(train.chs, y.train))
    test.chs.svm <- data.frame(test.chs)
    
    train.efsis.form.sam.geode.svm <- data.frame(cbind(train.efsis.form.sam.geode, y.train))
    test.efsis.form.sam.geode.svm <- data.frame(test.efsis.form.sam.geode)
    
    train.efsis.form.ref.chs.svm <- data.frame(cbind(train.efsis.form.ref.chs, y.train))
    test.efsis.form.ref.chs.svm <- data.frame(test.efsis.form.ref.chs)
    
    train.efsis.form.sam.geode.ref.chs.svm <- data.frame(cbind(train.efsis.form.sam.geode.ref.chs, y.train))
    test.efsis.form.sam.geode.ref.chs.svm <- data.frame(test.efsis.form.sam.geode.ref.chs)
    
    train.efsis.consensus.sam.geode.svm <- data.frame(cbind(train.efsis.consensus.sam.geode, y.train))
    test.efsis.consensus.sam.geode.svm <- data.frame(test.efsis.consensus.sam.geode)
    
    train.efsis.consensus.ref.chs.svm <- data.frame(cbind(train.efsis.consensus.ref.chs, y.train))
    test.efsis.consensus.ref.chs.svm <- data.frame(test.efsis.consensus.ref.chs)
    
    train.efsis.consensus.sam.geode.ref.chs.svm <- data.frame(cbind(train.efsis.consensus.sam.geode.ref.chs, y.train))
    test.efsis.consensus.sam.geode.ref.chs.svm <- data.frame(test.efsis.consensus.sam.geode.ref.chs)
    
    ## Train SVM classifiers
    
    model.sam.svm <- svm(formula = factor(y.train)~., data = train.sam.svm, probability = TRUE)
    model.geode.svm <- svm(formula = factor(y.train)~., data = train.geode.svm, probability = TRUE)
    model.ref.svm <- svm(formula = factor(y.train)~., data = train.ref.svm, probability = TRUE)
    model.chs.svm <- svm(formula = factor(y.train)~., data = train.chs.svm, probability = TRUE)
    model.efsis.form.sam.geode.svm <- svm(formula = factor(y.train)~., data = train.efsis.form.sam.geode.svm, probability = TRUE)
    model.efsis.form.ref.chs.svm <- svm(formula = factor(y.train)~., data = train.efsis.form.ref.chs.svm, probability = TRUE)
    model.efsis.form.sam.geode.ref.chs.svm <- svm(formula = factor(y.train)~., data = train.efsis.form.sam.geode.ref.chs.svm, probability = TRUE)
    model.efsis.consensus.sam.geode.svm <- svm(formula = factor(y.train)~., data = train.efsis.consensus.sam.geode.svm, probability = TRUE)
    model.efsis.consensus.ref.chs.svm <- svm(formula = factor(y.train)~., data = train.efsis.consensus.ref.chs.svm, probability = TRUE)
    model.efsis.consensus.sam.geode.ref.chs.svm <- svm(formula = factor(y.train)~., data = train.efsis.consensus.sam.geode.ref.chs.svm, probability = TRUE)
    
    ## Prediction on the test set
    
    pred.sam.svm <- attr(predict(model.sam.svm, test.sam.svm, probability = TRUE), 'probabilities')[, 2]
    pred.geode.svm <- attr(predict(model.geode.svm, test.geode.svm, probability = TRUE), 'probabilities')[, 2]
    pred.ref.svm <- attr(predict(model.ref.svm, test.ref.svm, probability = TRUE), 'probabilities')[, 2]
    pred.chs.svm <- attr(predict(model.chs.svm, test.chs.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.form.sam.geode.svm <- attr(predict(model.efsis.form.sam.geode.svm, test.efsis.form.sam.geode.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.form.ref.chs.svm <- attr(predict(model.efsis.form.ref.chs.svm, test.efsis.form.ref.chs.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.form.sam.geode.ref.chs.svm <- attr(predict(model.efsis.form.sam.geode.ref.chs.svm, test.efsis.form.sam.geode.ref.chs.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.consensus.sam.geode.svm <- attr(predict(model.efsis.consensus.sam.geode.svm, test.efsis.consensus.sam.geode.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.consensus.ref.chs.svm <- attr(predict(model.efsis.consensus.ref.chs.svm, test.efsis.consensus.ref.chs.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.consensus.sam.geode.ref.chs.svm <- attr(predict(model.efsis.consensus.sam.geode.ref.chs.svm, test.efsis.consensus.sam.geode.ref.chs.svm, probability = TRUE), 'probabilities')[, 2]
    
    ## Test the classification accuracy using the same test set
    
    auc.sam <- auc(y.test, pred.sam.svm)
    auc.geode <- auc(y.test, pred.geode.svm)
    auc.ref <- auc(y.test, pred.ref.svm)
    auc.chs <- auc(y.test, pred.chs.svm)
    auc.efsis.form.sam.geode <- auc(y.test, pred.efsis.form.sam.geode.svm)
    auc.efsis.form.ref.chs <- auc(y.test, pred.efsis.form.ref.chs.svm)
    auc.efsis.form.sam.geode.ref.chs <- auc(y.test, pred.efsis.form.sam.geode.ref.chs.svm)
    auc.efsis.consensus.sam.geode <- auc(y.test, pred.efsis.consensus.sam.geode.svm)
    auc.efsis.consensus.ref.chs <- auc(y.test, pred.efsis.consensus.ref.chs.svm)
    auc.efsis.consensus.sam.geode.ref.chs <- auc(y.test, pred.efsis.consensus.sam.geode.ref.chs.svm)
    auc.all <- c(auc.sam, auc.geode, auc.ref, auc.chs, auc.efsis.form.sam.geode, auc.efsis.form.ref.chs, auc.efsis.form.sam.geode.ref.chs, 
                 auc.efsis.consensus.sam.geode, auc.efsis.consensus.ref.chs, auc.efsis.consensus.sam.geode.ref.chs)
    auc.all.folds <- cbind.all(auc.all.folds, auc.all)
  }
  
  # save the auc for this #-sel-fea to file
  write.table(auc.all.folds, paste(path.data, 'num', num.sel.fea, '-', 'auc-stab10-noStab.txt', sep = ''), quote = F, col.names = T, row.names = T)
  
  # calculate the stability after 10 folds
  stab.sam <- stab(sel.fea.sam.folds)
  stab.geode <- stab(sel.fea.geode.folds)
  stab.ref <- stab(sel.fea.ref.folds)
  stab.chs <- stab(sel.fea.chs.folds)
  stab.efsis.form.sam.geode <- stab(sel.fea.efsis.form.sam.geode.folds)
  stab.efsis.form.ref.chs <- stab(sel.fea.efsis.form.ref.chs.folds)
  stab.efsis.form.sam.geode.ref.chs <- stab(sel.fea.efsis.form.sam.geode.ref.chs.folds)
  stab.efsis.consensus.sam.geode <- stab(sel.fea.efsis.consensus.sam.geode.folds)
  stab.efsis.consensus.ref.chs <- stab(sel.fea.efsis.consensus.ref.chs.folds)
  stab.efsis.consensus.sam.geode.ref.chs <- stab(sel.fea.efsis.consensus.sam.geode.ref.chs.folds)
  
  stab.all <- c(stab.sam, stab.geode, stab.ref, stab.chs, stab.efsis.form.sam.geode, stab.efsis.form.ref.chs, stab.efsis.form.sam.geode.ref.chs,
                stab.efsis.consensus.sam.geode, stab.efsis.consensus.ref.chs, stab.efsis.consensus.sam.geode.ref.chs)
  # save the stability for this #-sel-fea to file
  write.table(stab.all, paste(path.data, 'num', num.sel.fea, '-', 'stab-stab10-noStab.txt', sep = ''), quote = F, col.names = F, row.names = F)
}


