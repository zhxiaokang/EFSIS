# EFSIS (Ensemble of Feature Selection Integrating Stability
# Ensemble of Chi-Squared and ReliefF, using both formula and consensus ranking to integrate them. Compare with ReliefF, Chi-Squared, RGIFE
# The number of selected features is decided by RGIFE
# Use 10-fold cross-validation to get 10 prediction accuracies

remove(list = ls())

# library(rPython)  # to use RGIFE which is written in Python
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
path.script <- setwd('/export/jonassenfs/xiaokangz/projects/EFSIS/scripts')
path.data <- '../data/ProstateSingh/'  # path to the data
data.file <- 'prostate.singh.arff'
system(paste('cp', '-r', '../rgife/*', path.data))  # the code of RGIFE must be in the same directory as the data
k.folds <- 10  # k-fold cross validation
num.round <- 5  # number of rounds of resampling
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
  
  fea.rank.merge.ref <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.ref) <- fea.name
  
  fea.rank.merge.chs <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.chs) <- fea.name
  
  # the dataframe collecting the top features for all rounds
  
  fea.top.ref <- data.frame()
  fea.top.chs <- data.frame()
  
  ## score the features based on each round of resampling
  
  seed.resample <- 1234
  for (i in c(1:num.round)){
    seed.resample <- seed.resample + 1
    set.seed(seed.resample)
    id.train.resample <- c(sample(which(label.train == levels(label.train)[1]), num.resample.control), sample(which(label.train == levels(label.train)[2]), num.resample.treat))
    x.train.resample <- x.train[, id.train.resample]
    y.train.resample <- y.train[id.train.resample]
    
    # use ReliefF to rank the features
    data.ref <- data.frame(t(x.train.resample), y.train.resample, check.names = F)  # add the param to avoid changing '-' to '.'
    estReliefF <- attrEval('y.train.resample', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
    names(estReliefF) <- fea.name  # This command needs to be added because it's very annoying that 'attrEval' will change the '-' in the names to '.'
    fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
    fea.rank.ref <- data.frame(importance = fea.rank.ref)
    fea.rank.name.ref <- rownames(fea.rank.ref)  # the ranked feature list for this round
    fea.top.ref <- cbind.all(fea.top.ref, fea.rank.name.ref[1:num.sel.fea])  # add the top features for this round to the whole record
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
    fea.top.chs <- cbind.all(fea.top.chs, fea.rank.name.chs[1:num.sel.fea])  # add the top features for this round to the whole record
    colnames(fea.rank.chs) <- 'rank'
    fea.rank.chs['rank'] <- seq(1, num.fea)
    fea.rank.merge.chs <- merge(fea.rank.merge.chs, fea.rank.chs, by = 'row.names')
    row.names(fea.rank.merge.chs) <- fea.rank.merge.chs$Row.names
    fea.rank.merge.chs <- fea.rank.merge.chs[, -1]
  }
  
  # =============== Stability Calculation ===============
  
  stab.ref <- stab(fea.top.ref)
  stab.chs <- stab(fea.top.chs)
  
  # =============== Sub-Final Score Calculation ===============
  ## calculate the sub-final score of each feature
  
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
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.ref.order[, 'final.rank', drop = F], by = 'row.names', all = T)
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.chs.order[, 'final.rank', drop = F], by = 'row.names')
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  colnames(fea.rank.merge.efsis) <- c('final.rank', 'ref.rank', 'chs.rank')
  
  for (line in c(1:num.fea)){
    fea.rank.merge.efsis$final.rank[line] <- (1 - stab.ref^(2)) * fea.rank.merge.ref$final.rank[line] + 
      (1 - stab.chs^(2)) * fea.rank.merge.chs$final.rank[line]
  }
  ## Pick top features as feature subset
  fea.order.efsis <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank), ]
  sel.fea.efsis.form <- row.names(fea.order.efsis[1:num.sel.fea, ])
  
  # use RankAggreg to get the efsis ranking list
  rank.list.ref <- row.names(fea.rank.merge.ref[order(fea.rank.merge.ref[, 'final.rank']), ])
  rank.list.chs <- row.names(fea.rank.merge.chs[order(fea.rank.merge.chs[, 'final.rank']), ])
  rank.matrix <- rbind(rank.list.ref, rank.list.chs)
  rank.list.efsis <- RankAggreg(rank.matrix, num.sel.fea, NULL, method = 'GA', importance = c(stab.ref, stab.chs))
  sel.fea.efsis.consensus <- rank.list.efsis$top.list
  
  # =============== END of Selecting Features using efsis ===============
  output.list <- list('sel.fea.efsis.form' = sel.fea.efsis.form, 'sel.fea.efsis.consensus' = sel.fea.efsis.consensus)
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
set.seed(seed.10fold)
pos.control.train.list <- createFolds(index.control, k.folds, T, T)  # the function gives the position of samples based on the 1st parameter
pos.treat.train.list <- createFolds(index.treat, k.folds, T, T)

# =============== k-folds CV scheme ===============

# data frame to record the predictions (not probabilities)
pred.rgife.svm.folds <- data.frame()
pred.ref.svm.folds <- data.frame()
pred.chs.svm.folds <- data.frame()
pred.efsis.form.folds <- data.frame()
pred.efsis.consensus.folds <- data.frame()

sel.fea.rgife.folds <- data.frame()
sel.fea.ref.folds <- data.frame()
sel.fea.chs.folds <- data.frame()
sel.fea.efsis.form.folds <- data.frame()
sel.fea.efsis.consensus.folds <- data.frame()

for (i in c(1:k.folds)){
  
  # =============== Data Preparation ===============
  
  index.train <- c(index.control[pos.control.train.list[[i]]], index.treat[pos.treat.train.list[[i]]])
  index.train.control <- index.control[pos.control.train.list[[i]]]
  index.train.treat <- index.treat[pos.treat.train.list[[i]]]
  index.test <- setdiff(c(index.control, index.treat), index.train)  # index of all minus index of train gives the index of test
  data.train <- data.raw[index.train, ]  # row -> sample, column -> feature
  data.test <- data.raw[index.test, ]
  
  # save the data
  data.train.file.name <- paste(i, 'fold', '-', 'train', '-', data.file, sep = '')
  data.test.file.name <- paste(i, 'fold', '-', 'test', '-', data.file, sep = '')
  write.arff(data.train, paste(path.data, data.train.file.name, sep = ''))
  write.arff(data.test, paste(path.data, data.test.file.name, sep = ''))
  
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
  
  # =============== Feature Selection using RGIFE ===============
  
  setwd(path.data)
  system(paste('python', ' ', 'rgife.py', ' ', 'configuration.conf', ' ', data.train.file.name, sep = ''))  # python 2.7.13 is needed to execute RGIFE
  rgife.file.name <- paste(i, 'fold', '-', 'selected-rgife.arff', sep = '')  # file name of the selected features by RGIFE
  system(paste('cp', './BestIteration/selected_best_data.arff', rgife.file.name))
  data.sel.fea.rgife <- read.arff(rgife.file.name)
  sel.fea.rgife <- colnames(data.sel.fea.rgife[, -ncol(data.sel.fea.rgife)])  # the selected features by RGIFE
  setwd(path.script)
  num.sel.fea <- length(sel.fea.rgife)
  
  # =============== Feature Selection using ReliefF ===============
  
  print('Start selecting features using ReliefF')
  data.ref <- data.frame(t(x.train), y.train, check.names = F)  # add the param to avoid changing '-' to '.'
  estReliefF <- attrEval('y.train', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
  names(estReliefF) <- fea.name  # It's very annoying that 'attrEval' will change the '-' in the names to '.'
  fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
  sel.fea.ref <- names(fea.rank.ref[1:num.sel.fea])
  
  # =============== Feature Selection using Chi-Squared ===============
  
  print('Start selecting features using Chi-Squared')
  data.chs <- data.frame(t(x.train), y.train, check.names = F)
  weights <- chi.squared(y.train~., data.chs)
  fea.rank.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]
  sel.fea.chs <- row.names(fea.rank.chs[1:num.sel.fea, , drop = F])
  
  # =============== Feature Selection using EFSIS ===============
  
  num.resample.control <- ceiling(length(index.train.control) / num.round)  # number of sampled samples in each round
  num.resample.treat <- ceiling(length(index.train.treat) / num.round)
  output.list.efsis <- efsis()
  sel.fea.efsis.form <- output.list.efsis$sel.fea.efsis.form
  sel.fea.efsis.consensus <- output.list.efsis$sel.fea.efsis.consensus
  
  # save and record the selected features
  
  write.table(sel.fea.rgife, file = paste(path.data, i, 'fold', '-', 'sel-rgife.txt', sep = ''), quote = F, col.names = F, row.names = F)
  write.table(sel.fea.ref, file = paste(path.data, i, 'fold', '-', 'sel-ref.txt', sep = ''), quote = F, col.names = F, row.names = F)
  write.table(sel.fea.chs, file = paste(path.data, i, 'fold', '-', 'sel-chs.txt', sep = ''), quote = F, col.names = F, row.names = F)
  write.table(sel.fea.efsis.form, file = paste(path.data, i, 'fold', '-', 'sel-efsis-form.txt', sep = ''), quote = F, col.names = F, row.names = F)
  write.table(sel.fea.efsis.consensus, file = paste(path.data, i, 'fold', '-', 'sel-efsis-consensus.txt', sep = ''), quote = F, col.names = F, row.names = F)
  
  sel.fea.rgife.folds <- cbind.all(sel.fea.rgife.folds, sel.fea.rgife)
  sel.fea.ref.folds <- cbind.all(sel.fea.ref.folds, sel.fea.ref)
  sel.fea.chs.folds <- cbind.all(sel.fea.chs.folds, sel.fea.chs)
  sel.fea.efsis.form.folds <- cbind.all(sel.fea.efsis.form.folds, sel.fea.efsis.form)
  sel.fea.efsis.consensus.folds <- cbind.all(sel.fea.efsis.consensus.folds, sel.fea.efsis.consensus)
  
  # =============== END of Selecting Features ===============
  
  # =============== Prediction ===============
  
  # Prepare the data sets for classificaion, usually row --> sample, col --> feature
  
  train.rgife <- t(x.train[sel.fea.rgife, ])
  train.ref <- t(x.train[sel.fea.ref, ])
  train.chs <- t(x.train[sel.fea.chs, ])
  train.efsis.form <- t(x.train[sel.fea.efsis.form, ])
  train.efsis.consensus <- t(x.train[sel.fea.efsis.consensus, ])
  
  test.rgife <- t(x.test[sel.fea.rgife, ])
  test.ref <- t(x.test[sel.fea.ref, ])
  test.chs <- t(x.test[sel.fea.chs, ])
  test.efsis.form <- t(x.test[sel.fea.efsis.form, ])
  test.efsis.consensus <- t(x.test[sel.fea.efsis.consensus, ])
  
  # Use SVM to test the performance
  
  ## Prepare the data sets for SVM
  
  train.rgife.svm <- data.frame(cbind(train.rgife, y.train))
  test.rgife.svm <- data.frame(test.rgife)
  
  train.ref.svm <- data.frame(cbind(train.ref, y.train))
  test.ref.svm <- data.frame(test.ref)
  
  train.chs.svm <- data.frame(cbind(train.chs, y.train))
  test.chs.svm <- data.frame(test.chs)
  
  train.efsis.form.svm <- data.frame(cbind(train.efsis.form, y.train))
  test.efsis.form.svm <- data.frame(test.efsis.form)
  
  train.efsis.consensus.svm <- data.frame(cbind(train.efsis.consensus, y.train))
  test.efsis.consensus.svm <- data.frame(test.efsis.consensus)
  
  ## Train SVM classifiers
  
  model.rgife.svm <- svm(formula = factor(y.train)~., data = train.rgife.svm, probability = TRUE)
  model.ref.svm <- svm(formula = factor(y.train)~., data = train.ref.svm, probability = TRUE)
  model.chs.svm <- svm(formula = factor(y.train)~., data = train.chs.svm, probability = TRUE)
  model.efsis.form.svm <- svm(formula = factor(y.train)~., data = train.efsis.form.svm, probability = TRUE)
  model.efsis.consensus.svm <- svm(formula = factor(y.train)~., data = train.efsis.consensus.svm, probability = TRUE)
  
  ## Prediction on the test set
  
  pred.rgife.svm <- attr(predict(model.rgife.svm, test.rgife.svm, probability = TRUE), 'probabilities')[, 2]
  pred.ref.svm <- attr(predict(model.ref.svm, test.ref.svm, probability = TRUE), 'probabilities')[, 2]
  pred.chs.svm <- attr(predict(model.chs.svm, test.chs.svm, probability = TRUE), 'probabilities')[, 2]
  pred.efsis.form.svm <- attr(predict(model.efsis.form.svm, test.efsis.form.svm, probability = TRUE), 'probabilities')[, 2]
  pred.efsis.consensus.svm <- attr(predict(model.efsis.consensus.svm, test.efsis.consensus.svm, probability = TRUE), 'probabilities')[, 2]
  
  ## Test the classification accuracy using the same test set
  
  auc.rgife <- auc(y.test, pred.rgife.svm)
  auc.ref <- auc(y.test, pred.ref.svm)
  auc.chs <- auc(y.test, pred.chs.svm)
  auc.efsis.form <- auc(y.test, pred.efsis.form.svm)
  auc.efsis.consensus <- auc(y.test, pred.efsis.consensus.svm)
  auc.all <- c(auc.rgife, auc.ref, auc.chs, auc.efsis.form, auc.efsis.consensus)
  write.table(auc.all, paste(path.data, i, 'fold', '-', 'auc.txt', sep = ''), quote = F, col.names = F, row.names = F)
  
  # Using ACC
  ## Train SVM classifiers
  
  model.rgife.svm <- svm(formula = factor(y.train)~., data = train.rgife.svm)
  model.ref.svm <- svm(formula = factor(y.train)~., data = train.ref.svm)
  model.chs.svm <- svm(formula = factor(y.train)~., data = train.chs.svm)
  model.efsis.form.svm <- svm(formula = factor(y.train)~., data = train.efsis.form.svm)
  model.efsis.consensus.svm <- svm(formula = factor(y.train)~., data = train.efsis.consensus.svm)
  
  ## Prediction on the test set
  
  pred.rgife.svm <- predict(model.rgife.svm, test.rgife.svm)
  pred.ref.svm <- predict(model.ref.svm, test.ref.svm)
  pred.chs.svm <- predict(model.chs.svm, test.chs.svm)
  pred.efsis.form.svm <- predict(model.efsis.form.svm, test.efsis.form.svm)
  pred.efsis.consensus.svm <- predict(model.efsis.consensus.svm, test.efsis.consensus.svm)
  
  pred.rgife.svm.folds <- cbind.all(pred.rgife.svm.folds, pred.ref.svm)
  pred.ref.svm.folds <- cbind.all(pred.ref.svm.folds, pred.ref.svm)
  pred.chs.svm.folds <- cbind.all(pred.chs.svm.folds, pred.chs.svm)
  pred.efsis.form.folds <- cbind.all(pred.efsis.form.folds, pred.efsis.form.svm)
  pred.efsis.consensus.folds <- cbind.all(pred.efsis.consensus.folds, pred.efsis.consensus.svm)
  
  ## Test the classification accuracy using the same test set
  
}

