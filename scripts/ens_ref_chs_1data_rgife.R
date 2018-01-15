# Ensemble of FS methods & sample resampling
# Ensemble of Chi-Squared and ReliefF. Compare with ReliefF, Chi-Squared
# The data is not divided into training and test sets yet, so needs to be done in the code
# The number of selected features is decided by RGIFE

remove(list = ls())

# library(rPython)  # to use RGIFE which is written in Python
library(foreign)
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

path <- '../data/CNS'  # path to the data
sep.file <- ' '  # the seperator used in the files loaded
num.sel.fea <- 7  # number of selected features
num.round <- 5  # number of rounds of resampling
label.control <- '1'
label.treat <- '0'
num.resample.control <- 4  # number of sampled samples in each round
num.resample.treat <- 6
num.train.control <- 16  # number of control samples in training set
num.train.treat <- 29
num.fea <- 7129  # number of features
pos.label <- num.fea + 1  # the position of the label marked in the dataset (num.fea + 1 if after all attributes)

wt <- 10  # the weight to be assigned to stability

# ============ function definition ===========

# define the function of calculating stability
stab <- function(dataFrame, num.round, num.sel.fea){
  num.round <- num.round  # number of resampling rounds
  num.select.fea <- num.sel.fea
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

# =============== Data Preparation ===============

# Load the data
data.raw <- read.arff(paste(path, 'cns_60x7129.arff', sep = '/'))  # row -> sample, column -> feature
data.sel.fea.rgife <- read.arff(paste(path, 'selected_best_data_cns_train.arff', sep = '/'))
sel.fea.rgife <- colnames(data.sel.fea.rgife[, -ncol(data.sel.fea.rgife)])

# Split into training and test sets
labels <- data.raw[, pos.label]
index.control <- which(labels == label.control)
index.treat <- which(labels == label.treat)
set.seed(1234)
index.train <- c(sample(index.control, num.train.control), sample(index.treat, num.train.treat))
index.test <- c(index.control, index.treat)[-index.train]

data.train <- data.raw[index.train, ]  # row -> sample, column -> feature
data.test <- data.raw[index.test, ]

fea.name <- colnames(data.train[1:num.fea])
label.train <- data.train[, pos.label]

# Prepare for the training and test data
if (pos.label != 1){  # the label is after all attributes
  x.train.temp <- t(data.train[, 1:num.fea])  # usually row -> features, transform if not
  x.test.temp <- t(data.test[, 1:num.fea])
} else if (pos.label == 1){  # the label is before all attrbutes
  x.train.temp <- t(data.train[, 2:(1+num.fea)])  # usually row -> features, transform if not
  x.test.temp <- t(data.test[, 2:(1+num.fea)])
} else {
  print('ERROR!! Where is the label in the dataset?!!')
}

# Normalization (Avoid using info from test set)

x.mean <- apply(x.train.temp, 1, mean)
x.sd <- apply(x.train.temp, 1, sd)
x.train <- (x.train.temp - x.mean) / x.sd  # row --> feature, column --> sample
x.test <- (x.test.temp - x.mean) / x.sd
x.train[is.na(x.train)] <- 0
x.test[is.na(x.test)] <- 0

y.train <- data.train[, pos.label]
y.test <- data.test[, pos.label]

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

# =============== Select Features using Ensemble ===============

print('Start selecting features using Ensemble')
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

seed <- 1234
for (i in c(1:num.round)){
  seed <- seed + 1
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

stab.ref <- stab(fea.top.ref, num.round, num.sel.fea)
stab.chs <- stab(fea.top.chs, num.round, num.sel.fea)

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

# use formula to get the ensemble ranking list
fea.rank.merge.ensemble <- data.frame(final.rank = rep(0, num.fea))
row.names(fea.rank.merge.ensemble) <- fea.name
fea.rank.merge.ensemble <- merge(fea.rank.merge.ensemble, fea.rank.merge.ref.order[, 'final.rank', drop = F], by = 'row.names', all = T)
row.names(fea.rank.merge.ensemble) <- fea.rank.merge.ensemble$Row.names
fea.rank.merge.ensemble <- fea.rank.merge.ensemble[, -1]
fea.rank.merge.ensemble <- merge(fea.rank.merge.ensemble, fea.rank.merge.chs.order[, 'final.rank', drop = F], by = 'row.names')
row.names(fea.rank.merge.ensemble) <- fea.rank.merge.ensemble$Row.names
fea.rank.merge.ensemble <- fea.rank.merge.ensemble[, -1]
colnames(fea.rank.merge.ensemble) <- c('final.rank', 'ref.rank', 'chs.rank')

for (line in c(1:num.fea)){
  fea.rank.merge.ensemble$final.rank[line] <- (1 - stab.ref^(2)) * fea.rank.merge.ref$final.rank[line] + 
    (1 - stab.chs^(2)) * fea.rank.merge.chs$final.rank[line]
}
## Pick top features as feature subset
fea.order.ensemble <- fea.rank.merge.ensemble[order(fea.rank.merge.ensemble$final.rank), ]
sel.fea.ensemble.form <- row.names(fea.order.ensemble[1:num.sel.fea, ])

# use RankAggreg to get the ensemble ranking list
rank.list.ref <- row.names(fea.rank.merge.ref[order(fea.rank.merge.ref[, 'final.rank']), ])
rank.list.chs <- row.names(fea.rank.merge.chs[order(fea.rank.merge.chs[, 'final.rank']), ])
rank.matrix <- rbind(rank.list.ref, rank.list.chs)
rank.list.ensemble <- RankAggreg(rank.matrix, num.sel.fea, NULL, method = 'GA', importance = c(stab.ref, stab.chs))
sel.fea.ensemble.consensus <- rank.list.ensemble$top.list

# =============== END of Selecting Features using Ensemble ===============

# =============== Feature Selection using ReliefF ===============

data.ref <- data.frame(t(x.train), y.train, check.names = F)  # add the param to avoid changing '-' to '.'
estReliefF <- attrEval('y.train', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
names(estReliefF) <- fea.name  # It's very annoying that 'attrEval' will change the '-' in the names to '.'
fea.order.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]

# =============== Feature Selection using Chi-Squared ===============

data.chs <- data.frame(t(x.train), y.train, check.names = F)
weights <- chi.squared(y.train~., data.chs)
fea.order.chs <- weights[order(weights$attr_importance, decreasing = T), , drop = F]

# =============== END of Selecting Features ===============

# =============== Prediction ===============

# Prepare the data sets for classificaion, usually row --> sample, col --> feature

train.rgife <- t(x.train[sel.fea.rgife, ])
train.ref <- t(x.train[sel.fea.ref, ])
train.chs <- t(x.train[sel.fea.chs, ])
train.ensemble.form <- t(x.train[sel.fea.ensemble.form, ])
train.ensemble.consensus <- t(x.train[sel.fea.ensemble.consensus, ])

test.rgife <- t(x.test[sel.fea.rgife, ])
test.ref <- t(x.test[sel.fea.ref, ])
test.chs <- t(x.test[sel.fea.chs, ])
test.ensemble.form <- t(x.test[sel.fea.ensemble.form, ])
test.ensemble.consensus <- t(x.test[sel.fea.ensemble.consensus, ])

# Use SVM to test the performance

## Prepare the data sets for SVM

train.rgife.svm <- data.frame(cbind(train.rgife, y.train))
test.rgife.svm <- data.frame(test.rgife)

train.ref.svm <- data.frame(cbind(train.ref, y.train))
test.ref.svm <- data.frame(test.ref)

train.chs.svm <- data.frame(cbind(train.chs, y.train))
test.chs.svm <- data.frame(test.chs)

train.ensemble.form.svm <- data.frame(cbind(train.ensemble.form, y.train))
test.ensemble.form.svm <- data.frame(test.ensemble.form)

train.ensemble.consensus.svm <- data.frame(cbind(train.ensemble.consensus, y.train))
test.ensemble.consensus.svm <- data.frame(test.ensemble.consensus)

## Train SVM classifiers

model.rgife.svm <- svm(formula = factor(y.train)~., data = train.rgife.svm, probability = TRUE)
model.ref.svm <- svm(formula = factor(y.train)~., data = train.ref.svm, probability = TRUE)
model.chs.svm <- svm(formula = factor(y.train)~., data = train.chs.svm, probability = TRUE)
model.ensemble.form.svm <- svm(formula = factor(y.train)~., data = train.ensemble.form.svm, probability = TRUE)
model.ensemble.consensus.svm <- svm(formula = factor(y.train)~., data = train.ensemble.consensus.svm, probability = TRUE)

## Prediction on the test set

pred.rgife.svm <- attr(predict(model.rgife.svm, test.rgife.svm, probability = TRUE), 'probabilities')[, 2]
pred.ref.svm <- attr(predict(model.ref.svm, test.ref.svm, probability = TRUE), 'probabilities')[, 2]
pred.chs.svm <- attr(predict(model.chs.svm, test.chs.svm, probability = TRUE), 'probabilities')[, 2]
pred.ensemble.form.svm <- attr(predict(model.ensemble.form.svm, test.ensemble.form.svm, probability = TRUE), 'probabilities')[, 2]
pred.ensemble.consensus.svm <- attr(predict(model.ensemble.consensus.svm, test.ensemble.consensus.svm, probability = TRUE), 'probabilities')[, 2]

## Test the classification accuracy using the same test set

auc.rgife <- auc(y.test, pred.rgife.svm)
auc.ref <- auc(y.test, pred.ref.svm)
auc.chs <- auc(y.test, pred.chs.svm)
auc.ensemble.form <- auc(y.test, pred.ensemble.form.svm)
auc.ensemble.consensus <- auc(y.test, pred.ensemble.consensus.svm)


