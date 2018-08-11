# EFSIS (Ensemble of Feature Selection Integrating Stability
# Ensemble of SAM, GeoDE, InformationGain and ReliefF, using formula to integrate them. Compare with SAM, GeoDE, InformationGain, ReliefF, and Function Perturbation (using product of rankings for integration)
# Use 10-fold cross-validation to get 10 prediction accuracies
# Inherited from 'efsis.R', use 10-folds cross validation, and 50 rounds for bootstrap resampling, when calculating the stab, use the original # sel fea, no 10 times
# new updates compared with 'efsis.R': 
# 1. use StabPerf to calculate stab
# 2. But normalize stabPerf to [0,1]
# 3. Use the percents of selected features from paper <Pes B. et al., 2017>

remove(list = ls())

library(foreign)
library(caret)
library(samr)
library(GeoDE)
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
percent.sel.fea <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5) / 100
k.folds <- 10  # k-fold cross validation
num.round <- 50  # number of rounds of resampling for EFSIS
seed.10fold <- 12345

# ============ function definition ===========

# define the function of calculating stability
stab <- function(dataFrame){
  num.round <- ncol(dataFrame)  # number of resampling rounds
  union <- Reduce(union, dataFrame)  # F in the formula: the list of all features, which have been selected in at least one of n sampling steps
  whole <- unlist(list(dataFrame))  # all the features in List (with replicate)
  sum <- 0
  for(fea in union){
    freq <- length(which(whole == fea))
    sum <- sum + freq
  }
  stab.na <- sum / (num.round * length(union))
  stab.na <- (stab.na - (1/num.round))/(1 - (1/num.round))
  return(stab.na)
}

# # Kuncheva consistency index
# stab <- function(dataFrame){
#   m <- ncol(dataFrame)  # number of vectors to be compared
#   d <- num.fea  # number of total features
#   k <- length(dataFrame[ ,1])  # # selected features
#   index.sum <- 0
#   for (i in c(1:(m - 1))){
#     for (j in c((i + 1):m)){
#       x <- dataFrame[ ,i]
#       y <- dataFrame[ ,j]
#       r <- length(intersect(x,y))
#       index <- (r - ((k^2)/d))/(k - ((k^2)/d))
#       # scale to [0,1] range
#       # (x - xmin)/(xmax - xmin)
#       index <- (index - -1)/(1 - -1)
#       index.sum <- index.sum + index
#     }
#   }
#   stab.kuncheva <- 2 * index.sum / (m * (m - 1))
#   return(stab.kuncheva)
# }

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
  
  fea.rank.merge.infog <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.infog) <- fea.name
  
  # the dataframe collecting the top features for all rounds
  
  fea.top.sam <- data.frame()
  fea.top.geode <- data.frame()
  fea.top.ref <- data.frame()
  fea.top.infog <- data.frame()
  
  ## score the features based on each round of resampling
  
  seed.resample <- 1234
  for (i in c(1:num.round)){
    seed.resample <- seed.resample + 1
    set.seed(seed = seed.resample)
    id.train.resample <- c(sample(x = which(label.train == levels(label.train)[1]), size = num.resample.control, replace = T), sample(x = which(label.train == levels(label.train)[2]), size = num.resample.treat, replace = T))
    x.train.resample <- x.train[, id.train.resample]
    y.train.resample <- y.train[id.train.resample]
    
    # use SAM to rank the features
    label.sam <- factor(y.train.resample)
    levels(label.sam) <- c('1', '2')
    data.sam <- list(x=x.train.resample, y=factor(label.sam), genenames=paste("gene",as.character(fea.name),sep=""), geneid=as.character(fea.name), logged2=TRUE)  
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
    fea.top.sam <- cbind.all(fea.top.sam, fea.rank.name.sam[1:(num.sel.fea)])  # add the top features from this round to the whole record
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
    fea.top.geode <- cbind.all(fea.top.geode, fea.rank.name.geode[1:(num.sel.fea)])  # all the top features for this round to the whole record
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
    fea.top.ref <- cbind.all(fea.top.ref, fea.rank.name.ref[1:(num.sel.fea)])  # add the top features for this round to the whole record
    colnames(fea.rank.ref) <- 'rank'
    fea.rank.ref['rank'] <- seq(1, num.fea)
    fea.rank.merge.ref <- merge(fea.rank.merge.ref, fea.rank.ref, by = 'row.names')
    row.names(fea.rank.merge.ref) <- fea.rank.merge.ref$Row.names
    fea.rank.merge.ref <- fea.rank.merge.ref[, -1]
    
    # use Information Gain to rank the features
    data.infog <- data.frame(t(x.train.resample), y.train.resample, check.names = F)
    weights <- information_gain(y.train.resample~., data.infog)
    fea.rank.infog <- weights[order(weights$importance, decreasing = T), , drop = F]
    fea.rank.name.infog <- fea.rank.infog$attributes  # the ranked feature list for this round
    fea.top.infog <- cbind.all(fea.top.infog, fea.rank.name.infog[1:(num.sel.fea)])  # add the top features for this round to the whole record
    fea.rank.infog <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.infog) <- fea.rank.name.infog
    fea.rank.merge.infog <- merge(fea.rank.merge.infog, fea.rank.infog, by = 'row.names')
    row.names(fea.rank.merge.infog) <- fea.rank.merge.infog$Row.names
    fea.rank.merge.infog <- fea.rank.merge.infog[, -1]
  }
  
  # =============== Stability Calculation ===============
  
  stab.sam <- stab(fea.top.sam)
  stab.geode <- stab(fea.top.geode)
  stab.ref <- stab(fea.top.ref)
  stab.infog <- stab(fea.top.infog)
  
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
  
  ### sub-final score for infog
  for (line in c(1:num.fea)){
    fea.rank.merge.infog$final.rank[line] <- prod(fea.rank.merge.infog[line, ][-1])
  }
  fea.rank.merge.infog.order <- fea.rank.merge.infog[order(fea.rank.merge.infog$final.rank), , drop = F]
  fea.rank.merge.infog.order$final.rank <- seq(1, num.fea)
  
  # =============== Final Score Calculation ===============
  
  # use formula to get the EFSIS ranking list
  fea.rank.merge.efsis <- data.frame(final.rank = rep(0, num.fea))
  row.names(fea.rank.merge.efsis) <- fea.name
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.sam.order[, 'final.rank', drop = F], by = 'row.names', all = T)
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.geode.order[, 'final.rank', drop = F], by = 'row.names')
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.ref.order[, 'final.rank', drop = F], by = 'row.names', all = T)
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  fea.rank.merge.efsis <- merge(fea.rank.merge.efsis, fea.rank.merge.infog.order[, 'final.rank', drop = F], by = 'row.names')
  row.names(fea.rank.merge.efsis) <- fea.rank.merge.efsis$Row.names
  fea.rank.merge.efsis <- fea.rank.merge.efsis[, -1]
  colnames(fea.rank.merge.efsis) <- c('final.rank', 'sam.rank', 'geode.rank', 'ref.rank', 'infog.rank')
  
  for (line in c(1:num.fea)){
    # integrage SAM, GeoDE, Ref, infog
    rank.stab <- rank(-c(stab.sam, stab.geode, stab.ref, stab.infog))
    rank.stab.sam <- rank.stab[1]
    rank.stab.geode <- rank.stab[2]
    rank.stab.ref <- rank.stab[3]
    rank.stab.infog <- rank.stab[4]
    fea.rank.merge.efsis$final.rank.sam.geode.ref.infog[line] <- (fea.rank.merge.efsis$sam.rank[line])^rank.stab.sam * 
      (fea.rank.merge.efsis$geode.rank[line])^rank.stab.geode * (fea.rank.merge.efsis$ref.rank[line])^rank.stab.ref * 
      (fea.rank.merge.efsis$infog.rank[line])^rank.stab.infog
  }
  ## Pick top features as feature subset
  # integrage SAM, GeoDE, Ref, infog
  fea.order.efsis.sam.geode.ref.infog <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank.sam.geode.ref.infog), ]
  sel.fea.efsis <- row.names(fea.order.efsis.sam.geode.ref.infog[1:num.sel.fea, ])
  
  # =============== END of Selecting Features using efsis ===============
  output.list <- list('sel.fea.efsis' = sel.fea.efsis, 
                      'stab.sam' = stab.sam, 'stab.geode' = stab.geode, 'stab.ref' = stab.ref, 'stab.infog' = stab.infog)
}

# =============== Data Preparation ===============

# Load the data
data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature

# Get the general information about this dataset
num.sample <- nrow(data.raw)  # number of samples
num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
nums.sel.fea <- ceiling(num.fea * percent.sel.fea)
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
  sel.fea.infog.folds <- data.frame()
  sel.fea.func.folds <- data.frame()
  sel.fea.efsis.folds <- data.frame()
  
  auc.all.folds <- data.frame(row.names = c('sam', 'geode', 'ref', 'infog', 'func', 'EFSIS'))
  time.all.folds <- data.frame(row.names = c('sam', 'geode', 'ref', 'infog', 'func', 'EFSIS'))
  
  for (i in c(1:k.folds)){
    
    # =============== Data Preparation ===============
    
    index.train <- c(index.control[pos.control.train.list[[i]]], index.treat[pos.treat.train.list[[i]]])
    index.train.control <- index.control[pos.control.train.list[[i]]]
    index.train.treat <- index.treat[pos.treat.train.list[[i]]]
    index.test <- setdiff(c(index.control, index.treat), index.train)
    data.train <- data.raw[index.train, ]  # row -> sample, column -> feature
    data.test <- data.raw[index.test, ]
    
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
    
    start.time.func <- proc.time()  # Starting time for function perturbation (include the time taken by each individual method)
    # =============== Feature Selection using SAM ===============
    
    # print('Start selecting features using SAM')
    start.time <- proc.time()
    label.sam <- factor(y.train)
    levels(label.sam) <- c('1', '2')
    data.sam <- list(x=x.train, y=factor(label.sam), genenames=paste("gene",as.character(fea.name),sep=""), geneid=as.character(fea.name), logged2=TRUE)  
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
    time.taken.sam <- end.time - start.time
    #     cat('Time taken by SAM:', '\n')
    #     cat(time.taken[1], '\n')
    
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
    time.taken.geode <- end.time - start.time
    #     cat('Time taken by GeoDE:', '\n')
    #     cat(time.taken[1], '\n')
    
    # =============== Feature Selection using ReliefF ===============
    
    # print('Start selecting features using ReliefF')
    start.time <- proc.time()
    data.ref <- data.frame(t(x.train), y.train, check.names = F)  # add the param to avoid changing '-' to '.'
    estReliefF <- attrEval('y.train', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = 30)
    names(estReliefF) <- fea.name  # It's very annoying that 'attrEval' will change the '-' in the names to '.'
    fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
    fea.rank.name.ref <- names(fea.rank.ref)
    sel.fea.ref <- fea.rank.name.ref[1:num.sel.fea]
    end.time <- proc.time()
    time.taken.ref <- end.time - start.time
    #     cat('Time taken by ReliefF:', '\n')
    #     cat(time.taken[1], '\n')
    
    # =============== Feature Selection using Information Gain ===============
    
    # print('Start selecting features using Information Gain')
    start.time <- proc.time()
    data.infog <- data.frame(t(x.train), y.train, check.names = F)
    weights <- information_gain(y.train~., data.infog)
    fea.rank.infog <- weights[order(weights$importance, decreasing = T), , drop = F]
    fea.rank.name.infog <- fea.rank.infog$attributes
    sel.fea.infog <- fea.rank.name.infog[1:num.sel.fea]
    end.time <- proc.time()
    time.taken.infog <- end.time - start.time
    #     cat('Time taken by Information Gain:', '\n')
    #     cat(time.taken[1], '\n')
    
    # =============== Feature Selection using Function Perturbation ===============
    # print('Start selecting features using Function Perturbation')
    start.time <- start.time.func
    ##  The ranking lists for 4 individual ranking methods
    fea.rank.list.sam <- data.frame(rank.sam = seq(1, num.fea))
    row.names(fea.rank.list.sam) <- fea.rank.name.sam
    fea.rank.list.geode <- data.frame(rank.geode = seq(1, num.fea))
    row.names(fea.rank.list.geode) <- fea.rank.name.geode
    fea.rank.list.ref <- data.frame(rank.ref = seq(1, num.fea))
    row.names(fea.rank.list.ref) <- fea.rank.name.ref
    fea.rank.list.infog <- data.frame(rank.infog = seq(1, num.fea))
    row.names(fea.rank.list.infog) <- fea.rank.name.infog
    ## Merge the 4 lists together
    fea.rank.merge.func <- data.frame(final.rank = rep(0, num.fea))
    row.names(fea.rank.merge.func) <- fea.name
    fea.rank.merge.func <- merge(fea.rank.merge.func, fea.rank.list.sam[, 'rank.sam', drop = F], by = 'row.names', all = T)
    row.names(fea.rank.merge.func) <- fea.rank.merge.func$Row.names
    fea.rank.merge.func <- fea.rank.merge.func[, -1]
    fea.rank.merge.func <- merge(fea.rank.merge.func, fea.rank.list.geode[, 'rank.geode', drop = F], by = 'row.names', all = T)
    row.names(fea.rank.merge.func) <- fea.rank.merge.func$Row.names
    fea.rank.merge.func <- fea.rank.merge.func[, -1]
    fea.rank.merge.func <- merge(fea.rank.merge.func, fea.rank.list.ref[, 'rank.ref', drop = F], by = 'row.names', all = T)
    row.names(fea.rank.merge.func) <- fea.rank.merge.func$Row.names
    fea.rank.merge.func <- fea.rank.merge.func[, -1]
    fea.rank.merge.func <- merge(fea.rank.merge.func, fea.rank.list.infog[, 'rank.infog', drop = F], by = 'row.names', all = T)
    row.names(fea.rank.merge.func) <- fea.rank.merge.func$Row.names
    fea.rank.merge.func <- fea.rank.merge.func[, -1]
    ## Calculate the product of 4 ranking lists
    for (line in c(1:num.fea)){
      # integrage SAM, GeoDE, Ref, infog
      fea.rank.merge.func$final.rank[line] <- (fea.rank.merge.func$rank.sam[line]) / 100 * (fea.rank.merge.func$rank.geode[line]) / 100 * (fea.rank.merge.func$rank.ref[line]) /100 * (fea.rank.merge.func$rank.infog[line]) / 100
    }
    ## Pick top features as feature subset
    fea.order.func <- fea.rank.merge.func[order(fea.rank.merge.func$final.rank), ]
    sel.fea.func <- row.names(fea.order.func[1:num.sel.fea, ])
    end.time <- proc.time()
    time.taken.func <- end.time - start.time
    #     cat('Time taken by Function Perturbation:', '\n')
    #     cat(time.taken[1], '\n')
    
    # =============== Feature Selection using EFSIS ===============
    
    # print('Start selecting features using EFSIS')
    start.time <- proc.time()
    num.resample.control <- length(index.train.control)  # using bootsrap, so number of sampled samples in each round equals the original #
    num.resample.treat <- length(index.train.treat)
    output.list.efsis <- efsis()
    sel.fea.efsis <- output.list.efsis$sel.fea.efsis
    end.time <- proc.time()
    time.taken.efsis <- end.time - start.time
    #     cat('Time taken by EFSIS:', '\n')
    #     cat(time.taken[1], '\n')
    stab.sam <- output.list.efsis$stab.sam
    stab.geode <- output.list.efsis$stab.geode
    stab.ref <- output.list.efsis$stab.ref
    stab.infog <- output.list.efsis$stab.infog
    
    # save the selected features
    
    sel.fea.sam.folds <- cbind.all(sel.fea.sam.folds, sel.fea.sam)
    sel.fea.geode.folds <- cbind.all(sel.fea.geode.folds, sel.fea.geode)
    sel.fea.ref.folds <- cbind.all(sel.fea.ref.folds, sel.fea.ref)
    sel.fea.infog.folds <- cbind.all(sel.fea.infog.folds, sel.fea.infog)
    sel.fea.func.folds <- cbind.all(sel.fea.func.folds, sel.fea.func)
    sel.fea.efsis.folds <- cbind.all(sel.fea.efsis.folds, sel.fea.efsis)
    
    time.all <- c(time.taken.sam, time.taken.geode, time.taken.ref, time.taken.infog, time.taken.func, time.taken.efsis)
    time.all.folds <- cbind.all(time.all.folds, time.all)
    
    print('Finished selecting features')
    
    # =============== END of Selecting Features ===============
    
    # =============== Prediction ===============
    print('Start prediction')
    
    # Prepare the data sets for classificaion, usually row --> sample, col --> feature
    
    train.sam <- t(x.train[sel.fea.sam, ])
    train.geode <- t(x.train[sel.fea.geode, ])
    train.ref <- t(x.train[sel.fea.ref, ])
    train.infog <- t(x.train[sel.fea.infog, ])
    train.func <- t(x.train[sel.fea.func, ])
    train.efsis <- t(x.train[sel.fea.efsis, ])
    
    test.sam <- t(x.test[sel.fea.sam, ])
    test.geode <- t(x.test[sel.fea.geode, ])
    test.ref <- t(x.test[sel.fea.ref, ])
    test.infog <- t(x.test[sel.fea.infog, ])
    test.func <- t(x.test[sel.fea.func, ])
    test.efsis <- t(x.test[sel.fea.efsis, ])
    
    # Use SVM to test the performance
    
    ## Prepare the data sets for SVM
    
    train.sam.svm <- data.frame(cbind(train.sam, y.train))
    test.sam.svm <- data.frame(test.sam)
    
    train.geode.svm <- data.frame(cbind(train.geode, y.train))
    test.geode.svm <- data.frame(test.geode)
    
    train.ref.svm <- data.frame(cbind(train.ref, y.train))
    test.ref.svm <- data.frame(test.ref)
    
    train.infog.svm <- data.frame(cbind(train.infog, y.train))
    test.infog.svm <- data.frame(test.infog)
    
    train.func.svm <- data.frame(cbind(train.func, y.train))
    test.func.svm <- data.frame(test.func)
    
    train.efsis.svm <- data.frame(cbind(train.efsis, y.train))
    test.efsis.svm <- data.frame(test.efsis)
    
    
    ## Train SVM classifiers
    
    model.sam.svm <- svm(formula = factor(y.train)~., data = train.sam.svm, probability = TRUE)
    model.geode.svm <- svm(formula = factor(y.train)~., data = train.geode.svm, probability = TRUE)
    model.ref.svm <- svm(formula = factor(y.train)~., data = train.ref.svm, probability = TRUE)
    model.infog.svm <- svm(formula = factor(y.train)~., data = train.infog.svm, probability = TRUE)
    model.func.svm <- svm(formula = factor(y.train)~., data = train.func.svm, probability = TRUE)
    model.efsis.svm <- svm(formula = factor(y.train)~., data = train.efsis.svm, probability = TRUE)
    
    ## Prediction on the test set
    
    pred.sam.svm <- attr(predict(model.sam.svm, test.sam.svm, probability = TRUE), 'probabilities')[, 2]
    pred.geode.svm <- attr(predict(model.geode.svm, test.geode.svm, probability = TRUE), 'probabilities')[, 2]
    pred.ref.svm <- attr(predict(model.ref.svm, test.ref.svm, probability = TRUE), 'probabilities')[, 2]
    pred.infog.svm <- attr(predict(model.infog.svm, test.infog.svm, probability = TRUE), 'probabilities')[, 2]
    pred.func.svm <- attr(predict(model.func.svm, test.func.svm, probability = TRUE), 'probabilities')[, 2]
    pred.efsis.svm <- attr(predict(model.efsis.svm, test.efsis.svm, probability = TRUE), 'probabilities')[, 2]
    
    ## Test the classification accuracy using the same test set
    
    auc.sam <- auc(y.test, pred.sam.svm)
    auc.geode <- auc(y.test, pred.geode.svm)
    auc.ref <- auc(y.test, pred.ref.svm)
    auc.infog <- auc(y.test, pred.infog.svm)
    auc.func <- auc(y.test, pred.func.svm)
    auc.efsis <- auc(y.test, pred.efsis.svm)
    auc.all <- c(auc.sam, auc.geode, auc.ref, auc.infog, auc.func, auc.efsis)
    auc.all.folds <- cbind.all(auc.all.folds, auc.all)
    
    print(paste('Finished prediction and the', i, '/10 fold'))
  }
  
  # save the selected features for this #-sel-fea to file
  sel.fea.all.folds <- list(sam = sel.fea.sam.folds, geode = sel.fea.geode.folds, ref = sel.fea.ref.folds,
                            infog = sel.fea.infog.folds, func = sel.fea.func.folds, efsis = sel.fea.efsis.folds)
  save(sel.fea.all.folds, file = paste(path.data, 'num', num.sel.fea, '-', 'sel-fea-efsis_stabPerf.RData', sep = ''))
  
  # save the auc for this #-sel-fea to file
  write.table(auc.all.folds, paste(path.data, 'num', num.sel.fea, '-', 'auc-efsis_stabPerf.txt', sep = ''), quote = F, col.names = T, row.names = T)
  
  # save the time for this #-sel-fea to file
  write.table(time.all.folds, paste(path.data, 'num', num.sel.fea, '-', 'time-efsis_stabPerf.txt', sep = ''), quote = F, col.names = T, row.names = T)
  
  # calculate the stability after 10 folds
  stab.sam <- stab(sel.fea.sam.folds)
  stab.geode <- stab(sel.fea.geode.folds)
  stab.ref <- stab(sel.fea.ref.folds)
  stab.infog <- stab(sel.fea.infog.folds)
  stab.func <- stab(sel.fea.func.folds)
  stab.efsis <- stab(sel.fea.efsis.folds)
  
  stab.all <- c(stab.sam, stab.geode, stab.ref, stab.infog, stab.func, stab.efsis)
  # save the stability for this #-sel-fea to file
  write.table(stab.all, paste(path.data, 'num', num.sel.fea, '-', 'stab-efsis_stabPerf.txt', sep = ''), quote = F, col.names = F, row.names = F)
}

