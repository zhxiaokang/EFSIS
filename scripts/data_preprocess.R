# Data preprocess to convert all data into .arff

library(foreign)
# Define the function to add columns to empty dataframe
cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# Dlbcl
dlbcl.raw <- read.table('../data/DLBCL/dlbcl_preprocessed.txt', row.names = 1)
dlbcl.df <- as.data.frame.matrix(t(dlbcl.raw))
dlbcl.df$ydlbcl <- factor(dlbcl.df$ydlbcl)
write.arff(dlbcl.df, '../data/DLBCL/dlbcl.arff')

# CNS
cns.raw <- read.arff('../data/CNS/cns')
write.arff(cns.raw, '../data/CNS/cns.arff')

# Leukemia
leukemia.train <- read.arff('../data/Leukemia/train.arff')
leukemia.test <- read.arff('../data/Leukemia/test.arff')
leukemia.raw <- rbind(leukemia.train, leukemia.test)
write.arff(leukemia.raw, '../data/Leukemia/leukemia.arff')


# Prostate Singh
prostate.singh.feature <- read.table('../data/ProstateSingh/prostate_preprocessed.txt', row.names = 1, nrows = 2135)
prostate.singh.class <- read.table('../data/ProstateSingh/prostate_preprocessed.txt', row.names = 1, skip = 2135)
prostate.singh.feature.df <- as.data.frame.matrix(t(prostate.singh.feature))
prostate.singh.class.df <- as.data.frame.matrix(t(prostate.singh.class))
prostate.singh.df <- cbind(prostate.singh.feature.df, prostate.singh.class.df)
write.arff(prostate.singh.df, '../data/ProstateSingh/prostate.singh.arff')

# Breast
breast.feature <- read.table('../data/Breast/breast_preprocessed.txt', row.names = 1, nrows = 47293)
breast.class <- read.table('../data/Breast/breast_preprocessed.txt', row.names = 1, skip = 47293)
breast.feature.df <- as.data.frame.matrix(t(breast.feature))
breast.class.df <- as.data.frame.matrix(t(breast.class))
breast.df <- cbind(breast.feature.df, breast.class.df)
write.arff(breast.df, '../data/Breast/breast.arff')

# AML
aml.header <- read.table('../data/AML/AMLGSE2191.tab', nrows = 1, header = FALSE, stringsAsFactors = FALSE)
aml.raw <- read.table('../data/AML/AMLGSE2191.tab', fill = T, skip = 3, header = FALSE)
colnames(aml.raw) <- unlist(aml.header)
aml.clean <- aml.raw[, -ncol(aml.raw)]
write.arff(aml.clean, '../data/AML/aml.arff')

# ColonBreast
colonbreast.header <- read.table('../data/ColonBreast/BC_CCGSE3726_frozen.txt', nrows = 1, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colonbreast.raw <- read.table('../data/ColonBreast/BC_CCGSE3726_frozen.txt', sep = '\t', fill = T, skip = 3, header = FALSE)
colnames(colonbreast.raw) <- unlist(colonbreast.header)
colonbreast.clean <- colonbreast.raw[, -ncol(colonbreast.raw)]
write.arff(colonbreast.clean, '../data/ColonBreast/colonbreast.arff')

# ProstateSboner
file.list <- list.files(path = '../data/ProstateSboner/E-GEOD-16560.processed.1', pattern = '.txt')
prostate.sboner.feature <- data.frame()
for (sample in file.list) {
  expression <- read.table(paste('../data/ProstateSboner/E-GEOD-16560.processed.1/', sample, sep = ''), header = T, row.names = 1, sep = '\t')
  colnames(expression) <- sample
  prostate.sboner.feature <- cbind.all(prostate.sboner.feature, expression)
}
prostate.sboner.feature <- t(prostate.sboner.feature)
prostate.sboner.class <- read.table('../data/ProstateSboner/name.label.uniq', header = F, row.names = 1)
colnames(prostate.sboner.class) <- 'class'
rownames(prostate.sboner.feature) <- rownames(prostate.sboner.class)
prostate.sboner <- cbind(prostate.sboner.feature, prostate.sboner.class)
write.arff(prostate.sboner, '../data/ProstateSboner/prostate.sboner.arff')

