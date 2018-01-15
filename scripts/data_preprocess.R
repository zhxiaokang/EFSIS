# Data preprocess to convert all data into .arff

library(foreign)

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

# AML

