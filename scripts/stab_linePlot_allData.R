# Arrange all plots of stability from all datasets in one figure

library(foreign)
library(reshape2)
library(ggplot2)
library(scmamp)
library(ggpubr)

cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

plot.stab.line <- function(data.set, script.version) {
  num.folds <- 10
  
  # Load the data
  path.data <- paste('../data/', data.set, '/', sep = '')  # path to the data
  data.file <- list.files(path = path.data, pattern = '.arff')
  data.raw <- read.arff(paste(path.data, data.file, sep = ''))  # row -> sample, column -> feature
  percent.sel.fea <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5) / 100
  
  # Get the general information about this dataset
  
  num.fea <- ncol(data.raw) - 1  # number of features, but notice that the last column is the label
  nums.sel.fea <- ceiling(num.fea * percent.sel.fea)  # the numbers of selected features
  
  df.merge <- data.frame()  # to record the final merged dataframe
  names.methods <- c('SAM', 'GeoDE', 'ReliefF', 'Information Gain', 'Function Perturbation','EFSIS')
  stabs <- matrix(nrow = length(nums.sel.fea), ncol = length(names.methods))
  j <- 0
  for (i in nums.sel.fea) {
    j = j + 1
    file.name <- paste(path.data, 'num', i, '-stab-', script.version, '.txt', sep = '')
    stab <- read.table(file.name)[c(1:5, 7), ]  # read the file, skip the line for 'EFSOS'
    stabs[j, ] <- unlist(stab)
  }
  
  df.stabs <- as.data.frame.matrix(stabs)
  colnames(df.stabs) <- names.methods
  # df.stabs$num.sel.fea <- nums.sel.fea
  df.stabs$percent.sel.fea <- percent.sel.fea * 100
  stabs.long <- melt(df.stabs, id = 'percent.sel.fea')
  cols <- c('grey75', 'grey69', 'grey62', 'grey55', 'grey30', 'grey0')
  linetype <- c('dotted', 'dotdash', 'dashed', 'twodash', 'longdash', 'solid')
  pic <- ggplot(data = stabs.long, aes(x = percent.sel.fea, y = value, linetype = variable, shape = variable, colour = variable)) + 
    scale_linetype_manual(values = linetype) +
    scale_color_manual(values = cols) + 
    geom_line(size = 1, aes(colour = factor(variable))) + geom_point(size = 2.5) + labs(x = 'Percentage of selected features (%)', y = 'Stability', title = data.set) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
                       axis.text=element_text(size=13), axis.title=element_text(size=13, face = 'bold'), 
                       legend.text=element_text(size=12), legend.title = element_blank())
  return(pic)
}

script.version <- 'efsosParalStab'

pic.AML <- plot.stab.line('AML', script.version)
pic.CNS <- plot.stab.line('CNS', script.version)
pic.DLBCL <- plot.stab.line('DLBCL', script.version)
pic.ProstateSingh <- plot.stab.line('ProstateSingh', script.version)
pic.Leukemia <- plot.stab.line('Leukemia', script.version)
pic.ColonBreast <- plot.stab.line('ColonBreast', script.version)

pic <- ggarrange(pic.AML, pic.CNS, pic.DLBCL, pic.ProstateSingh, 
          pic.Leukemia, pic.ColonBreast, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")

pdf(file = paste('../fig/stab-AllData-', script.version, '.pdf', sep = ''), width = 10, height = 12)
pic
dev.off()
