# Compare the two stability calsulation methods

cbind.all <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

stab1 <- function(dataFrame){
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
  # stab <- max(0, stab.na - 10 * (num.sel.fea/num.fea))
  return(stab.na)
}

stab2 <- function(dataFrame){
  m <- ncol(dataFrame)  # number of vectors to be compared
  d <- num.fea  # number of total features
  k <- length(dataFrame[ ,1])  # # selected features
  index.sum <- 0
  for (i in c(1:(m - 1))){
    for (j in c((i + 1):m)){
      x <- dataFrame[ ,i]
      y <- dataFrame[ ,j]
      r <- length(intersect(x,y))
      index <- (r - ((k^2)/d))/(k - ((k^2)/d))
      # scale to [0,1] range
      # (x - xmin)/(xmax - xmin)
      index <- (index - -1)/(1 - -1)
      index.sum <- index.sum + index
    }
  }
  stab.kuncheva <- 2 * index.sum / (m * (m - 1))
  return(stab.kuncheva)
}


full <- seq(1,10000)

num.fea <- length(full)
num.vector <- 10
# num.sel.fea <- 100

stabs1 <- c()
stabs2 <- c()

for (num.sel.fea in seq(10,500,10)){
  df.sel <- data.frame()
  for (i in c(1:num.vector)){
    df.sel <- cbind.all(df.sel, sample(full, num.sel.fea))
  }
  
  stabs1 <- c(stabs1, stab1(df.sel))
  stabs2 <- c(stabs2, stab2(df.sel))
}

