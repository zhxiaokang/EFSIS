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