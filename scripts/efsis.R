# define the function of EFSIS
efsis <- function(num.fea, fea.name, num.round, num.cores){
  
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
  rounds.btsp <- seq(1, num.round)  # sequence of numbers for bootstrap, used for parallel computing
  seed.resample <- 1234
  
  source('boot.strap.R')  # the function of bootstrap
  
  # Use parallel computing for 50 rounds
  ## Use foreach
  #   cl <- makeForkCluster(num.cores)
  #   registerDoParallel(cl)
  #   results.btsp <- foreach(round.btsp = rounds.btsp) %dopar% boot.strap(round.btsp)
  #   stopCluster(cl)
  #   registerDoSEQ()
  ## Use 'parallel' - parLapply
  #   cl <- makeCluster(num.cores, type = 'FORK')
  #   results.btsp <- parLapply(cl, rounds.btsp, boot.strap)
  #   stopCluster(cl)
  print('start parallel computing')
  results.btsp <- mclapply(rounds.btsp, boot.strap, mc.preschedule = FALSE, mc.cores = num.cores)
  
  # Check whether the output is correct, otherwise repeat this multi-thread processing
  round.count <- 1
  num.redo.paral <- 0
  while (round.count <= num.round){
    if (length(results.btsp[[round.count]]) == 4) {
      round.count <- round.count + 1
    }
    else {
      round.count <- 1
      num.redo.paral <- num.redo.paral + 1
      results.btsp <- mclapply(rounds.btsp, boot.strap, mc.cores = num.cores, mc.preschedule = FALSE)
    }
  }
  
  print(paste('Times of redo parallel computing:', num.redo.paral))
  
  for (i in c(1:num.round)){
    fea.rank.name.sam <- results.btsp[[i]]$sam
    # print(paste('sam', str(fea.rank.name.sam)))
    fea.rank.name.geode <- results.btsp[[i]]$geode
    # print(paste('geode', str(fea.rank.name.geode)))
    fea.rank.name.ref <- results.btsp[[i]]$ref
    # print(paste('ref', str(fea.rank.name.ref)))
    fea.rank.name.infog <- results.btsp[[i]]$infog
    # print(paste('infog', str(fea.rank.name.infog)))
    
    fea.top.sam <- cbind.all(fea.top.sam, fea.rank.name.sam[1:(num.sel.fea)])  # add the top features from this round to the whole record
    fea.rank.sam <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.sam) <- fea.rank.name.sam
    fea.rank.merge.sam <- merge(fea.rank.merge.sam, fea.rank.sam, by = 'row.names')
    row.names(fea.rank.merge.sam) <- fea.rank.merge.sam$Row.names
    fea.rank.merge.sam <- fea.rank.merge.sam[, -1]
    
    fea.top.geode <- cbind.all(fea.top.geode, fea.rank.name.geode[1:(num.sel.fea)])  # all the top features for this round to the whole record
    fea.rank.geode <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.geode) <- fea.rank.name.geode
    fea.rank.merge.geode <- merge(fea.rank.merge.geode, fea.rank.geode, by = 'row.names')
    row.names(fea.rank.merge.geode) <- fea.rank.merge.geode$Row.names
    fea.rank.merge.geode <- fea.rank.merge.geode[, -1]
    
    fea.top.ref <- cbind.all(fea.top.ref, fea.rank.name.ref[1:(num.sel.fea)])  # add the top features for this round to the whole record
    fea.rank.ref <- data.frame(rank = seq(1, num.fea))
    row.names(fea.rank.ref) <- fea.rank.name.ref
    fea.rank.merge.ref <- merge(fea.rank.merge.ref, fea.rank.ref, by = 'row.names')
    row.names(fea.rank.merge.ref) <- fea.rank.merge.ref$Row.names
    fea.rank.merge.ref <- fea.rank.merge.ref[, -1]
    
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
  
  fea.rank.merge.efsos <- fea.rank.merge.efsis
  
  for (line in c(1:num.fea)){
    # integrage SAM, GeoDE, Ref, infog
    
    fea.rank.merge.efsis$final.rank.sam.geode.ref.infog[line] <- (fea.rank.merge.efsis$sam.rank[line])^(1-stab.sam)/100 * 
      (fea.rank.merge.efsis$geode.rank[line])^(1-stab.geode)/100 * (fea.rank.merge.efsis$ref.rank[line])^(1-stab.ref)/100 * 
      (fea.rank.merge.efsis$infog.rank[line])^(1-stab.infog)/100
    
    fea.rank.merge.efsos$final.rank.sam.geode.ref.infog[line] <- (fea.rank.merge.efsos$sam.rank[line])/100 * 
      (fea.rank.merge.efsos$geode.rank[line])/100 * (fea.rank.merge.efsos$ref.rank[line])/100 * 
      (fea.rank.merge.efsos$infog.rank[line])/100
  }
  ## Pick top features as feature subset
  # integrage SAM, GeoDE, Ref, infog
  fea.order.efsis.sam.geode.ref.infog <- fea.rank.merge.efsis[order(fea.rank.merge.efsis$final.rank.sam.geode.ref.infog), ]
  sel.fea.efsis <- row.names(fea.order.efsis.sam.geode.ref.infog[1:num.sel.fea, ])
  
  fea.order.efsos.sam.geode.ref.infog <- fea.rank.merge.efsos[order(fea.rank.merge.efsos$final.rank.sam.geode.ref.infog), ]
  sel.fea.efsos <- row.names(fea.order.efsos.sam.geode.ref.infog[1:num.sel.fea, ])
  
  # =============== END of Selecting Features using efsis ===============
  output.list <- list('sel.fea.efsis' = sel.fea.efsis, 'sel.fea.efsos' = sel.fea.efsos, 
                      'stab.sam' = stab.sam, 'stab.geode' = stab.geode, 'stab.ref' = stab.ref, 'stab.infog' = stab.infog)
}