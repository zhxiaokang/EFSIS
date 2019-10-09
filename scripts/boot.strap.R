# define the function of boot.strap

boot.strap <- function(seed.resample, round.btsp, label.train, num.resample.control, num.resample.treat, x.train, y.train, fea.name) {
  seed.resample <- seed.resample + round.btsp
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
  
  # use GeoDE to rank the features
  gammas <- 1
  data.geode <- data.frame(fea.name, x.train.resample)
  label.geode <- factor(y.train.resample)
  levels(label.geode) <- c('1', '2')  # GeoDE::chdirAnalysis only accepts '1' and '2' as label factors
  chdir.analysis <- chdirAnalysis(data.geode, factor(label.geode), gammas, CalculateSig=FALSE, nnull=10)
  fea.rank.name.geode <- names(chdir.analysis$results[[1]])  # the ranked feature list for this round
  # print(paste('geode', str(fea.rank.name.geode)))
  
  # use ReliefF to rank the features
  data.ref <- data.frame(t(x.train.resample), y.train.resample, check.names = F)  # add the param to avoid changing '-' to '.'
  estReliefF <- attrEval('y.train.resample', data.ref, estimator = 'ReliefFexpRank', ReliefIterations = -2, maxThreads = 1)
  names(estReliefF) <- fea.name  # This command needs to be added because it's very annoying that 'attrEval' will change the '-' in the names to '.'
  fea.rank.ref <- estReliefF[order(abs(estReliefF), decreasing = T)]
  fea.rank.ref <- data.frame(importance = fea.rank.ref)
  fea.rank.name.ref <- rownames(fea.rank.ref)  # the ranked feature list for this round
  
  # use Information Gain to rank the features
  data.infog <- data.frame(t(x.train.resample), y.train.resample, check.names = F)
  weights <- information_gain(y.train.resample~., data.infog)
  fea.rank.infog <- weights[order(weights$importance, decreasing = T), , drop = F]
  fea.rank.name.infog <- fea.rank.infog$attributes  # the ranked feature list for this round
  # print(paste('infog', str(fea.rank.name.infog)))
  result.this.round <- list(sam = fea.rank.name.sam, geode = fea.rank.name.geode, ref = fea.rank.name.ref, infog = fea.rank.name.infog)
  return(result.this.round)
}