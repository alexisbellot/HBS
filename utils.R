adaptive_jump = function(current_jump, jump.mat,wei){
  # take 50 last jump probabilities
  jump_prob <- jump.mat[(wei-48):wei]
  new_jump  <- ifelse(mean(jump_prob)<0.40, max(current_jump - min(0.05,1/(wei %/% 50)),0.01),
                      current_jump + min(0.05,1/(wei %/% 50)))
  return(new_jump)
}


xval = function(data,n.group=2,k = 5) {
  
  Y = subset(data,select=c(Survival,Status))
  X = subset(data,select=-c(Survival,Status))
  n = nrow(data)
  folds = split(sample(n), seq_len(k))
  
  xval.fold = function(fold) {
    
    unique.times <- unique(sort(data$Survival[-fold]))
    
    ggITE <- HBS(X.train=X[-fold,],Y.train=Y[-fold,],X.test=X[fold,],n.iter=2000,ntree=50,burn.in=500,
                 thinning=5,n.group=n.group)
    
    ## Create gg survival matrix
    ggamma.surv = matrix(NA, nrow(X[fold,]),length(unique.times))
    
    for (j in 1:n.group){
      for (i in which(X$group[fold] == j)){
        ggamma.surv[i,] = 1 - pgengamma(unique.times,mu=median(ggITE$beta.test[,i]),
                                        sigma=exp(median(ggITE$logsig.mat[,j])),
                                        Q=median(ggITE$loglam.mat[,j])) }}
    
    ## C indeces
    cindex.ggamma = cindex(object=ggamma.surv,formula=Surv(Survival,Status)~., data=data[fold,],eval.times=unique.times)
    brier.ggamma = pec(object=ggamma.surv,formula=Surv(Survival,Status)~1,data=Y[fold,],start=unique.times[1],maxtime=unique.times[length(unique.times)],
                       exact=FALSE, times=unique.times,reference=FALSE,splitMethod = "none")
    
    return(data.frame('cindex.ggamma'=cindex.ggamma$AppCindex,'brier.gammma'=brier.ggamma$AppErr,'unique.times'=unique.times))
  }#XVAL.FOLD
  
  results = lapply(folds, xval.fold)
  return(do.call("rbind", results))
  
}#XVAL

perf.est = function(data,n.group,k){
  xval.out = xval(data,n.group=n.group,k=k)
  return(data.frame(gg.cindex=tapply(xval.out[,1],as.factor(xval.out[,3]),mean),
                    sd.gg.cindex=tapply(xval.out[,1],as.factor(xval.out[,3]),sd),
                    gg.brier=tapply(xval.out[,2],as.factor(xval.out[,3]),mean),
                    sd.gg.cindex=tapply(xval.out[,2],as.factor(xval.out[,3]),sd)))
}#PERFORMANCE ESTIMATE

perf.table = function(perf.est,quant=c(0.1,0.25,0.5,0.75,0.9)){
  time = unique(sort(Y$Survival))
  quantiles <- as.numeric(sapply(quantile(time, probs=quant), function(y)which.min(abs(time - y))))
  perf.table = perf.est[quantiles,]
  rownames(perf.table) <- quant
  t(perf.table)
}#PERFORMANCE TABLE