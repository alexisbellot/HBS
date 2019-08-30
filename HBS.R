### Nonparametric hierchical Bayes model for survival

HBS <-function(X.train,Y.train,X.test,n.iter,sigma.prior.bart=10,ntree, n.group,
                           burn.in=1,thinning=1){
  
  ### import libraries
  library(BayesTree)
  library(flexsurv)
  
  ### initialize space allocation
  delta                 <- Y.train$Status
  time                  <- as.vector(Y.train$Survival) #event time
  n.samples             <- nrow(X.train)
  n.samples.test        <- nrow(X.test)
  # prior means for beta
  sigma.error           <- matrix(NA,n.iter,1)
  g.of.x.train          <- matrix(NA,n.iter,n.samples)
  g.of.x.test           <- matrix(NA,n.iter,n.samples.test)
  # beta matrices
  beta.mat              <- matrix(NA,(n.iter+1),n.samples)
  beta.mat.test         <- matrix(NA,(n.iter+1),n.samples.test)
  beta.jump.mat         <- matrix(NA,(n.iter+1),n.samples) #matrix for storing beta proposals 
  beta.jump.mat.test    <- matrix(NA,(n.iter+1),n.samples.test)
  sigma.jump.beta.mat   <- matrix(NA,(n.iter+1),n.samples)
  # logsig matrices
  logsig.mat            <- matrix(NA,(n.iter+1),n.group)
  logsig.jump.mat       <- matrix(NA,(n.iter+1),n.group) 
  sigma.jump.logsig.mat <- matrix(NA,(n.iter+1),n.group)
  # lambda matrices
  loglam.mat            <- matrix(NA,(n.iter+1),n.group)
  loglam.jump.mat       <- matrix(NA,(n.iter+1),n.group)
  sigma.jump.loglam.mat <- matrix(NA,(n.iter+1),n.group)
  
  ### initialize variables
  initial               <- flexsurvreg(Surv(Y.train$Survival, Y.train$Status) ~ 1, dist = "gengamma")
  beta.mat[1,]          <- rnorm(n.samples,0,1) + initial$coeff[[1]]
  logsig.mat[1,]        <- rnorm(n.group,0,1) + initial$coeff[[2]] 
  loglam.mat[1,]        <- rnorm(n.group,0,1) + initial$coeff[[3]] 
  sigma.jump.beta       <- rep(2.5,n.samples)
  sigma.jump.logsig     <- rep(1,n.group)
  sigma.jump.loglam     <- rep(1,n.group)

  group                 <- X.train$group
  group.mat             <- matrix(nrow=nrow(X.train),ncol=n.group)
  for(i in 1:n.group){group.mat[,i] <- as.numeric(X.train$group == i)}
  
  ndpost                <- 200 #BART posterior draws after burn in 
  nskip                 <- 100 #BART burn in
  
  ####################################################################################################
  ## functions on for MCMC
  #log of posterior (likelihood*prior) of beta=location to be used in MH ratio
  log.post.beta<-function(beta,loglam,logsig,delta,time,sigma.error,g.of.x,group.mat){ #group.mat and logsig must be vectors
    summation = 0
      for(i in 1:n.group){
         temp      <- sum(group.mat[i]*(delta*log(0.001+dgengamma(time,mu=beta,sigma=exp(logsig[i]),Q=loglam[i])) + 
                         (1-delta)*log(1.001-pgengamma(time,mu=beta,sigma=exp(logsig[i]),Q=loglam[i]))))
         summation <- summation + temp}
    summation + sum(dnorm(beta,g.of.x,sigma.error,log=TRUE))}
  
  #MCMC update for beta=location
  beta.update <- function(sigma.jump.beta){
    beta.star     <- rnorm(1,beta.mat[wei,ii],sigma.jump.beta) 
    log.post.old  <- log.post.beta(beta.mat[wei,ii],loglam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                                sigma.error[wei],g.of.x.train[wei,ii],group.mat[ii,])
    log.post.star <- log.post.beta(beta.star,loglam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                                 sigma.error[wei],g.of.x.train[wei,ii],group.mat[ii,])
    r             <- exp(log.post.star-log.post.old) 
    beta          <- ifelse(runif(1)<r,beta.star,beta.mat[wei,ii])
    
    if(is.na(beta)==TRUE||is.infinite(beta)==TRUE||is.nan(beta)==TRUE) {beta<-1;print('BETA NA')}
    p.jump<-min(r,1) # probability of jumping
    
    list(beta=beta,p.jump=p.jump)
  }
  
  #log of posterior (likelihood*prior) of loglam = log(lambda=shape) to be used in MH ratio
  log.post.lambda<-function(beta,loglam,logsig,delta,time,group_ind){
    summation <- 0
    for (ii in group_ind){
      temp      <- delta[ii]*log(dgengamma(time[ii],mu=beta[ii],sigma=exp(logsig),Q=loglam)+0.001) + 
                   (1-delta[ii])*log(1.001-pgengamma(time[ii],mu=beta[ii],sigma=exp(logsig),Q=loglam))
      summation <- summation+temp
      
    }
    sum(summation) + sum(dnorm(loglam,mean=0,sd=1,log=TRUE))
  }
  
  #MCMC update for log(lambda=shape)
  loglam.update <- function(sigma.jump.loglam,group_ind){
    loglam.star   <- rnorm(1,loglam.mat[wei,kk],sigma.jump.loglam)
    log.post.old  <- log.post.lambda(beta.mat[wei+1,],loglam.mat[wei,kk],logsig.mat[wei+1,kk],delta,time,group_ind)
    log.post.star <- log.post.lambda(beta.mat[wei+1,],loglam.star,logsig.mat[wei+1,kk],delta,time,group_ind)
    r             <- exp(log.post.star-log.post.old)
    loglam        <- ifelse(runif(1)<r,loglam.star,loglam.mat[wei,kk])
    if(is.na(loglam)==TRUE||is.nan(loglam)==TRUE||is.infinite(loglam)==TRUE) {loglam<-loglam.mat[wei,kk];print('loglam NA')}
    p.jump.loglam <- min(r,1)
    
    list(loglam=loglam,p.jump.loglam=p.jump.loglam)
  }
   
  #log of posterior (likelihood*prior) of sigma to be used in MH ratio
  log.post.sigma <- function(beta,loglam,logsig,delta,time,group_ind){
    summation <- 0
    for (ii in group_ind){
      temp      <- delta[ii]*log(dgengamma(time[ii],mu=beta[ii],sigma=exp(logsig),Q=loglam)+0.001) + 
                  (1-delta[ii])*log(1.001-pgengamma(time[ii],mu=beta[ii],sigma=exp(logsig),Q=loglam))
      summation <- summation+temp
      
    }
    sum(summation) + sum(dnorm(logsig,mean=0,sd=1,log=TRUE))
  }
  
  
  
  #MCMC update for log(sigma=scale)
  logsig.update <- function(sigma.jump.logsig,group_ind ){
    logsig.star   <- rnorm(1,logsig.mat[wei,jj],sigma.jump.logsig)
    log.post.old  <- log.post.sigma(beta.mat[wei+1,],loglam.mat[wei,jj],logsig.mat[wei,jj],delta,time,group_ind)
    log.post.star <- log.post.sigma(beta.mat[wei+1,],loglam.mat[wei,jj],logsig.star,delta,time,group_ind)
    r             <- exp(log.post.star-log.post.old)
    logsig        <- ifelse(runif(1)<r,logsig.star,logsig.mat[wei,jj])
    
    if(is.na(logsig)==TRUE) {logsig<-logsig.mat[wei,jj];print('logsig NA')}
    p.jump        <- min(r,1)
    
    list(logsig=logsig,p.jump=p.jump)
  }
  
  ##############################################################################################
  #FULL MCMC procedure
  for (wei in 1:n.iter){
    if (wei %% 20 == 1){
    bart1               <- bart(X.train,beta.mat[wei,],x.test=X.test,sigest=sigma.prior.bart,ntree=ntree,
                                ndpost=ndpost,nskip=nskip,verbose=FALSE) 
    sigma.error[wei]    <- sqrt(bart1$sigma[1]) 
    g.of.x.train[wei,]  <- bart1$yhat.train.mean
    beta.mat.test[wei,] <- bart1$yhat.test.mean
    print("Done BART")
    }
    else{
      sigma.error[wei]    <- sigma.error[wei-1]
      g.of.x.train[wei,]  <- g.of.x.train[wei-1,]
      beta.mat.test[wei,] <- beta.mat.test[wei-1,]
    }
    
    
    for (ii in 1:n.samples){
      if (wei %% 50 == 0) {sigma.jump.beta[ii] <- adaptive_jump(sigma.jump.beta[ii],beta.jump.mat[,ii],wei)}
      temp                          <- beta.update(sigma.jump.beta[ii])
      beta.mat[wei+1,ii]            <- temp$beta
      beta.jump.mat[wei+1,ii]       <- temp$p.jump
      sigma.jump.beta.mat[wei+1,ii] <- sigma.jump.beta[ii]
    }
    print('DONE BETA UPDATE')
    
    for (jj in 1:n.group){# num treatment is 
      group_ind                       <- which(X.train$group == jj)
      if (wei %% 50 == 0) {sigma.jump.logsig[jj] <- adaptive_jump(sigma.jump.logsig[jj],logsig.jump.mat[,jj],wei)}
      temp                            <- logsig.update(sigma.jump.logsig[jj],group_ind)
      logsig.mat[wei+1,jj]            <- temp$logsig
      logsig.jump.mat[wei+1,jj]       <- temp$p.jump
      sigma.jump.logsig.mat[wei+1,jj] <- sigma.jump.logsig[jj]
    }

    print('DONE LOGSIG UPDATE')
    
    # lambda update
    for (kk in 1:n.group){# num treatment is 
      group_ind                       <- which(X.train$group == kk)
      if (wei %% 50 == 0) {sigma.jump.loglam[kk] <- adaptive_jump(sigma.jump.loglam[kk],loglam.jump.mat[,kk],wei)}
      temp                            <- loglam.update(sigma.jump.loglam[kk],group_ind)
      loglam.mat[wei+1,kk]            <- temp$loglam
      loglam.jump.mat[wei+1,kk]       <- temp$p.jump
      sigma.jump.loglam.mat[wei+1,kk] <- sigma.jump.loglam[kk]
    }
    
    print(sprintf("DRAW_GGAMMA.TREE_%i",wei))
  }
  
  
  #post-processing
  draws                <- seq(burn.in,n.iter,thinning)
  beta.mat             <- as.data.frame(beta.mat)
  names(beta.mat)      <- row.names(X.train)
  beta.mat.test        <- as.data.frame(beta.mat.test)
  names(beta.mat.test) <- row.names(X.test)
  
  
  list(beta.train=beta.mat[draws,],beta.jump.mat=beta.jump.mat[draws,],beta.test=beta.mat.test[draws,],
       logsig.mat=logsig.mat[draws,],logsig.jump.mat=logsig.jump.mat[draws,],loglam.mat=loglam.mat[draws,],
       loglam.jump.mat=loglam.jump.mat[draws,],sigma.jump.beta.mat=sigma.jump.beta.mat[draws,],
       sigma.jump.logsig.mat=sigma.jump.logsig.mat[draws,],sigma.jump.loglam.mat=sigma.jump.loglam.mat[draws,],
       sigma.error=sigma.error[draws])
  
}
