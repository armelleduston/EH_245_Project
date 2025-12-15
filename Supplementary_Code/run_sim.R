### run_sim


### load other libraries 
library(GIGrvg)
library(refund)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(splines)
library(MASS)
library(bkmr) ## BKMR
library(NLinteraction) ## hierarchical method
library(spikeSlabGAM) ## ssGAM
library(bsmim2)
library(bayesSurv)
library(CKMR) ## install_github("glenmcgee/CKMR")
source("MixSelect-master/MixSelect.R") ## available at https://github.com/fedfer/MixSelect
#### MixSelect requires 2 archived packages not on CRAN
###### first install smoothSurv_2.5.tar.gz, then install bayesSurv_3.6.tar.gz

## set to TRUE to run locally ## FALSE is on cluster
runLOCAL=FALSE 


### set up cluster 
if(runLOCAL==TRUE){
  path <- "" ## path for results
  index=0
  interaction=0
  iter_no=1
  suffix <- paste0("_index",index,"_interac",interaction,"_iter",iter_no) ## append output file names with the iteration number, to be combined later
} else{
  path <- "Results/" ## path for results
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  index <- as.integer(args[1])        ## 0 = bkmr, 1 = index structure
  interaction <- as.integer(args[2])  ## 0 = main effects only; 1 = interaction; 2 = supplementary
  iter_no <- as.integer(args[3])      ## get iteration number
  suffix <- paste0("_index",index,"_interac",interaction,"_iter",iter_no) ## append output file names with the iteration number, to be combined later
}

## set other params
nn <- 200 ##  sample size
RR <- 30000
burnpct <- 0.5
thinnum <- 10
sel <- seq(burnpct*RR+1,RR,by=thinnum)  ## for bkmr
if(index==FALSE){
  ncov <- 10
  nindex <- 10
}else{
  ncov <- 16
  nindex <- 10
}
NLdf <- 6 # d.o.f for splines in NLinter method
xcor <- 0.25 # exposure correlation
p <- ncov

#################################
## Generate data

gendat <- function(n,p,interac=FALSE,ind=FALSE){
  
  ## covariates
  Xmat <- scale(matrix(runif(p*n,0,1),ncol=p,nrow=n))
  colnames(Xmat) <- paste0("x",1:p)
  
  z <- matrix(rnorm(n*2,0,1),ncol=2) ## matrix of 2 columns corresponding to 2 covariates (no column of ones since we dont need an intercept)
  colnames(z) <- paste0("z",1:2)
  
  if(ind==FALSE){ ## non-index (standard bkmr case)
    
    Xlist <- list() ## for fitting
    for(jj in 1:p){
      Xlist[[jj]] <- as.matrix(Xmat[,jj],ncol=1)
    }
    Xmattheta <- Xmat ## for generating data 

  }else{ ## indexwise
    Xlist <- list() ## for fitting
    theta1 <- c(3,2,1,0); theta1 <- theta1/sqrt(sum(theta1^2))
    theta2 <- c(1,1,1,1); theta2 <- theta2/sqrt(sum(theta2^2))
    Xlist[[1]] <- as.matrix(Xmat[,1:4],ncol=4)
    Xlist[[2]] <- as.matrix(Xmat[,5:8],ncol=4)
    for(jj in 9:p){ ## all others are individual
      Xlist[[jj-6]] <- as.matrix(Xmat[,jj],ncol=1)
    }
    Xmattheta <- cbind(Xlist[[1]]%*%theta1,
                       Xlist[[2]]%*%theta2,
                       Xmat[,9:p]) ## for generating data 
    
  }
  
  
  ## surface
  h1 <- 2*cos(2*Xmattheta[,1])
  h2 <- Xmattheta[,2]
  h3 <- 4*dt(2*Xmattheta[,3],df=10) 
  h4 <- sin(2*Xmattheta[,4])  
  h5 <- Xmattheta[,5]^2 
  h6 <- -Xmattheta[,6] 

  h <- h1+h2+h3+h4+h5+h6
  if(interac==TRUE){
    h <- h+0.5*h1*h5
  }
  if(interac==2){ ## supplementary sim without "main" effect
    h <- h2+h3+h4+h6 ## no main effects
    h <- h+h1*h5
  }

  ## generate outcomes
  y <- h+z%*%c(0.5,0)+rnorm(n,0,sqrt(1))

  ## for two-way interaction
  if(ind==FALSE){
    XInter <- matrix(0,nrow=60,ncol=ncol(Xmat))
    XInter[1:30,1] <- rep(seq(quantile(Xmat[,1],0.05),quantile(Xmat[,1],0.95),length=10),3)
    XInter[1:30,5] <- rep(quantile(Xmat[,5],c(0.1,0.5,0.9)),each=10)
    XInter[31:60,1] <- rep(quantile(Xmat[,1],c(0.1,0.5,0.9)),each=10)
    XInter[31:60,5] <- rep(seq(quantile(Xmat[,5],0.05),quantile(Xmat[,5],0.95),length=10),3)
    XmatthetaInter <- XInter
  }else{
    XInter <- lapply(Xlist,function(x) return(matrix(0,nrow=60,ncol=ncol(x))))
    for(jj in 1:4){
      XInter[[1]][,jj] <- c(rep(seq(quantile(Xlist[[1]][,jj],0.05),quantile(Xlist[[1]][,jj],0.95),length=10),3),
                            rep(quantile(Xlist[[1]][,jj],c(0.1,0.5,0.9)),each=10))
    }
    XInter[[5]] <- matrix(c(rep(quantile(Xmattheta[,5],c(0.1,0.5,0.9)),each=10),
                          rep(seq(quantile(Xmattheta[,5],0.05),quantile(Xmattheta[,5],0.95),length=10),3))
                          ,ncol=1,nrow=60)
    XmatthetaInter <- cbind(XInter[[1]]%*%theta1,
                            XInter[[2]]%*%theta2,
                            Reduce("cbind",XInter[3:nindex]))  
  }

  
  return(list(y=y,
              Xmat=data.frame(Xmat),
              Xlist=Xlist,
              XInter=XInter,
              XmatthetaInter=XmatthetaInter,
              z=z,
              h=h)) ## for testing
}

set.seed(1000+iter_no)
dat <- gendat(nn,p=ncov,interac=interaction,ind=index) ## generate data
Xmat <- dat$Xmat ## matrix of all mixture components
Xlist <- dat$Xlist ## list of mixture components
XInter <- dat$XInter ## mat/list of mixture components for interactions
XmatthetaInter <- dat$XmatthetaInter ## matrix of mixture components for interactions
covariates <- dat$z ## just a matrix of exposures
response <- dat$y ## just a vector of outcomes


#################################
## Model fitting

## CBMIM
cbmim_mod <- cbmim(y=response,x=Xlist,z=covariates, ## y is outcome; ## x is list of exposure index matrices; ## z is matrix of covariates (need at least 1)
             niter=RR, ## number of iterations 
             nburn=burnpct*RR, ## burn-in fraction
             nthin=thinnum, ## thinning number
             nchains=1, ## number of chains ## not yet implemented
             ## prior hyperparameters
             prior_pi_additive=c(1,1), ## prior for pi in spike & slab on additive component
             prior_tau2_additive=c(3,50), ## shape and rate for inverse gamma on spline penalty terms
             prior_pi_nonadditive=c(1,1), ## prior for pi in spike & slab on non-additive component
             prior_rho_nonadditive=c(5,5), ## shape and rate for gamma
             prior_tau2_nonadditive=c(3,20),## for the inverse gamma that matches moments with bkmr, which useswith mean 10 sd 10. #c(0.001,0.001), ## shape and rate for inverse gamma on nonadditive GP component
             prior_sigma2=c(1,1), ## shape and rate for inverse gamma on \sigma^2
             ## MH tuning
             stepsize_tau2_nonadditive=5, ##jumpsize/sd for gamma proposal on tau2_nonadditive
             stepsize_rho_nonadditive=0.5, ##jumpsize/sd for gamma proposal on rho
             ## approximating large inverse via Sherman-Morrison-Woodbury
             invmethod=c("exact"), ## "exact" for full matrix inversion, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
             rank=100, ## rank of approximation for lz/rp approximations of PKP (SMW) or I+PKP (DIRECT)
             kernelfun="gaussian", ## choice of kernel function
             approxProj = FALSE, ## dont update P 
             draw_h=FALSE,
             centering=TRUE, ## should covariates be centered to have mean=0
             scaling=TRUE, 
             hierarchical=TRUE,
             oversample=TRUE,
             polar=FALSE)

cbmim_pred <- pred_surface(cbmim_mod,Xnew=Xlist)
cbmim_pred <- summarize_pred(cbmim_pred,assoc=FALSE,centered=FALSE)

cbmim_bias <- mean(cbmim_pred$mean-dat$h)
cbmim_mse <- sqrt(mean((cbmim_pred$mean-dat$h)^2))
cbmim_width <- mean(cbmim_pred$uci-cbmim_pred$lci)
cbmim_cvg <- mean((cbmim_pred$lci<=dat$h) & (dat$h<=cbmim_pred$uci))

cbmim_PIPs <- get_PIPs(cbmim_mod)
cbmim_mainPIPs <- cbmim_PIPs$mainPIPs
cbmim_interPIPs <- lower_tri(t(cbmim_PIPs$twowayPIPs)) ## in order of (1,2); (1,3); (1,4);.... (2,3); (2,4); (3,4)


cbmim_results <- c(cbmim_bias,
                   cbmim_mse,
                   cbmim_width,
                   cbmim_cvg,
                   cbmim_mainPIPs,
                   cbmim_interPIPs)


## get two-way interactions for plotting
if(interaction==1 & index==0){
  cbmim_pred_inter_raw <- pred_surface_indexwise_interac(cbmim_mod,
                                                         includeInt=TRUE,
                                                         gridlen = 10,
                                                         whichids=c(1,5),restrict = FALSE)
  cbmim_pred_inter <- summarize_pred(cbmim_pred_inter_raw,assoc=F,centered=F)
  
  ## get postmeans
  cbmim_pred_intermeans <- Reduce(c,lapply(1:2,function(jj){
    Reduce(c,lapply(1:3,function(ll){
      cbmim_pred_inter[[jj]][[1]][[ll]]$mean
    }))
  }))
  
  cbmim_inter <- data.frame(iter_no=iter_no,
                            var1=c(rep("X1",30),rep("X5",30)),
                            qntl=c(rep(c(0.1,0.5,0.9),each=10),rep(c(0.1,0.5,0.9),each=10)),
                            X1=XInter[,1],
                            X5=XInter[,5],
                            pred=cbmim_pred_intermeans)
  
  write.csv(cbmim_inter,file=paste0(path,"inter_sim",suffix,".csv"),row.names=F)
  
  
}else if(interaction==1){
  cbmim_pred_inter_raw <- pred_surface(cbmim_mod,Xnew=XInter,includeInt=TRUE)
  cbmim_pred_intermeans <- summarize_pred(cbmim_pred_inter_raw,assoc=F,centered=F)$mean
  
  cbmim_inter <- data.frame(iter_no=iter_no,
                            var1=c(rep("X1",30),rep("X5",30)),
                            qntl=c(rep(c(0.1,0.5,0.9),each=10),rep(c(0.1,0.5,0.9),each=10)),
                            X1=XmatthetaInter[,1],
                            X5=XmatthetaInter[,5],
                            pred=cbmim_pred_intermeans)
  
  write.csv(cbmim_inter,file=paste0(path,"inter_sim",suffix,".csv"),row.names=F)
  
}



## CKMR with nonadaptive projection
nonadapt_mod <- cbmim(y=response,x=Xlist,z=covariates, ## y is outcome; ## x is list of exposure index matrices; ## z is matrix of covariates (need at least 1)
                      niter=RR, ## number of iterations 
                      nburn=burnpct*RR, ## burn-in fraction
                      nthin=thinnum, ## thinning number
                      nchains=1, ## number of chains ## not yet implemented
                      ## prior hyperparameters
                      prior_pi_additive=c(1,1), ## prior for pi in spike & slab on additive component
                      prior_tau2_additive=c(3,50), ## shape and rate for inverse gamma on spline penalty terms
                      prior_pi_nonadditive=c(1,1), ## prior for pi in spike & slab on non-additive component
                      prior_rho_nonadditive=c(5,5), ## shape and rate for gamma
                      prior_tau2_nonadditive=c(3,20),## for the inverse gamma that matches moments with bkmr, which useswith mean 10 sd 10. #c(0.001,0.001), ## shape and rate for inverse gamma on nonadditive GP component
                      prior_sigma2=c(1,1), ## shape and rate for inverse gamma on \sigma^2
                      ## MH tuning
                      stepsize_tau2_nonadditive=5, ##jumpsize/sd for gamma proposal on tau2_nonadditive
                      stepsize_rho_nonadditive=0.5, ##jumpsize/sd for gamma proposal on rho
                      ## approximating large inverse via Sherman-Morrison-Woodbury
                      invmethod=c("exact"), ## "exact" for full matrix inversion, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
                      rank=100, ## rank of approximation for lz/rp approximations of PKP (SMW) or I+PKP (DIRECT)
                      kernelfun="gaussian", ## choice of kernel function
                      approxProj = TRUE, ## dont update P 
                      draw_h=FALSE,
                      centering=TRUE, ## should covariates be centered to have mean=0
                      scaling=TRUE, 
                      hierarchical=TRUE,
                      polar=FALSE)

nonadapt_pred <- pred_surface(nonadapt_mod,Xnew=Xlist)
nonadapt_pred <- summarize_pred(nonadapt_pred,assoc=FALSE,centered=FALSE)

nonadapt_bias <- mean(nonadapt_pred$mean-dat$h)
nonadapt_mse <- sqrt(mean((nonadapt_pred$mean-dat$h)^2))
nonadapt_width <- mean(nonadapt_pred$uci-nonadapt_pred$lci)
nonadapt_cvg <- mean((nonadapt_pred$lci<=dat$h) & (dat$h<=nonadapt_pred$uci))

nonadapt_PIPs <- get_PIPs(nonadapt_mod)
nonadapt_mainPIPs <- nonadapt_PIPs$mainPIPs
nonadapt_interPIPs <- lower_tri(t(nonadapt_PIPs$twowayPIPs)) ## in order of (1,2); (1,3); (1,4);.... (2,3); (2,4); (3,4)

nonadapt_results <- c(nonadapt_bias,
                      nonadapt_mse,
                      nonadapt_width,
                      nonadapt_cvg,
                      nonadapt_mainPIPs,
                      nonadapt_interPIPs)

if(index==FALSE){ ## run BKMR, NLinter, ssGAM , nonadaptProj,
  
  ### BKMR via kmbayes
  bkmr_mod <-  kmbayes(y=response, Z=Xmat, X=covariates, 
                       iter=RR, verbose=FALSE, varsel=TRUE)
  
  bkmr_pred <- ComputePostmeanHnew(bkmr_mod, Znew = Xmat, sel = sel, method = "exact") ## using approx method as exact method errors out due to memory
  bkmr_pred <- data.frame(mean=bkmr_pred$postmean,
                          lci=bkmr_pred$postmean-1.96*sqrt(diag(bkmr_pred$postvar)),
                          uci=bkmr_pred$postmean+1.96*sqrt(diag(bkmr_pred$postvar)) )
  
  bkmr_bias <- mean(bkmr_pred$mean-dat$h)
  bkmr_mse <- sqrt(mean((bkmr_pred$mean-dat$h)^2))
  bkmr_width <- mean(bkmr_pred$uci-bkmr_pred$lci)
  bkmr_cvg <- mean((bkmr_pred$lci<=dat$h) & (dat$h<=bkmr_pred$uci))
  
  bkmr_mainPIPs <- ExtractPIPs(bkmr_mod,sel=sel)$PIP
  bkmr_interPIPs <- rep(NA,ncov*(ncov-1)/2)
  
  bkmr_results <- c(bkmr_bias,
                    bkmr_mse,
                    bkmr_width,
                    bkmr_cvg,
                    bkmr_mainPIPs,
                    bkmr_interPIPs)
  
  
  ### Non-linear interactions via Antonelli et al
  NL_mod <- NLint(Y=response,X=Xmat,C=covariates,
                  nIter = RR/2,nBurn=burnpct*RR/2,thin=thinnum,
                  nChains=2,ns=NLdf) ## errors with 1 chain
  
  NL_pred <- NLpredict(NL_mod,X=Xmat,Xnew=Xmat,Cnew=0*covariates) ## 0*covariates since we only want h
  NL_pred <- data.frame(Reduce(cbind,NL_pred[1:3]))
  colnames(NL_pred) <- c("mean","lci","uci")
  
  NL_bias <- mean(NL_pred$mean-dat$h)
  NL_mse <- sqrt(mean((NL_pred$mean-dat$h)^2))
  NL_width <- mean(NL_pred$uci-NL_pred$lci)
  NL_cvg <- mean((NL_pred$lci<=dat$h) & (dat$h<=NL_pred$uci))
  
  NL_mainPIPs <- NL_mod$MainPIP
  NL_interPIPs <- lower_tri(t(NL_mod$InteractionPIP)) ## in order of (1,2); (1,3); (1,4); (2,3); (2,4); (3,4)
  
  NL_results <- c(NL_bias,
                  NL_mse,
                  NL_width,
                  NL_cvg,
                  NL_mainPIPs,
                  NL_interPIPs)
  
  
  ### ssGAM
  dd <- data.frame(response,Xmat,covariates)
  newdd <- data.frame(response,Xmat,0*covariates)
  
  GAMfrmla <- as.formula(paste0("response~",
                                paste(paste0("(x",1:ncov,")+"),collapse=""),
                                "lin(z1)+lin(z2)"))
  ssGAM_mod <- spikeSlabGAM(GAMfrmla, data = dd,mcmc=list(nChains=1))
  
  
  
  ## code for estimates of h
  ssGAM_pred <- (ssGAM_mod$X[,1:(ncol(ssGAM_mod$X)-2)]%*%t(ssGAM_mod$samples$beta[[1]][,1:(ncol(ssGAM_mod$X)-2)])) ## -2 excludes the covariates
  ssGAM_pred <- as.data.frame(t(apply(ssGAM_pred,1,function(x) return(c(mean(x),
                                                                        quantile(x,c(0.025,0.975)))))))
  colnames(ssGAM_pred) <- c("mean","lci","uci")
  
  ssGAM_bias <- mean(ssGAM_pred$mean-dat$h)
  ssGAM_mse <- sqrt(mean((ssGAM_pred$mean-dat$h)^2))
  ssGAM_width <- mean(ssGAM_pred$uci-ssGAM_pred$lci)
  ssGAM_cvg <- mean((ssGAM_pred$lci<=dat$h) & (dat$h<=ssGAM_pred$uci))
  
  ## need to combine the linear and smooth components for each variable
  ssGAM_mainPIPs <- sapply(1:ncov,function(jj){
    mean(apply(as.matrix(ssGAM_mod$samples$gamma[,2*(jj-1)+(1:2)]),1,max))
  })
  ssGAM_interPIPs <- rep(NA,ncov*(ncov-1)/2)
  
  ssGAM_results <- c(ssGAM_bias,
                     ssGAM_mse,
                     ssGAM_width,
                     ssGAM_cvg,
                     ssGAM_mainPIPs,
                     ssGAM_interPIPs)
  
  ## mix select
  #### requires 2 archived packages not on CRAN
  #### first install smoothSurv_2.5.tar.gz
  #### then install bayesSurv_3.6.tar.gz
  try( ## error catching because it often doesnt work
    mixselect_mod <- MixSelect(y=response,X=as.matrix(Xmat),Z=covariates,
                               nrun = RR, 
                               burn = burnpct*RR,
                               thin = thinnum)
  )
  if(!exists("mixselect_mod")){ ## if it didnt work, just continue with other methods
    mixselect_results <- rep(NA,length(cbmim_results))
    print("mixselect FAILED. ")
  }else{
    mixselect_pred <- mixselect_mod$y_hat - (mixselect_mod$beta_z%*%t(covariates)) ## subtracting off the Z\beta from yhat (that got added in automatically)
    mixselect_pred <- data.frame(mean=apply(mixselect_pred,2,mean),
                                 lci=apply(mixselect_pred,2,function(x) quantile(x,0.025)),
                                 uci=apply(mixselect_pred,2,function(x) quantile(x,0.975)))
    
    mixselect_bias <- mean(mixselect_pred$mean-dat$h)
    mixselect_mse <- sqrt(mean((mixselect_pred$mean-dat$h)^2))
    mixselect_width <- mean(mixselect_pred$uci-mixselect_pred$lci)
    mixselect_cvg <- mean((mixselect_pred$lci<=dat$h) & (dat$h<=mixselect_pred$uci))
    
    ## need to add interactions AND non linear effects to be true PIPs(otherwise it only picks up on linear effects)
    compute_mixPIPs <- function(mod){
      
      indmat <- matrix(NA,ncol=p,nrow=p)
      indmat[lower.tri(indmat)] <- 1:(p*(p-1)/2)
      return(sapply(1:p,function(j) {
              whichcols <- c(indmat[j,0:(j-1)],indmat[min(p,(j+1)):p,j])
              whichcols <- whichcols[!is.na(whichcols)] # exclude diagonals
              mean(apply(cbind(mod$gamma_int[,whichcols],mod$gamma_beta[,j],mod$gamma_int[,j]),2,max))}
              ) )
      
    }
    mixselect_mainPIPs <- compute_mixPIPs(mixselect_mod)
    mixselect_interPIPs <- apply(mixselect_mod$gamma_int,2,mean) ## in order of (1,2); (1,3); (1,4);.... (2,3); (2,4); (3,4)
    
    mixselect_results <- c(mixselect_bias,
                           mixselect_mse,
                           mixselect_width,
                           mixselect_cvg,
                           mixselect_mainPIPs,
                           mixselect_interPIPs)
  }
  
  
  ## save results
  results <- data.frame(rbind(cbmim_results,
                              bkmr_results,
                              NL_results,
                              ssGAM_results,
                              nonadapt_results,
                              mixselect_results))
  colnames(results) <- c("bias","mse","width","cvg",paste0("mainPIPs",1:ncov),paste0("interPIPs",1:(ncov*(ncov-1)/2)))
  
  results$method <- c("CBMIM","BKMR","NLInter","ssGAM","NonAdaptive","mixselect")
  results$iter_no <- iter_no
  
  
}else{ ## indexwise: run standard BMIM
  
  bmim_mod <- bsmim2(y=response,x=Xlist,z=as.matrix(covariates),
                     niter=RR,nburn=burnpct*RR,nthin=thinnum,
                     centering=F,scaling=F,
                     prior_sigma=c(0.001,0.001),
                     prior_lambda_shaperate=c(1,0.1),
                     gaussian=TRUE,
                     spike_slab=TRUE,
                     gauss_prior=TRUE,
                     prior_theta_slab_sd=0.5,
                     stepsize_theta=0.1,#0.35
                     draw_h=FALSE,
                     num_theta_steps = 10)

  bmim_pred <- predict_hnew_X2(bmim_mod,newX=Xlist)
  bmim_pred <- bmim_pred$fits[,c("mean","lower","upper")]
  colnames(bmim_pred) <- c("mean","lci","uci")
  
  bmim_bias <- mean(bmim_pred$mean-dat$h)
  bmim_mse <- sqrt(mean((bmim_pred$mean-dat$h)^2))
  bmim_width <- mean(bmim_pred$uci-bmim_pred$lci)
  bmim_cvg <- mean((bmim_pred$lci<=dat$h) & (dat$h<=bmim_pred$uci))
  
  ## need to combine the linear and smooth components for each variable
  bmim_mainPIPs <- apply(bmim_mod$rho!=0,2,mean) 
  bmim_interPIPs <- rep(NA,nindex*(nindex-1)/2)
  
  bmim_results <- c(bmim_bias,
                     bmim_mse,
                     bmim_width,
                     bmim_cvg,
                     bmim_mainPIPs,
                     bmim_interPIPs)
  
  
  
  
  ## save results
  results <- data.frame(rbind(cbmim_results,
                              bmim_results,
                              nonadapt_results))
  colnames(results) <- c("bias","mse","width","cvg",paste0("mainPIPs",1:nindex),paste0("interPIPs",1:(nindex*(nindex-1)/2)))
  
  results$method <- c("CBMIM","BMIM","NonAdaptive")
  results$iter_no <- iter_no
  
}


## save results
write.csv(results,file=paste0(path,"results_sim",suffix,".csv"),row.names=F)


