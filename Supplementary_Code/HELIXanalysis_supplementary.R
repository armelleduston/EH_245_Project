## set up
options(echo=TRUE)
options(stringsAsFactors = FALSE)
exp_set="postnatal" 

library(tidyverse)
# devtools::install_github("glenmcgee/bsmim2")
library(bsmim2)
library(MASS)
library(parallel)
library(CKMR)
load("exposome_v2.RData")
load("covariates.rda") ## from Make_Covariates.R
load("index_exposure_list.rda") ## from Make_Exposures.R

################
## exposures ###
################
exposure_list <- index_exposure_list[grepl("Post",names(index_exposure_list))]
exposure_list <- exposure_list[-c(2,14)]
exposure_list <- lapply(exposure_list,scale)

#### SUPPLEMENTARY ANALYSIS
##### Include metals separately:
for(cc in 1:ncol(exposure_list$Metals_Postnatal)){
  exposure_list[[length(exposure_list)+1]] <- matrix(exposure_list$Metals_Postnatal[,cc],ncol=1)
  names(exposure_list)[length(exposure_list)] <- colnames(exposure_list$Metals_Postnatal)[cc]
}
exposure_list <- exposure_list[-4] # remove metals index





################
## covariates ##
################
z <- cbind(z_base_preg,z_base_post)
# remove covariates that are already adjusted for in the outcome or related to the outcome
z <- z[,-which(colnames(z)%in%c("hs_child_age_None","e3_sex_Nonemale","hs_c_weight_None","hs_c_height_None"))]






#####################
## Full CBMIM Fit  ##
#####################
nit <- 10000
thin = 5
burnpct <- 0.5
set.seed(1234)
  
  cbmim_mod <- cbmim(y=scale(phenotype$hs_zbmi_who),x=exposure_list,z=z, ## y is outcome; ## x is list of exposure index matrices; ## z is matrix of covariates (need at least 1)
                     niter=nit, ## number of iterations 
                     nburn=burnpct*nit,## burn-in fraction
                     nthin=thin, ## thinning number
                     nchains=4, ## number of chains 
                     ncores=4,
                     ## prior hyperparameters
                     prior_theta_kappa=5, 
                     prior_pi_additive=c(1,1), ## prior for pi in spike & slab on additive component
                     prior_tau2_additive=c(10,10), ## shape and rate for inverse gamma on spline penalty terms
                     prior_pi_nonadditive=c(1,1), ## prior for pi in spike & slab on non-additive component
                     prior_rho_nonadditive=c(5,5), ## shape and rate for gamma
                     prior_tau2_nonadditive=c(20,20),## this was for the inverse gamma that matches moments with bkmr, which useswith mean 10 sd 10. #c(0.001,0.001), ## shape and rate for inverse gamma on nonadditive GP component
                     prior_sigma2=c(1,1),## shape and rate for inverse gamma on \sigma^2
                     ## MH tuning
                     stepsize_theta_kappa=100, # larger means smaller steps 
                     stepsize_tau2_nonadditive=1,##jumpsize/sd for gamma proposal on tau2_nonadditive
                     stepsize_rho_nonadditive=0.2,##jumpsize/sd for gamma proposal on rho
                     oversample=TRUE,
                     ## approximating large inverse via Sherman-Morrison-Woodbury
                     invmethod=c("GPP"), ## "exact" for full matrix inversion, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
                     rank=200, ## rank of approximation 
                     kernelfun="gaussian", ## choice of kernel function
                     approxProj = FALSE, 
                     draw_h=FALSE,
                     centering=TRUE, ## should covariates be centered to have mean=0
                     scaling=TRUE, 
                     hierarchical=TRUE)
  
  cbmim_pred <- pred_surface(cbmim_mod,Xnew=exposure_list,includeInt=TRUE)
  cbmim_pred_summary <- summarize_pred(cbmim_pred,assoc=FALSE,centered=FALSE)
  cbmim_pred_PSR <- summarize_pred_PSR(cbmim_pred,nchains=cbmim_mod$nchains)
  
  cbmim_chain <- list(cbmim_mod,cbmim_pred,cbmim_pred_summary,cbmim_pred_PSR)
  save(cbmim_chain,file=paste0("cbmim_bmi_base_post_supplementary.Rdata"))
  
  ## labels
  index_names <- names(exposure_list)


  ## indexwise curves
  cbmim_pred_ind_raw <- try(pred_surface_indexwise(cbmim_mod,includeInt=FALSE)) ## catching errors
  cbmim_pred_ind <- summarize_pred(cbmim_pred_ind_raw,assoc=F)
  cbmim_assoc_ind <- summarize_pred(cbmim_pred_ind_raw,assoc=T)


  ## indexwise PIPs (main effects+interactions)
  cbmim_PIPs <- get_PIPs(cbmim_mod)
  
  # ## save all results
  cbmim_chain <- list(mod=cbmim_mod, 
                      index_names=index_names, 
                      pred=cbmim_pred_summary,
                      pred_ind=cbmim_pred_ind, 
                      PIPs=cbmim_PIPs)

  save(cbmim_chain,file=paste0("cbmim_bmi_base_post_supplementary",".Rdata"))
  
  
  
  
  
  
  
  ##########################
  ## RESULTS AND FIGURES ###
  ##########################
  
  library(bsmim2)
  library(tidyverse)
  library(reshape)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  
  # load results from full model
  load("cbmim_bmi_base_post_supplementary.Rdata")
  ls()
  
  ##
  index_names <- cbmim_chain$index_names
  index_names[13:21] <- c("As","Cd","Co","Cs","Cu","Hg","Mn","Mo","Pb")
  cbmim_PIPs <- cbmim_chain$PIPs
  cbmim_pred_ind <- cbmim_chain$pred_ind
  cbmim_pred_inter <- cbmim_chain$pred_inter
  
  #------------------------------
  # Clean names
  #------------------------------
  index_names <- str_replace_all(index_names,"Pregnancy","Preg")
  index_names <- str_replace_all(index_names,"_Postnatal","")
  
  index_names <- str_replace_all(index_names,"Organophosphatepesticides","Organophosphates")
  index_names <- str_replace_all(index_names,"Polybrominateddiphenylethers\\(PBDE\\)","PBDE")
  index_names <- str_replace_all(index_names,"Per-andpolyfluoroalkylsubstances\\(PFAS\\)","PFAS")
  index_names <- str_replace_all(index_names,"Socialandeconomiccapital","SECapital")
  
  #------------------------------
  # Plot PIPs
  #------------------------------
  
  table_PIPs <- data.frame(Index=str_sub(index_names, end=-6),
                           Lm=cbmim_chain$mod$Lm,
                           PIP=round(100*cbmim_PIPs$mainPIPs)) 
  row.names(table_PIPs) <- NULL
  print(xtable(table_PIPs, include.rownames=FALSE))
  
  
  ## data for scatterplots
  dfPIPs <- data.frame(Index=1:length(index_names),
                       PIP=cbmim_PIPs$mainPIPs,
                       interPIP=cbmim_PIPs$interPIPs)
  
  ## Main effects
  plot_mainPIPs <- ggplot(dfPIPs,aes(x=Index,y=PIP))+
    geom_point()+
    theme_classic() +
    scale_x_discrete(labels=index_names)+
    ylim(c(0,1))+
    guides(x =  guide_axis(angle = 45)) +
    ggtitle("Main Effects")
  plot_mainPIPs
  
  
  ## Non-additive PIPs
  plot_interPIPs <- ggplot(dfPIPs,aes(x=Index,y=interPIP))+
    geom_point()+
    theme_classic() +
    scale_x_discrete(labels=index_names)+
    ylim(c(0,1))+
    labs(x="Index", y="PIP")+
    guides(x =  guide_axis(angle = 45)) +
    ggtitle("Non-Additive Effects")
  plot_interPIPs
  
  
  ### HEATPLOT OF PIPs
  heat <- cbmim_PIPs$twowayPIPs
  diag(heat) <- cbmim_PIPs$mainPIPs
  colnames(heat) <- rownames(heat) <- index_names
  roundPIPs <- function(val){
    lab <- rep("",length(val)) 
    lab[!is.na(val)&val>0.5] <- paste0(round(100*val[!is.na(val)&val>0.5]))
    return(lab)
  }
  heatplot <- ggplot(data = melt(heat), aes(x=X1, y=X2,
                                            fill=value)) +
    geom_tile()+
    geom_text(aes(X1, X2, label = roundPIPs(value)),
              color = "white", size = 4)+
    # theme_void() +
    labs(x="", y="")+
    guides(x =  guide_axis(angle = 45)) +
    guides(fill=guide_legend(title="PIP"))+
    ggtitle("")
  heatplot
  
  
  #------------------------------
  # Plot index-wise curves
  #------------------------------
  
  ## indexwise (smoothed)
  for (jj in 1:length(cbmim_pred_ind)){
    df <- data.frame(
      gridx=seq(-2,2,length=21),
      fitted=predict(loess(cbmim_pred_ind[[jj]]$mean~seq(-2,2,length=21))),
      uci=predict(loess(cbmim_pred_ind[[jj]]$uci~seq(-2,2,length=21))),
      lci=predict(loess(cbmim_pred_ind[[jj]]$lci~seq(-2,2,length=21)))
    )
    
    pp <- ggplot(df,aes(x=gridx,y=fitted))+
      geom_line(linetype=3)+ ## 3 is dotted
      geom_ribbon(aes_string(ymin="lci",ymax="uci"),alpha=0.2)+
      ylim(min(df$lci)-0.05,max(df$uci)+0.05)+
      ylab(bquote(f[.(jj)](E[.(jj)])))+
      xlab(bquote(E[.(jj)]))+
      # ggtitle(str_sub(index_names, end=-6)[jj])+
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    ggsave(filename=paste0("index",jj,"_supplementary.pdf"),plot=pp,width=4,height=4)
  }
  
  
  ## Plot index weight posteriors
  ## weights for index 6
  id6 <- (cumsum(cbmim_chain$mod$Lm)[6-1]+1):cumsum(cbmim_chain$mod$Lm)[6]
  nms6 <- c("DDE","DDT","HCB","PCB118","PCB138","PCB153","PCB170","PCB180")
  ylims <- c(-0.8,1)
  boxplot(cbmim_chain$mod$theta[,id6],ylim=ylims,ylab="Index Weight",main="Organochlorines",names=nms6,las=2);abline(h=0,lty=2)
  
  library(xtable)
  xtable(cbind(nms6,t(apply(cbmim_chain$mod$theta[,id6],2,function(x)round(quantile(x,c(0.025,0.5,0.975)),3)))))
  
  
  