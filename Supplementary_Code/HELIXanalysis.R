## set up
options(echo=TRUE)
options(stringsAsFactors = FALSE)
exp_set="postnatal" 


library(tidyverse)
devtools::install_github("glenmcgee/bsmim2")
library(bsmim2)
library(MASS)
library(parallel)
devtools::install_github("glenmcgee/CKMR")
library(CKMR) ## install_github("glenmcgee/CKMR")
load("exposome.RData")
load("covariates.rda") ## from Make_Covariates.R
load("index_exposure_list.rda") ## from Make_Exposures.R

################
## exposures ###
################
exposure_list <- index_exposure_list[grepl("Post",names(index_exposure_list))]
exposure_list <- exposure_list[-c(2,14)]
exposure_list <- lapply(exposure_list,scale)

################
## covariates ##
################
z <- cbind(z_base_preg,z_base_post)
# remove covariates that are already adjusted for in the outcome or related to the outcome
z <- z[,-which(colnames(z)%in%c("hs_child_age_None","e3_sex_Nonemale","hs_c_weight_None","hs_c_height_None"))]


###############
##  Fit CMIM ##
###############
nit <- 20000
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
                     prior_tau2_additive=c(10,10),## shape and rate for inverse gamma on spline penalty terms
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
  save(cbmim_chain,
       file=paste0("CBMIM/Output/cbmim_bmi_base_post","_newlongchain.Rdata"))
  
  ## labels
  index_names <- names(exposure_list)


  ## indexwise curves
  cbmim_pred_ind_raw <- try(pred_surface_indexwise(cbmim_mod,includeInt=FALSE)) ## catching errors
  cbmim_pred_ind <- summarize_pred(cbmim_pred_ind_raw,assoc=F)

  ## indexwise PIPs (main effects+interactions)
  cbmim_PIPs <- get_PIPs(cbmim_mod)
  
  # ## save all results
  cbmim_chain <- list(mod=cbmim_mod, 
                      index_names=index_names, 
                      pred=cbmim_pred_summary,
                      pred_ind=cbmim_pred_ind, 
                      PIPs=cbmim_PIPs)

  save(cbmim_chain,file=paste0("cbmim_bmi_base_post.Rdata"))
  
  
  
  ########################
  ## Standard BMIM Fit  ##
  ######################## 
  nit <- 20000
  thin = 5
  burnpct <- 0.5
  set.seed(1234)
  
  bmim_mod <- bsmim2(y=scale(phenotype$hs_zbmi_who), 
                     x=exposure_list, 
                     z=z, 
                     niter=nit,nburn=burnpct*nit,nthin=thin,
                     centering=T,scaling=T,
                     prior_sigma=c(0.001,0.001),
                     prior_lambda_shaperate=c(1,0.1),
                     gaussian=TRUE,
                     spike_slab=TRUE,
                     gauss_prior=TRUE,
                     prior_theta_slab_sd=0.25,
                     stepsize_theta=0.1,
                     draw_h=FALSE,
                     num_theta_steps = 10) 
 
  index_names <- colnames(index_exposure_list)
  
  bmim_pred <- predict_hnew_X2(bmim_mod,newX=exposure_list)#,newY=y_TEST,newZ=dat$covariates_TEST)
  bmim_pred_summary <- bmim_pred$fits[,c("mean","lower","upper")]
  colnames(bmim_pred_summary) <- c("mean","lci","uci")
  
  bmim_pred_ind <- predict_hnew_indexwise2(bmim_mod)
  
  bmim_PIPs <- summarize_thetas(bmim_mod)

  
  pred_twoway <- function(obj,whichindex=c(4,5,7),qtls=c(0.1,0.9),qtl_lims=c(0.01,0.99),pts=20,include_median=TRUE){
    
    ## get predictions at each level (skip 0.5 since it is implicitly computed)
    res <- list()
    for(mm in whichindex){
      res_mm <- list()
      for(qq in 1:length(qtls)){
        res_mm[[qq]] <- predict_hnew_indexwise2(obj,crossM=mm,qtl=qtls[[qq]],qtl_lims=qtl_lims,points=pts)
      }
      names(res_mm) <- qtls  ## label
      res[[mm]] <- res_mm
    }
    
    ## combine predictions into dataframe for plotting
    df_var1 <- df_var2 <- df_grid <- df_quantile <- df_est <- df_est_centered <- c()
    df_lower <- df_lower_centered <- df_upper <- df_upper_centered <- c()
    for(xx in (whichindex)){
      for(yy in (whichindex)){
        if(xx==yy){
          next
        }
        for(qq in 1:length(qtls)){
          df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
          df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
          df_grid <- c(df_grid,res[[yy]][[qq]]$grid[(xx-1)*pts+1:pts])
          df_est <- c(df_est,res[[yy]][[qq]]$mean[(xx-1)*pts+1:pts])
          df_est_centered <- c(df_est_centered,res[[yy]][[qq]]$mean_centered[(xx-1)*pts+1:pts])
          df_lower <- c(df_lower,res[[yy]][[qq]]$lower[(xx-1)*pts+1:pts])
          df_lower_centered <- c(df_lower_centered,res[[yy]][[qq]]$lower_centered[(xx-1)*pts+1:pts])
          df_upper <- c(df_upper,res[[yy]][[qq]]$upper[(xx-1)*pts+1:pts])
          df_upper_centered <- c(df_upper_centered,res[[yy]][[qq]]$upper_centered[(xx-1)*pts+1:pts])
          df_quantile <- c(df_quantile,rep(qtls[qq],pts))
        }
        if(include_median==TRUE){ ## implicitly predicted above
          df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
          df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
          df_grid <- c(df_grid,res[[xx]][[qq]]$grid[(xx-1)*pts+1:pts]) ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
          df_est <- c(df_est,res[[xx]][[qq]]$mean[(xx-1)*pts+1:pts])   
          df_est_centered <- c(df_est_centered,res[[xx]][[qq]]$mean_centered[(xx-1)*pts+1:pts])
          df_lower <- c(df_lower,res[[xx]][[qq]]$lower[(xx-1)*pts+1:pts])   
          df_lower_centered <- c(df_lower_centered,res[[xx]][[qq]]$lower_centered[(xx-1)*pts+1:pts])
          df_upper <- c(df_upper,res[[xx]][[qq]]$upper[(xx-1)*pts+1:pts])   
          df_upper_centered <- c(df_upper_centered,res[[xx]][[qq]]$upper_centered[(xx-1)*pts+1:pts])
          df_quantile <- c(df_quantile,rep(0.5,pts))
        }
        
      }
    }
    
    pred_df <- data.frame(var1=df_var1,var2=df_var2,grid=df_grid,quantile=df_quantile,est=df_est,est_centered=df_est_centered,
                          lower=df_lower,lower_centered=df_lower_centered,upper=df_upper,upper_centered=df_upper_centered)
    
    return(pred_df)
  }
  
  bmim_pred_inter <- pred_twoway(bmim_mod)
  
  bmim_chain <- list(mod=bmim_mod, 
                      index_names=index_names, 
                      pred=bmim_pred_summary,
                      pred_ind=bmim_pred_ind, 
                      PIPs=bmim_PIPs,
                      pred_inter=bmim_pred_inter)

  save(bmim_chain,
       file=paste0("bmim_bmi_base_select.Rdata"))
  
  
  
  ## PLOT TWO WAY INTERACTIONS
  range <- c(min(bmim_pred_inter$est_centered)-0.1,max(bmim_pred_inter$est_centered)) ## ylims
  red <- "#F21A00"
  ylw <- "#EBCC2A"
  blu <- "#3B9AB2"
  ryb <- c(ylw,red,blu)
  
  df_4 <- data.frame(est=bmim_pred_ind$mean_centered[bmim_pred_ind$m==4],
            lower=bmim_pred_ind$lower_centered[bmim_pred_ind$m==4],
            upper=bmim_pred_ind$upper_centered[bmim_pred_ind$m==4],
            grid=bmim_pred_ind$grid[bmim_pred_ind$m==4])
  pp_index_4 <-  ggplot(df_4, aes(grid)) +
    geom_smooth(aes(y=est),col="black",size=0.4,linetype="dotted") +
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+   
    labs(x="",y="Metals",title="Metals")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  df_5 <- data.frame(est=bmim_pred_ind$mean_centered[bmim_pred_ind$m==5],
                     lower=bmim_pred_ind$lower_centered[bmim_pred_ind$m==5],
                     upper=bmim_pred_ind$upper_centered[bmim_pred_ind$m==5],
                     grid=bmim_pred_ind$grid[bmim_pred_ind$m==5])
  pp_index_5 <-  ggplot(df_5, aes(grid)) +
    geom_smooth(aes(y=est),col="black",size=0.4,linetype="dotted") +
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+   
    labs(x="",y="")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  df_7 <- data.frame(est=bmim_pred_ind$mean_centered[bmim_pred_ind$m==7],
                     lower=bmim_pred_ind$lower_centered[bmim_pred_ind$m==7],
                     upper=bmim_pred_ind$upper_centered[bmim_pred_ind$m==7],
                     grid=bmim_pred_ind$grid[bmim_pred_ind$m==7])
  pp_index_7 <-  ggplot(df_7, aes(grid)) +
    geom_smooth(aes(y=est),col="black",size=0.4,linetype="dotted") +
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2)+   
    labs(x="",y="")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  
  pp_inter_45 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 4" & bmim_pred_inter$var2=="Index 5",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+
    labs(x="",y="",title="Meteorological")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
 
  pp_inter_47 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 4" & bmim_pred_inter$var2=="Index 7",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+
    labs(x="",y="",title="Organochlorines")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  pp_inter_57 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 5" & bmim_pred_inter$var2=="Index 7",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+    
    labs(x="",y="",col="Quantile")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  
  pp_inter_54 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 5" & bmim_pred_inter$var2=="Index 4",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+  
    labs(x="",y="Meteorological",col="Quantile")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  pp_inter_74 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 7" & bmim_pred_inter$var2=="Index 4",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+   
    labs(x="",y="Organochlorines",col="Quantile")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  pp_inter_75 <-  ggplot(bmim_pred_inter[bmim_pred_inter$var1=="Index 7" & bmim_pred_inter$var2=="Index 5",], aes(grid)) +
    geom_smooth(aes(y=est_centered, col = as.factor(quantile)), stat = "identity",fill="white",size=0.5) +
    scale_color_manual(values=ryb,guide="none")+
    labs(x="",y="",col="Quantile")+
    ylim(range)+theme_light()+
    theme(strip.text = element_text(colour = 'black'),plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12,hjust = 0.5),
          axis.title.y = element_text(size=12))
  
  
  (pp_index_4+pp_inter_45+pp_inter_47)/
  (pp_inter_54+pp_index_5+pp_inter_57)/
  (pp_inter_74+pp_inter_75+pp_index_7)
  
  
  
  
  
  
  
  
  ##########################
  ## RESULTS AND FIGURES ###
  ##########################
  
  library(bsmim2)
  library(tidyverse)
  library(reshape)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  
  load("cbmim_bmi_base_post.Rdata")
  ls()
  
  ##
  index_names <- cbmim_chain$index_names
  cbmim_PIPs <- cbmim_chain$PIPs
  cbmim_pred_ind <- cbmim_chain$pred_ind
  
  #------------------------------
  # Clean names
  #------------------------------
  index_names <- str_replace_all(index_names,"Pregnancy","Preg")
  index_names <- str_replace_all(index_names,"Postnatal","Post")
  
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
  
  
  ## prep data for scatterplots
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
    labs(x="", y="")+
    guides(x =  guide_axis(angle = 45)) +
    guides(fill=guide_legend(title="PIP"))+
    ggtitle("")
  heatplot
  
  #------------------------------
  # Plot index-wise curves
  #------------------------------
  
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
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    ggsave(filename=paste0("index",jj,".pdf"),plot=pp,width=4,height=4)
  }
  
  
  
  ## Plot index weight posteriors
  ## weights for index 4 and 7
  id4 <- (cumsum(cbmim_chain$mod$Lm)[4-1]+1):cumsum(cbmim_chain$mod$Lm)[4]
  id7 <- (cumsum(cbmim_chain$mod$Lm)[7-1]+1):cumsum(cbmim_chain$mod$Lm)[7]
  nms4 <- c("As","Cd","Co","Cs","Cu","Hg","Mn","Mo","Pb")
  nms7 <- c("DDE","DDT","HCB","PCB118","PCB138","PCB153","PCB170","PCB180")
  ylims <- c(-0.8,1)
  par(mfrow=c(1,2))
  boxplot(cbmim_chain$mod$theta[,id4],ylim=ylims,ylab="Index Weight",main="Metals",names=nms4,las=2);abline(h=0,lty=2)
  boxplot(cbmim_chain$mod$theta[,id7],ylim=ylims,ylab="Index Weight",main="Organochlorines",names=nms7,las=2);abline(h=0,lty=2)
  
  library(xtable)
  xtable(cbind(nms4,t(apply(cbmim_chain$mod$theta[,id4],2,function(x)round(quantile(x,c(0.025,0.5,0.975)),3)))))
  xtable(cbind(nms7,t(apply(cbmim_chain$mod$theta[,id7],2,function(x)round(quantile(x,c(0.025,0.5,0.975)),3)))))
  
  
  
  #------------------------------
  # CORRPLOT
  #------------------------------
  
  library(reshape2)
  library(ggplot2)
  
  corr_mat <- round(100*(cor(Reduce("cbind",cbmim_chain$mod$x))))
  corr_mat[upper.tri(corr_mat)] <- NA
  melted_corr_mat <- melt(corr_mat) # reshape
  xnames <- c(
    "NO2",
    "PM10",
    "PM2.5",
    "PMAbs",
    
    "PMAbs",
    "Benzene",
    "NO2",
    "PM2.5",
    "BTEX",
    
    "KIDMED",
    "PhysActivity",
    "Sedentary",
    "Sleep",
    
    "Arsenic",
    "Cadmium",
    "Cobalt",
    "Cesium",
    "Copper",
    "Mercury",
    "Manganese",
    "Molybdenum",
    "Lead",
    
    "Humidity",
    "Temperature",
    "UV",
    
    "NDVI(home)",
    "NDVI(school)",
    
    "DDE",
    "DDT",
    "HCB",
    "PCB118",
    "PCB138",
    "PCB153",
    "PCB170",
    "PCB180",
    
    "DEP",
    "DETP",
    "DMP",
    "DMTP",
    
    "PBDE153",
    "PBDE47",
    
    "PFHXS",
    "PFNA",
    "PFOA",
    "PFOS",
    "PFUNDA",
    
    "BPA",
    "BUPA",
    "ETPA",
    "MEPA",
    "OXBE",
    "PRPA",
    "Triclosan",
    
    "MBzP",
    "MECPP",
    "MEHHP",
    "MEHP",
    "MEOHP",
    "MEP",
    "MiBP",
    "MnBP",
    "OHMiNP",
    "OXOMiNP",
    
    "TotalLoad",
    "Density"
  )
  
  indexnames <- c("AirPollution",
                  "Indoorair",
                  "Lifestyle",
                  "Metals",
                  "Meteorological",
                  "NaturalSpaces",
                  "Organochlorines",
                  "Organophosphates",
                  "PBDE",
                  "PFAS",
                  "Phenols",
                  "Phthalates",
                  "Traffic")
  
  ## heatplot
  corrplot <- ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "lightblue",
                         mid="white",
                         high = "red",
                         midpoint=0,
                         na.value = "white",
                         guide = "colorbar",name="")+
    xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text=element_text(size=8))+
    scale_x_discrete(labels=xnames)+
    scale_y_discrete(labels=xnames)
  corrplot
  
  
  ###########################################
  ### BMIM
  
  # load results from full model
  load("bmim_bmi_base_select.Rdata")
  ls()
  
  nms4 <- c("As","Cd","Co","Cs","Cu","Hg","Mn","Mo","Pb")
  nms7 <- c("DDE","DDT","HCB","PCB118","PCB138","PCB153","PCB170","PCB180")
  ylims <- c(-0.8,1)
  par(mfrow=c(1,2))
  boxplot(bmim_mod$theta$Metals_Postnatal,ylim=ylims,ylab="Index Weight",main="Metals",names=nms4,las=2);abline(h=0,lty=2)
  boxplot(bmim_mod$theta$Organochlorines_Postnatal,ylim=ylims,ylab="Index Weight",main="Organochlorines",names=nms7,las=2);abline(h=0,lty=2)
  par(mfrow=c(1,1))
  
  t(round(apply(bmim_mod$theta$Metals_Postnatal,2,function(x) quantile(x,c(0.025,0.5,0.975))),3))
  t(round(apply(bmim_mod$theta$Organochlorines_Postnatal,2,function(x) quantile(x,c(0.025,0.5,0.975))),3))
  
  