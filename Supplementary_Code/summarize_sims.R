
## First combine results:
# source("combineResults.R")
# combine_dat("results_sim_index0_interac0_iter*")
# combine_dat("results_sim_index0_interac1_iter*")
# combine_dat("results_sim_index1_interac0_iter*")
# combine_dat("results_sim_index1_interac1_iter*")

library(tidyverse)
library(ggplot2)
library(reshape)
library(patchwork)
library(xtable)
library(Rfast)

path <- ""

## load results
results_sim_index0_interac0 <- read_csv(paste0(path,"results_sim_index0_interac0.csv"))
results_sim_index0_interac1 <- read_csv(paste0(path,"results_sim_index0_interac1.csv"))
results_sim_index1_interac0 <- read_csv(paste0(path,"results_sim_index1_interac0.csv"))
results_sim_index1_interac1 <- read_csv(paste0(path,"results_sim_index1_interac1.csv"))

results_sim2_index0_interac0 <- read_csv(paste0(path,"results_sim2_index0_interac0.csv"))
results_sim2_index0_interac1 <- read_csv(paste0(path,"results_sim2_index0_interac1.csv"))
results_sim2_index1_interac0 <- read_csv(paste0(path,"results_sim2_index1_interac0.csv"))
results_sim2_index1_interac1 <- read_csv(paste0(path,"results_sim2_index1_interac1.csv"))

results_sim_index0_interac2 <- read_csv(paste0(path,"results_sim_index0_interac2.csv"))
results_sim2_index0_interac2 <- read_csv(paste0(path,"results_sim2_index0_interac2.csv"))

## relabel CKMR and CMIM
results_sim_index0_interac0$method[results_sim_index0_interac0$method=="CBMIM"] <- "CKMR"
results_sim_index0_interac1$method[results_sim_index0_interac1$method=="CBMIM"] <- "CKMR"
results_sim2_index0_interac0$method[results_sim2_index0_interac0$method=="CBMIM"] <- "CKMR"
results_sim2_index0_interac1$method[results_sim2_index0_interac1$method=="CBMIM"] <- "CKMR"

results_sim_index1_interac0$method[results_sim_index1_interac0$method=="CBMIM"] <- "CMIM"
results_sim_index1_interac1$method[results_sim_index1_interac1$method=="CBMIM"] <- "CMIM"
results_sim2_index1_interac0$method[results_sim2_index1_interac0$method=="CBMIM"] <- "CMIM"
results_sim2_index1_interac1$method[results_sim2_index1_interac1$method=="CBMIM"] <- "CMIM"

results_sim_index0_interac2$method[results_sim_index0_interac2$method=="CBMIM"] <- "CKMR"
results_sim2_index0_interac2$method[results_sim2_index0_interac2$method=="CBMIM"] <- "CKMR"


## summaries
summary_index0_interac0 <- results_sim_index0_interac0 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary_index0_interac1 <- results_sim_index0_interac1 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary_index1_interac0 <- results_sim_index1_interac0 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary_index1_interac1 <- results_sim_index1_interac1 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)

summary2_index0_interac0 <- results_sim2_index0_interac0 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary2_index0_interac1 <- results_sim2_index0_interac1 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary2_index1_interac0 <- results_sim2_index1_interac0 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary2_index1_interac1 <- results_sim2_index1_interac1 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)

summary_index0_interac2 <- results_sim_index0_interac2 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)
summary2_index0_interac2 <- results_sim2_index0_interac2 %>% group_by(method) %>% summarise_all(mean,na.rm=TRUE)


### summary tables
tab_index0_interac0 <- summary_index0_interac0 %>% mutate(index="No",interac="No")   %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab_index0_interac1 <- summary_index0_interac1 %>% mutate(index="No",interac="Yes")  %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab_index1_interac0 <- summary_index1_interac0 %>% mutate(index="Yes",interac="No")  %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab_index1_interac1 <- summary_index1_interac1 %>% mutate(index="Yes",interac="Yes") %>% dplyr::select(index,interac,method,mse,bias,width,cvg)

tab2_index0_interac0 <- summary2_index0_interac0 %>% mutate(index="No",interac="No")   %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab2_index0_interac1 <- summary2_index0_interac1 %>% mutate(index="No",interac="Yes")  %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab2_index1_interac0 <- summary2_index1_interac0 %>% mutate(index="Yes",interac="No")  %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab2_index1_interac1 <- summary2_index1_interac1 %>% mutate(index="Yes",interac="Yes") %>% dplyr::select(index,interac,method,mse,bias,width,cvg)

tab_index0_interac2 <- summary_index0_interac2 %>% mutate(index="Yes",interac="Supp") %>% dplyr::select(index,interac,method,mse,bias,width,cvg)
tab2_index0_interac2 <- summary2_index0_interac2 %>% mutate(index="Yes",interac="Supp") %>% dplyr::select(index,interac,method,mse,bias,width,cvg)


print(xtable(rbind(tab_index0_interac0,
                   tab_index0_interac1,
                   tab_index1_interac0,
                   tab_index1_interac1),dig=2),
      include.rownames=FALSE)


print(xtable(rbind(tab2_index0_interac0,
                   tab2_index0_interac1,
                   tab2_index1_interac0,
                   tab2_index1_interac1),dig=2),
      include.rownames=FALSE)

print(xtable(rbind(tab_index0_interac2,
                   tab2_index0_interac2),dig=2),
      include.rownames=FALSE)

### Plotting functions

## boxplots
boxplot_mse <- function(res,title=""){
  box <- ggplot(res, aes(x=method, y=mse,fill=method)) + 
    scale_fill_manual(values=c("#FF0000","#5BBCD6","#00A08A","#F2AD00","#F98400","Gray","White"),
                      guide=F)+
    geom_boxplot(outlier.alpha = 0.2)+
    theme_classic() +
    labs(x="", y="MSE")+
    ggtitle(title)
  return(box)
}


boxplot_mainPIPs <- function(res,title=""){
  
  box <- res %>% pivot_longer(
    cols = starts_with("mainPIPs"),
    names_to = "Component",
    names_prefix = "mainPIPs",
    values_to = "PIP"
  ) %>% ggplot(aes(x=factor(Component,levels=as.numeric(unique(Component))), y=PIP,fill=method))+
    geom_boxplot()+
    theme_classic() +
    labs(x="Component", y="PIP")+
    guides(fill=guide_legend(title="Method"))+
    ggtitle(title)
    
  return(box)
}


boxplot_interPIPs <- function(res,title=""){
  
  res <- res %>% pivot_longer(
    cols = starts_with("interPIPs"),
    names_to = "Component",
    names_prefix = "interPIPs",
    values_to = "PIP") 
  res$Component <- factor(res$Component,
                          levels=as.numeric(unique(res$Component)))

  box <- res %>% ggplot(aes(x=Component, y=PIP,fill=method))+
    geom_boxplot()+
    theme_classic() +
    labs(x="Components", y="PIP")+
    guides(fill=guide_legend(title="Method"))+
    ggtitle(title)
  
  return(box)
}

## heat plots for PIPs
heatplot_PIPs <- function(sumres,title="",axislab="X"){

  mainPIPs <- sumres %>% dplyr::select(starts_with("mainPIPs")) %>% as.numeric()
  interPIPs <- sumres %>% dplyr::select(starts_with("interPIPs")) %>% as.numeric()
  
  heat <- matrix(NA,ncol=length(mainPIPs),nrow=length(mainPIPs))
  heat <- t(lower_tri.assign(heat,interPIPs))
  diag(heat) <- mainPIPs
  melted_corr_mat <- melt(heat)
  box <- ggplot(data = melted_corr_mat, aes(x=factor(Var1), y=factor(Var2),
                                     fill=value)) +
    geom_tile()+
    geom_text(aes(Var1, Var2, label = round(value,2)),
              color = "white", size = 4)+
    theme_minimal() +
    theme(plot.margin=unit(c(0,0,0,0),"pt"))+
    labs(x="", y="")+
    guides(fill=guide_legend(title="PIP"))+
    ggtitle(title)+ guides(fill="none")+
    scale_x_discrete(labels=paste0(axislab,1:10))+
    scale_y_discrete(labels=paste0(axislab,1:10))

  
  return(box)
}


heatplot_PIPs(summary_index0_interac0 %>% dplyr::filter(method=="CKMR"),"CKMR -- Scenario A","X")+
  heatplot_PIPs(summary_index0_interac0 %>% dplyr::filter(method=="NLInter"),"NLInter -- Scenario A","X")

heatplot_PIPs(summary_index0_interac1 %>% dplyr::filter(method=="CKMR"),"CKMR -- Scenario B","X")+
  heatplot_PIPs(summary_index0_interac1 %>% dplyr::filter(method=="NLInter"),"NLInter -- Scenario B","X")

heatplot_PIPs(summary_index1_interac0 %>% dplyr::filter(method=="CMIM"),"CMIM -- Scenario C","E")+
  heatplot_PIPs(summary_index1_interac0 %>% dplyr::filter(method=="BMIM"),"BMIM -- Scenario C","E")

heatplot_PIPs(summary_index1_interac1 %>% dplyr::filter(method=="CMIM"),"CMIM -- Scenario D","E")+
  heatplot_PIPs(summary_index1_interac1 %>% dplyr::filter(method=="BMIM"),"BMIM -- Scenario D","E")


## nonadaptive
heatplot_PIPs(summary_index0_interac0 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptive -- Scenario A","X")+
  heatplot_PIPs(summary_index0_interac1 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptive -- Scenario B","X")

heatplot_PIPs(summary_index1_interac0 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptiveIndex -- Scenario A","E")+
  heatplot_PIPs(summary_index1_interac1 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptiveIndex -- Scenario B","E")



## sigma2=2
heatplot_PIPs(summary2_index0_interac0 %>% dplyr::filter(method=="CKMR"),"CKMR -- Scenario A","X")+
  heatplot_PIPs(summary2_index0_interac0 %>% dplyr::filter(method=="NLInter"),"NLInter -- Scenario A","X")

heatplot_PIPs(summary2_index0_interac1 %>% dplyr::filter(method=="CKMR"),"CKMR -- Scenario B","X")+
  heatplot_PIPs(summary2_index0_interac1 %>% dplyr::filter(method=="NLInter"),"NLInter -- Scenario B","X")

heatplot_PIPs(summary2_index1_interac0 %>% dplyr::filter(method=="CMIM"),"CMIM -- Scenario C","E")+
  heatplot_PIPs(summary2_index1_interac0 %>% dplyr::filter(method=="BMIM"),"BMIM -- Scenario C","E")

heatplot_PIPs(summary2_index1_interac1 %>% dplyr::filter(method=="CMIM"),"CMIM -- Scenario D","E")+
  heatplot_PIPs(summary2_index1_interac1 %>% dplyr::filter(method=="BMIM"),"BMIM -- Scenario D","E")


## nonadaptive
heatplot_PIPs(summary2_index0_interac0 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptive -- Scenario A","X")+
  heatplot_PIPs(summary2_index0_interac1 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptive -- Scenario B","X")

heatplot_PIPs(summary2_index1_interac0 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptiveIndex -- Scenario A","E")+
  heatplot_PIPs(summary2_index1_interac1 %>% dplyr::filter(method=="NonAdaptive"),"NonAdaptiveIndex -- Scenario B","E")


## supplementary
heatplot_PIPs(summary_index0_interac2 %>% dplyr::filter(method=="CKMR"),"CKMR -- Supplementary Sim","X")+
  heatplot_PIPs(summary_index0_interac2 %>% dplyr::filter(method=="NLInter"),"NLInter -- Supplementary Sim","X")

heatplot_PIPs(summary2_index0_interac2 %>% dplyr::filter(method=="CKMR"),"CKMR -- Supplementary Sim","X")+
  heatplot_PIPs(summary2_index0_interac2 %>% dplyr::filter(method=="NLInter"),"NLInter -- Supplementary Sim","X")



## exclude mixselect and nonadaptive for plots
results_sim_index0_interac0 <- results_sim_index0_interac0[!results_sim_index0_interac0$method%in%c("mixselect","NonAdaptive"),]
results_sim_index0_interac1 <- results_sim_index0_interac1[!results_sim_index0_interac1$method%in%c("mixselect","NonAdaptive"),]
results_sim_index1_interac0 <- results_sim_index1_interac0[!results_sim_index1_interac0$method%in%c("NonAdaptive"),]
results_sim_index1_interac1 <- results_sim_index1_interac1[!results_sim_index1_interac1$method%in%c("NonAdaptive"),]

results_sim2_index0_interac0 <- results_sim2_index0_interac0[!results_sim2_index0_interac0$method%in%c("mixselect","NonAdaptive"),]
results_sim2_index0_interac1 <- results_sim2_index0_interac1[!results_sim2_index0_interac1$method%in%c("mixselect","NonAdaptive"),]
results_sim2_index1_interac0 <- results_sim2_index1_interac0[!results_sim2_index1_interac0$method%in%c("NonAdaptive"),]
results_sim2_index1_interac1 <- results_sim2_index1_interac1[!results_sim2_index1_interac1$method%in%c("NonAdaptive"),]

results_sim_index0_interac2 <- results_sim_index0_interac2[!results_sim_index0_interac2$method%in%c("NonAdaptive"),]
results_sim2_index0_interac2 <- results_sim2_index0_interac2[!results_sim2_index0_interac2$method%in%c("NonAdaptive"),]


## make plots
boxplot_mse(results_sim_index0_interac0,"MSE -- Scenario A")+
  boxplot_mse(results_sim_index0_interac1,"MSE -- Scenario B")

boxplot_mse(results_sim_index1_interac0,"MSE -- Scenario C")+
  boxplot_mse(results_sim_index1_interac1,"MSE -- Scenario D")


boxplot_mse(results_sim2_index0_interac0,"MSE -- Scenario A")+
  boxplot_mse(results_sim2_index0_interac1,"MSE -- Scenario B")

boxplot_mse(results_sim2_index1_interac0,"MSE -- Scenario C")+
  boxplot_mse(results_sim2_index1_interac1,"MSE -- Scenario D")


boxplot_mse(results_sim_index0_interac2,"MSE -- Supplementary Sim")

boxplot_mse(results_sim2_index0_interac2,"MSE -- Supplementary Sim")






### Plotting two-way interactions

inter_sim_index0_interac1 <- read_csv(paste0(path,"inter_sim_index0_interac1.csv"))
inter_sim2_index0_interac1 <- read_csv(paste0(path,"inter_sim2_index0_interac1.csv"))
inter_sim_index1_interac1 <- read_csv(paste0(path,"inter_sim_index1_interac1.csv"))
inter_sim2_index1_interac1 <- read_csv(paste0(path,"inter_sim2_index1_interac1.csv"))


spag_plots <- function(inter,ind=FALSE){
  xplots <- yplots <- vector(mode = "list", length = 3)
  
  ## baselines for true ERF
  h2 <- 0
  h3 <- 4*dt(2*0,df=10) 
  h4 <- sin(0)  
  h6 <- -0 
  
  for(jj in 1:3){
    qq <- c(0.1,0.5,0.9)[jj]
    grid <- seq(-1.5,1.5,by=0.01)
    
    dfx <- inter%>%dplyr::filter(var1=="X1",qntl==qq)
    h1 <- 2*cos(2*grid)
    h5 <- median(dfx$X5)^2 
    xtrue <- h1+h2+h3+h4+h5+h6+0.5*h1*h5
    dfxtrue <- data.frame(xtrue,grid)
    xplots[[jj]] <- ggplot2::ggplot(dfx,ggplot2::aes(x=X1,y=pred))+
      ggplot2::geom_line(aes(group=dfx$iter_no),alpha=0.15,col="lightblue")+
      ggplot2::geom_line(data=dfxtrue,ggplot2::aes(x=grid,y=xtrue))+
      ggplot2::scale_x_continuous(expand=c(0,0),limits=c(-1.5,1.5))+
      ggplot2::scale_y_continuous(limits=c(-3,10))+
      ggplot2::ylab("Est")+
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank())
    if(ind==TRUE){
      xplots[[jj]] <- xplots[[jj]]+ggplot2::scale_x_continuous(name="E1",expand=c(0,0),limits=c(-1.5,1.5))
    }
    
    dfy <- inter%>%dplyr::filter(var1=="X5",qntl==qq)
    h1 <- 2*cos(2*median(dfy$X1))
    h5 <- grid^2 
    ytrue <- h1+h2+h3+h4+h5+h6+0.5*h1*h5
    dfytrue <- data.frame(ytrue,grid)
    yplots[[jj]] <- ggplot2::ggplot(dfy,ggplot2::aes(x=X5,y=pred))+
      ggplot2::geom_line(aes(group=dfx$iter_no),alpha=0.15,col="lightblue")+
      ggplot2::geom_line(data=dfytrue,ggplot2::aes(x=grid,y=ytrue))+
      ggplot2::scale_x_continuous(expand=c(0,0),limits=c(-1.5,1.5))+
      ggplot2::scale_y_continuous(limits=c(-3,10))+
      ggplot2::ylab("Est")+
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank())
    if(ind==TRUE){
      yplots[[jj]] <- yplots[[jj]]+ggplot2::scale_x_continuous(name="E5",expand=c(0,0),limits=c(-1.5,1.5))
    }
  }
  
  return((xplots[[1]]+xplots[[2]]+xplots[[3]])/
           (yplots[[1]]+yplots[[2]]+yplots[[3]]))
}

spag_plots(inter_sim_index0_interac1)
spag_plots(inter_sim2_index0_interac1)

spag_plots(inter_sim_index1_interac1,ind=TRUE)
spag_plots(inter_sim2_index1_interac1,ind=TRUE)


