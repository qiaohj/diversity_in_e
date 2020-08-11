library(raster)
library(ggplot2)
library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

if (F){
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  VARs<-c("pr", "tasmax", "tasmin")
  Ys<-c(1850:2100)
  
  i=10
  COMs<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=Ys, stringsAsFactors = F)
  template<-"../../Raster/ENV/Monthly/%s/%s/%s/%d/%s_%d.tif"
  template_year<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s.tif"
  
  COMs<-COMs[sample(nrow(COMs)),]
  env_var_year<-NULL
  for (i in c(1:nrow(COMs))){
    print(paste(i, nrow(COMs)))
    com<-COMs[i,]
    label<-""
    if (com$VAR=="pr"){
      label<-"sum"
    }
    if (com$VAR=="tasmax"){
      label<-"max"
    }
    if (com$VAR=="tasmin"){
      label<-"min"
    }
    r<-raster(sprintf(template_year, com$GCM, com$SSP, com$VAR, com$Y, label))
    mean_v<-mean(values(r))
    com$V<-mean_v
    if (is.null(env_var_year)){
      env_var_year<-com
    }else{
      env_var_year<-bind_rows(env_var_year, com)
    }
  }
  saveRDS(env_var_year, "../../Objects/mean_env_year.rda")
}

mean_env_year<-readRDS("../../Objects/mean_env_year.rda")
mean_env_year$label<-paste(mean_env_year$GCM, mean_env_year$SSP)
p<-ggplot(mean_env_year, aes(x=Y, y=V, linetype=factor(GCM), color=factor(SSP)))+
  geom_line()+
  theme_bw()+
  geom_vline(xintercept = 2014)+
  facet_wrap(~VAR, nrow=3, scale="free")
ggsave(p, file="../../Figures/Env/GCM_Curves.png", width=10, height=10)




