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


head(mean_env_year)
mean_env_year_item<-mean_env_year%>%dplyr::filter(Y>2014)
env_se<-mean_env_year_item%>%dplyr::group_by(GCM, SSP, VAR, label)%>%
  dplyr::summarise(min_V=min(V),
                   max_V=max(V))
env_se$speed<-(env_se$max_V-env_se$min_V)/(2100-2014)
env_se_mean<-env_se%>%dplyr::group_by(SSP, VAR)%>%
  dplyr::summarise(mean_speed=mean(speed),
                   sd_speed=sd(speed))

#type                          direction mean_speed sd_speed
#<chr>                             <dbl>      <dbl>    <dbl>
#1 Maximum Monthly Precipitation        -1    0.00651  0.00422
#2 Maximum Monthly Precipitation         1    0.00561  0.00271
#3 Maximum Monthly Temperature          -1    0.0302   0.0228 
#4 Maximum Monthly Temperature           1    0.0202   0.00956
#5 Minimum Monthly Temperature          -1    0.0269   0.0137 
#6 Minimum Monthly Temperature           1    0.0259   0.0121
env_se_mean$history_speed<-0
env_se_mean[which(env_se_mean$VAR=="tasmax"), "history_speed"]<-0.0202/100
env_se_mean[which(env_se_mean$VAR=="tasmin"), "history_speed"]<-0.0259/100
env_se_mean[which(env_se_mean$VAR=="pr"), "history_speed"]<-0.00561/100 * 365
env_se_mean$times<-env_se_mean$mean_speed/env_se_mean$history_speed
env_se_mean$