library(raster)
#library(rgdal)
#library(rgeos)
library(MASS)
library(dplyr)
#library(cluster)
library(data.table)
library(batchtools)
#library(ggplot2)
rm(list=ls())
#in_Ellipsoid <- stats::qchisq(0.95, 2)
#setDTthreads(1)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")
start_env_layers<-readRDS("../../Objects_Full_species/stacked_layers_1850_2020_df.rda")
start_env_layers$diff<-start_env_layers$TEMP_MAX-start_env_layers$TEMP_MIN
if (F){
  
  ggplot(start_env_layers[year==1850])+geom_tile(aes(x=x, y=y, fill=TEMP_MIN))
}
df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
i=2

start_env_layers<-start_env_layers[year>=1970]
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
#group<-"Birds"
#item<-data.table(sp="Cinnyris_jugularis", area=1)

for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  
  target_folder<-sprintf("../../Objects_Full_species/Niche_Models_1850_1970/%s/%s", group, item$sp)
  
  if (dir.exists(target_folder)){
    next()
  }
  print(paste(i, nrow(df_list)))
  
  dir.create(target_folder, showWarnings = F, recursive = T)
  occ<-readRDS(sprintf("../../Objects_Full_species/IUCN_Distribution/%s/%s.rda", group, item$sp))
  v<-raster::extract(mask, occ[, c("x", "y")])
  v<-start_env_layers %>%dplyr::filter(mask_index %in% v)
  j=1
  
  all_v<-v%>%dplyr::group_by(year, x, y)%>%dplyr::summarise(PR=mean(PR),
                                                            TEMP_MAX=mean(TEMP_MAX),
                                                            TEMP_MIN=mean(TEMP_MIN))
  
  if (F){
    ggplot()+geom_tile(data=mask_p, aes(x=x, y=y))+
      geom_tile(data=all_v, aes(x=x, y=y, fill="red"))
  }
  #fit <- cov.rob(all_v[, c("PR", "TEMP_MAX", "TEMP_MIN")], 
  #               quantile.used=NDquntil(nrow(all_v), 0.95),  method = "mve")
  
  mean_TEMP_MAX<-mean(all_v$TEMP_MAX, na.rm=T)
  mean_TEMP_MIN<-mean(all_v$TEMP_MIN, na.rm=T)
  sd_TEMP_MAX<-sd(all_v$TEMP_MAX, na.rm=T)
  sd_TEMP_MIN<-sd(all_v$TEMP_MIN, na.rm=T)
  quantile_TEMP_MAX<-quantile(all_v$TEMP_MAX, c(0.25, 0.75), na.rm=T)
  IQR_TEMP_MAX<-quantile_TEMP_MAX[2]-quantile_TEMP_MAX[1]
  quantile_TEMP_MIN<-quantile(all_v$TEMP_MIN, c(0.25, 0.75), na.rm=T)
  IQR_TEMP_MIN<-quantile_TEMP_MIN[2]-quantile_TEMP_MIN[1]
  range_TEMP_sd_min<-mean_TEMP_MIN-3*sd_TEMP_MIN
  range_TEMP_sd_max<-mean_TEMP_MAX+3*sd_TEMP_MAX
  range_TEMP_IQR_min<-mean_TEMP_MIN-1.5*IQR_TEMP_MIN
  range_TEMP_IQR_max<-mean_TEMP_MAX+1.5*IQR_TEMP_MAX
  
  
  mean_PR<-mean(all_v$PR, na.rm=T)
  sd_PR<-sd(all_v$PR, na.rm=T)
  quantile_PR<-quantile(all_v$PR, c(0.25, 0.75), na.rm=T)
  IQR_PR<-quantile_PR[2]-quantile_PR[1]
  range_PR_sd_min<-mean_PR-3*sd_PR
  range_PR_sd_max<-mean_PR+3*sd_PR
  range_PR_IQR_min<-mean_PR-1.5*IQR_PR
  range_PR_IQR_max<-mean_PR+1.5*IQR_PR
  
  if (F){
    ggplot(all_v)+geom_point(aes(x=PR, y=TEMP_MAX))+
      geom_point(aes(x=PR, y=TEMP_MIN))+
      geom_hline(yintercept=range_TEMP_sd_min)+
      geom_hline(yintercept=range_TEMP_sd_max)+
      geom_hline(yintercept=range_TEMP_IQR_min, color="red")+
      geom_hline(yintercept=range_TEMP_IQR_max, color="red")+
      geom_vline(xintercept=range_PR_sd_min)+
      geom_vline(xintercept=range_PR_sd_max)+
      geom_vline(xintercept=range_PR_IQR_min, color="red")+
      geom_vline(xintercept=range_PR_IQR_max, color="red")
    
  }
  
  min_x<-min(all_v$x, na.rm=T)
  max_x<-max(all_v$x, na.rm=T)
  mean_x<-mean(all_v$x, na.rm=T)
  
  min_y<-min(all_v$y, na.rm=T)
  max_y<-max(all_v$y, na.rm=T)
  max_abs_y<-max(abs(all_v$y), na.rm=T)
  mean_y<-mean(all_v$y, na.rm=T)
  
  
  N_CELL<-nrow(all_v%>%dplyr::filter(year==1970))
  fit<-data.frame(
    mean_TEMP_MAX=mean_TEMP_MAX,
    mean_TEMP_MIN=mean_TEMP_MIN,
    sd_TEMP_MAX=sd_TEMP_MAX,
    sd_TEMP_MIN=sd_TEMP_MIN,
    quantile_TEMP_MAX_low=quantile_TEMP_MAX[1],
    quantile_TEMP_MAX_high=quantile_TEMP_MAX[2],
    IQR_TEMP_MAX=IQR_TEMP_MAX,
    quantile_TEMP_MIN_low=quantile_TEMP_MIN[1],
    quantile_TEMP_MIN_high=quantile_TEMP_MIN[2],
    IQR_TEMP_MIN=IQR_TEMP_MIN,
    range_TEMP_sd_min=range_TEMP_sd_min,
    range_TEMP_sd_max=range_TEMP_sd_max,
    range_TEMP_IQR_min=range_TEMP_IQR_min,
    range_TEMP_IQR_max=range_TEMP_IQR_max,
    mean_PR=mean_PR,
    sd_PR=sd_PR,
    quantile_PR_low=quantile_PR[1],
    quantile_PR_high=quantile_PR[2],
    IQR_PR=IQR_PR,
    range_PR_sd_min=range_PR_sd_min,
    range_PR_sd_max=range_PR_sd_max,
    range_PR_IQR_min=range_PR_IQR_min,
    range_PR_IQR_max=range_PR_IQR_max,
    min_x=min_x,
    max_x=max_x,
    mean_x=mean_x,
    min_y=min_y,
    max_y=max_y,
    mean_y=mean_y,
    max_abs_y=max_abs_y,
    N_CELL=N_CELL
  )
  saveRDS(fit, sprintf("%s/fit.rda", target_folder))
  
  all_v$in_out<-0
  all_v<-data.table(all_v)
  all_v[((PR %between% c(fit$range_PR_sd_min, fit$range_PR_sd_max))&
           (TEMP_MAX %between% c(fit$range_TEMP_sd_min, fit$range_TEMP_sd_max))&
           (TEMP_MIN %between% c(fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)))]$in_out<-1
  
  saveRDS(all_v, sprintf("%s/occ_with_env.rda", target_folder))
}

