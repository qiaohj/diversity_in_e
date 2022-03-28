library(raster)
library(dplyr)
library(ggplot2)
g<-"Mammals"
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/functions.r")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Mammals"
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
i=100
j=1
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")


predict_range<-c(2021:2100)
YYYY<-2021
for (YYYY in predict_range){
  print(paste(group, YYYY))
  target_rda<-sprintf("../../Objects/Min_distance_to_Dispersal/%s/%d_10km.rda", group, YYYY)
  if (file.exists(target_rda)){
    print("skip")
    next()
  }
  saveRDS(NULL, target_rda)
  result<-NULL
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$SP<-gsub(" ", "_", item$SP)
    
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, item$SP)
    if (!file.exists(sprintf("%s/initial_disp_exposure_0_dispersal_0.rda", target_folder))){
      next()
    }
    start_dis<-readRDS(sprintf("%s/initial_disp_exposure_0_dispersal_0.rda", target_folder))
    start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
    colnames(start_dis)<-c("x", "y")
    target<-sprintf("%s/dispersal", target_folder)
    print(paste(group, YYYY, item$SP, i, ":", nrow(df_list)))
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      
      env_item<-env_layers[[layer_item$LABEL]]
      
      env_item$in_out<-between(env_item$bio1, model$range_bio1_sd_min, model$range_bio1_sd_max)&
                          between(env_item$bio5, model$range_bio5_sd_min, model$range_bio5_sd_max)&
                          between(env_item$bio6, model$range_bio6_sd_min, model$range_bio6_sd_max)&
                          between(env_item$bio12, model$range_bio12_sd_min, model$range_bio12_sd_max)&
                          between(env_item$bio13, model$range_bio13_sd_min, model$range_bio13_sd_max)&
                          between(env_item$bio14, model$range_bio14_sd_min, model$range_bio14_sd_max)
      
      
      #env_item<-env_item%>%dplyr::filter(in_out_1&in_out_2)
      env_item<-env_item%>%dplyr::filter(in_out)
      
      if (nrow(env_item)>0){
        env_item<-env_item%>%dplyr::filter(year==YYYY)
      }
      end_dis<-env_item
      if (nrow(end_dis)==0){
        layer_item$dist_min<-Inf
      }else{
        end_dis<-end_dis%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, start_dis)/1e+03)
        layer_item$dist_min<-min(end_dis$dist, na.rm = T)
      }
      layer_item$SP<-item$SP
      layer_item$group<-item$group
      layer_item$year<-YYYY
      result<-bind_dplyr(result, layer_item)
    }
    
    
  }
  saveRDS(result, target_rda)
  
}
