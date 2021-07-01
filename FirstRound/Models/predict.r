rm(list=ls())
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=predict_range)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, layer_df$Y, sep="_")

print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  #target_folders<-c(sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp),
  #                 sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folders<-c(sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/predict", target_folder)
    if (dir.exists(target)){
      next()
    }
    dir.create(target, showWarnings = F, recursive = T)
    
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    print(paste(group, i, nrow(df_list), item$sp, target))
    
    j=1
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      
      env_item<-env_layers[[layer_item$LABEL]]
      
      env_item$in_out<-(between(env_item$PR, model$range_PR_min, model$range_PR_max)&
                          (between(env_item$TEMP, model$range_TEMP_min, model$range_TEMP_max)))
      
      
      #env_item<-env_item%>%dplyr::filter(in_out_1&in_out_2)
      env_item<-env_item%>%dplyr::filter(in_out)
      if (F){
        plot(env_layers[[layer_item$LABEL]]$TEMP, env_layers[[layer_item$LABEL]]$PR, pch=".", col="grey")
        colors<-c("blue", "red")
        plot(env_item$TEMP, env_item$PR, col=colors[as.numeric(env_item$in_out)+1], pch=".")
      }
      saveRDS(env_item, sprintf("%s/%s.rda", target, layer_item$LABEL))
    }
  }
  
}
