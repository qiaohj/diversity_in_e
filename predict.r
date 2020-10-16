rm(list=ls())
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
source("addEllipse.R")
source("genCircle.R")
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)

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
var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=predict_range)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, layer_df$Y, sep="_")
if (F){
  mask<-data.frame(rasterToPoints(raster("../../Raster/mask_index.tif")))
  env_layers<-list()
  for (i in c(1:nrow(layer_df))){
    print(paste("Init layer list:", i, nrow(layer_df)))
    item<-layer_df[i,]
    pr<-sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum")
    tasmax<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max")
    tasmin<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min")
    layers<-stack(c(pr, tasmax, tasmin))
    names(layers)<-c("PR", "TEMP_MAX", "TEMP_MIN")
    layers_table<-mask
    layers_table[,c("PR", "TEMP_MAX", "TEMP_MIN")]<-raster::extract(layers, mask[, c("x", "y")])
    env_layers[[item$LABEL]]<-layers_table
  }
  saveRDS(env_layers, "../../Objects/stacked_layers_2015_2100.rda")
}
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2015_2100.rda")

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
  target_folders<-c(sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/predict", target_folder)
    if (dir.exists(target)){
      next()
    }
    dir.create(target, showWarnings = F)
    
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    print(paste(group, i, nrow(df_list), item$sp, target))
    
    j=1
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      
      env_item<-env_layers[[layer_item$LABEL]]
      #env_item$dist1 <- stats::mahalanobis(env_item[, c("PR", "TEMP_MIN")], center = model$center, 
      #                                     cov = model$cov)
      env_item$dist2 <- stats::mahalanobis(env_item[, c("PR", "TEMP_MAX")], center = model$center, 
                                           cov = model$cov)
      
      #env_item$in_out_1<-F
      #env_item[which(env_item$dist1<in_Ellipsoid), "in_out_1"]<-T
      env_item$in_out_2<-F
      env_item[which(env_item$dist2<in_Ellipsoid), "in_out_2"]<-T
      
      if (F){
        df_ori<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
        colors<-c("red", "blue")
        plot(env_item$PR, env_item$TEMP_MAX, col=colors[env_item$in_out_2+1], 
             xlim=range(env_item$PR), ylim=range(c(env_item$TEMP_MAX)), pch=".")
             #xlim=range(env_item$PR), ylim=range(c(env_item$TEMP_MAX, env_item$TEMP_MIN)), pch=".")
        #points(env_item$PR, env_item$TEMP_MIN, col=colors[env_item$in_out_1+1], pch=".")
        
        addEllipse(model$center, model$cov, col="black", p.interval=0.95)
        points(df_ori$PR, df_ori$TEMP, col="purple")
        env_item_present<-env_item%>%dplyr::filter(in_out_1&in_out_2)
        plot(raster("../../Raster/mask.tif"), col="grey")
        points(env_item_present$x, env_item_present$y, pch=".", col="black")
        points(df_ori$X, df_ori$Y, col="purple")
      }
      #env_item<-env_item%>%dplyr::filter(in_out_1&in_out_2)
      env_item<-env_item%>%dplyr::filter(in_out_2)
      saveRDS(env_item, sprintf("%s/%s.rda", target, layer_item$LABEL))
    }
  }
  
}