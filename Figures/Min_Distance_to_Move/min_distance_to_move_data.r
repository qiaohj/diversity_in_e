library(raster)
library(dplyr)
library(ggplot2)
g<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/functions.r")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Amphibians"
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
j=1
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")


predict_range<-c(2021:2100)
year<-2021
for (year in predict_range){
  print(paste(group, year))
  target_rda<-sprintf("../../Objects/Min_distance_to_Dispersal/%s/%d.rda", group, year)
  if (file.exists(target_rda)){
    print("skip")
    next()
  }
  saveRDS(NULL, target_rda)
  result<-NULL
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area<=0){
      next()
    }
    target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
    start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
    start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
    colnames(start_dis)<-c("x", "y")
    target<-sprintf("%s/dispersal", target_folder)
    print(paste(group, year, item$sp, i, ":", nrow(df_list)))
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      
      env_item<-env_layers[[paste(layer_item$LABEL, year, sep="_")]]
      
      env_item$in_out<-(between(env_item$PR, model$range_PR_sd_min, model$range_PR_sd_max)&
                          (between(env_item$TEMP, model$range_TEMP_sd_min, model$range_TEMP_sd_max)))
      
      
      #env_item<-env_item%>%dplyr::filter(in_out_1&in_out_2)
      env_item<-env_item%>%dplyr::filter(in_out)
      
      
      end_dis<-env_item
      if (nrow(end_dis)==0){
        layer_item$dist_min<-Inf
      }else{
        end_dis<-end_dis%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, start_dis)/1e+03)
        layer_item$dist_min<-min(end_dis$dist, na.rm = T)
      }
      layer_item$sp<-item$sp
      layer_item$group<-item$group
      layer_item$year<-year
      result<-bind_dplyr(result, layer_item)
    }
    
    
  }
  saveRDS(result, target_rda)
  
}
