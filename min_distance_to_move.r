library(raster)
library(dplyr)
library(ggplot2)
g<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")


args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Amphibians"
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}

layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
source("functions.r")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
j=1
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

predict_range<-c(2015:2100)

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
    target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp)
    start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
    start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
    colnames(start_dis)<-c("x", "y")
    
    target<-sprintf("%s/dispersal", target_folder)
    model<-"Mean"
    print(paste(group, year, item$sp, i, ":", nrow(df_list)))
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      enm_folder<-sprintf("%s/predict", target_folder)
      end_dis<-readRDS(sprintf("%s/%s_%d.rda", enm_folder, layer_item$LABEL, year))
      if (nrow(end_dis)==0){
        layer_item$dist_min<-Inf
      }else{
        end_dis<-end_dis%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, start_dis)/1e+03)
        layer_item$dist_min<-min(end_dis$dist, na.rm = T)
      }
      layer_item$sp<-item$sp
      layer_item$group<-item$group
      layer_item$year<-year
      result<-bind(result, layer_item)
    }
    
    
  }
  saveRDS(result, target_rda)
  
}
