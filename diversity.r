library(raster)
#library(rgdal)
#library(rgeos)
#library(MASS)
library(data.table)
library(dplyr)
library(ggplot2)
rm(list=ls())
source("addEllipse.R")
source("genCircle.R")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Mammals"
}


GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")

predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
j=1
dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, -1), N=c(rep(1,5), c(2:5), 2, 1))
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

if (F){
  mask<-readRDS("../../Objects/stacked_layers_2015_2100.rda")
  xy<-unique(mask[[1]][, c("x", "y")])
  xy$label<-paste(xy$x, xy$y, sep=",")
  saveRDS(xy, "../../Objects/mask.rda")
}
mask<-readRDS("../../Objects/mask.rda")

for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:nrow(dispersals))){
    layer$M<-dispersals[k, "M"]
    layer$N<-dispersals[k, "N"]
    target_folder<-sprintf("../../Objects/Diversity/%s_%d_%d", layer$LABEL, layer$M, layer$N)
    if (dir.exists(target_folder)){
      next()
    }
    dir.create(target_folder, showWarnings = F)
    
    i=247
    diversity_df<-list()
    for (i in c(1:nrow(df_list))){
      item<-df_list[i,]
      item$sp<-gsub(" ", "_", item$sp)
      
      if (item$area<=0){
        next()
      }
      #print(paste(Sys.time(), 1))
      enm_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/Mammals/%s/dispersal", item$sp)
      env_item_all<-readRDS(sprintf("%s/%s_%s_%d_%d.rda", enm_folder, layer$GCM, layer$SSP, layer$M, layer$N))
      env_item_all<-env_item_all[, c("x", "y", "YEAR")]
      #print(paste(Sys.time(), 2))
      if (is.null(env_item_all)){
        next()
      }
      colnames(env_item_all)[3]<-"year"
      print(paste(layer$LABEL, layer$M, layer$N, i, nrow(df_list)))
      splited<-split(data.table(env_item_all), by="year")
      
      for (YYYY in c(2014:2100)){
        #print(paste(Sys.time(), 3, YYYY))
        diversity<-diversity_df[[as.character(YYYY)]]
        #print(paste(Sys.time(), 4, YYYY))
        if (YYYY==2014){
          source_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp)
          start_dis<-readRDS(sprintf("%s/occ_with_env.rda", source_folder))
          start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(X, Y)
          
          colnames(start_dis)<-c("x", "y")
          start_dis$year<-2014
          start_dis$sp<-item$sp
          if (is.null(diversity)){
            diversity<-list(start_dis)
          }else{
            diversity[[length(diversity)+1]]<-start_dis
          }
        }else{
          
          if (layer$M<0){
            enm_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/Mammals/%s/predict", item$sp)
            env_item<-readRDS(sprintf("%s/%s_%s_%d.rda", enm_folder, layer$GCM, layer$SSP, 2100))
            env_item<-env_item[, c("x", "y")]
            env_item$year<-2100
            
            if (F){
              plot(env_item$x, env_item$y)
              points(start_dis$x, start_dis$y, col="red")
            }
            
          }else{
            #print(paste(Sys.time(), 5, YYYY))
            env_item<-splited[[as.character(YYYY)]]
            #print(paste(Sys.time(), 6, YYYY))
            if (is.null(env_item)){
              next()
            }
            
          }
          if (nrow(env_item)==0){
            next()
          }
          #print(paste(Sys.time(), 7, YYYY))
          env_item$sp<-item$sp
          if (is.null(diversity)){
            diversity<-list(env_item)
          }else{
            diversity[[length(diversity)+1]]<-env_item
          }
          #print(paste(Sys.time(), 8, YYYY))
        }
        #print(paste(Sys.time(), 9))
        diversity_df[[as.character(YYYY)]]<-diversity
        #print(paste(Sys.time(), 10))
      }
    }
    saveRDS(diversity_df, sprintf("%s/diversity_df.rda", target_folder))
  }
  
 
  
}

if (F){
  tt<-rbindlist(diversity_df[[69]])
}