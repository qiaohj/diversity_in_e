library(raster)
#library(rgdal)
#library(rgeos)
#library(MASS)
#library(cluster)
library(dplyr)
#library(ggplot2)
source("addEllipse.R")
source("genCircle.R")
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}
bind<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-dplyr::bind_rows(df1, df2)
  }
  return(df1)
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

predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
dispersals<-data.frame(M=-1, N=1)
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  #target_folders<-c(sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp),
  #                  sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folders<-c(sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/dispersal", target_folder)
    model<-"Normal"
    if (grepl("Mean_GCM", target)){
      model<-"Mean"
    }
    
    j=1
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      year<-predict_range[1]
      enm_folder<-sprintf("%s/predict", target_folder)
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
      colnames(start_dis)<-c("x", "y")
      k=1
      
      
      for (k in c(1:nrow(dispersals))){
        dispersal<-dispersals[k,]
        prev_dis<-start_dis
        dispersal_log<-NULL
        target_save<-sprintf("%s/%s_%d_%d.rda", target, layer_item$LABEL, dispersal$M, dispersal$N)
        if (file.exists(target_save)){
          next()
        }
        saveRDS(NULL, target_save)
        print("Init all potential distributions")
        distributoins<-list()
        for (year in predict_range){
          #print(year)
          env_item<-readRDS(sprintf("%s/%s_%d.rda", enm_folder, layer_item$LABEL, year))
          distributoins[[as.character(year)]]<-env_item
        }
        for (year in predict_range){
          print(paste(item$sp, layer_item$LABEL, year, 
                      j, ":", nrow(layer_df), "/",
                      i, ":", nrow(df_list), "/",
                      k, ":", nrow(dispersals), "/",
                      model))
          prev_dis<-distributoins[[as.character(year)]]
          if (nrow(prev_dis)==0){
            next()
          }
          
          if (nrow(prev_dis)>0){
            prev_dis$M<-dispersal$M
            prev_dis$N<-dispersal$N
            prev_dis$N_Step<-0
            prev_dis$YEAR<-year
            dispersal_log<-bind(dispersal_log, prev_dis)
          }
        }
        saveRDS(dispersal_log, target_save)

      }
      
    }
  }
  
}
if (F){
  df<-readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Niche_Models_Mean_GCM/Mammals/Carollia_benkeithi/dispersal/UKESM1_SSP585_-1_1.rda")
  df2<-readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Niche_Models_Mean_GCM/Mammals/Carollia_benkeithi/dispersal/UKESM1_SSP585_0_1.rda")
  dim(df2)
}