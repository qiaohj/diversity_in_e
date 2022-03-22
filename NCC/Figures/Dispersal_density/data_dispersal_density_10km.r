library(ceramic)
library(raster)
library(data.table)
library(igraph)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(tidyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
rm(list=ls())
source("commonFuns/colors.r")
source("commonFuns/functions.r")

setDTthreads(threads=1)
print(sprintf("%d CPUs are using", getDTthreads()))
points_10km<-readRDS("../../Raster/points_10km.rda")
colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-data.table(rasterToPoints(mask_100km))
colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
points_10km<-merge(points_10km, points_100km, by="mask_100km")

GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1
df_sp_list<-list()
group<-"Mammals"
for (group in c("Mammals", "Birds")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
  df_list$group<-group
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list, fill=T)


#df_sp_list<-df_sp_list[area>ttt]
df_sp_list$SP<-gsub(" ", "_", df_sp_list$SP)
df_sp_list<-df_sp_list[sample(nrow(df_sp_list), nrow(df_sp_list)),]
j=1
mask<-raster("../../Raster/mask_100km.tif")
sp_i<-412
exposure<-5
l_i<-1
mask_p<-data.table(rasterToPoints(mask))
for (l_i in c(1:nrow(layer_df))){
  layer_item<-layer_df[l_i,]
  
  for (exposure in c(0, 5)){
    ppp<-mask_p
    ppp$GCM<-layer_item$GCM
    ppp$SSP<-layer_item$SSP
    ppp$exposure<-exposure
    target_folder<-sprintf("../../Objects/density_based_pathway_10km/%s_exposure_%d.rda", layer_item$LABEL, exposure)
    
    if (file.exists(target_folder)){
      print("skip")
      next()
    }
    saveRDS(NULL, target_folder)
    p_list<-list()
    for (year in c(2021:2100)){
      t_p<-ppp
      t_p$year<-year
      t_p$count<-0
      p_list[[as.character(year)]]<-t_p
    }
    for (sp_i in c(1:nrow(df_sp_list))){
      
      sp_test<-df_sp_list[sp_i,]
      print(sprintf("combination:%d/%d %s, exposure:%d, sp:%d/%d, %s, %s", 
                    l_i, nrow(layer_df), layer_item$LABEL, exposure, sp_i, nrow(df_sp_list),
                    sp_test$SP, sp_test$group))
      
      source_folder<-sprintf("../../Objects/Dispersal/%s/%s", sp_test$group, sp_test$SP)
      f<-sprintf("%s/initial_disp_10km_exposure_%d_dispersal_1.rda", 
                 source_folder,  exposure)
      res<-"10km"
      dis_details_f<-sprintf("%s/%s_%d_dispersal_1_10km.rda", 
                             source_folder, layer_item$LABEL, exposure)
      if ((!file.exists(f))|(!file.exists(dis_details_f))){
        f<-sprintf("%s/initial_disp_exposure_%d_dispersal_1.rda", 
                   source_folder,  exposure)
        dis_details_f<-sprintf("%s/%s_%d_dispersal_1.rda", 
                               source_folder, layer_item$LABEL, exposure)
        res<-"100km"
        
      }
      if (!file.exists(f)){
  
        next()
      }
      start_dis<-readRDS(f)
      
      dis_details<-readRDS(dis_details_f)
      dis_details<-rbindlist(dis_details)
      
      if (is.null(dis_details)){
        next()
      }
      if (nrow(dis_details)==0){
        next()
      }
      start_dis$YEAR<-2020
      dis_details<-rbindlist(list(dis_details, start_dis), fill=T)
      if (res=="10km"){
        dis_details<-merge(dis_details, points_10km, by.x="mask_10km", by.y="mask_10km")
        if (nrow(dis_details)==0){
          next()
        }
        
        dis_details<-dis_details[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR")]
        colnames(dis_details)[1:3]<-c("x", "y", "mask_100km")
      }
      year=2021
      for (year in c(2021:2100)){
        item_1<-dis_details[YEAR==(year-1)]
        item_2<-dis_details[YEAR==year]
        if ((nrow(item_1)==0)|nrow(item_2)==0){
          break()
        }
        new_item<-item_2[!(mask_100km %in% item_1$mask_100km)]
        
        p_list[[as.character(year)]][mask_100km %in% new_item$mask_100km]$count<-
          p_list[[as.character(year)]][mask_100km %in% new_item$mask_100km]$count+1
      }
      
    }
    saveRDS(p_list, target_folder)
  }
}
ggplot(p_list[["2100"]])+geom_tile(aes(x=x, y=y, fill=count))
