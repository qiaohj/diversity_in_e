library(raster)
library(dplyr)
library(concaveman)
library(sf)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
j_index<-as.numeric(args[1])
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
no_na<-!is.na(values(mask))
if (is.na(j_index)){
  j_index<-1
}

source("commonFuns/functions.r")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

if (T){
  threshold<-as.numeric(args[2])
  if (is.na(threshold)){
    threshold<-5
  }
  
  smooth_path<-NULL
  plot.new()
  #for (j in c(nrow(layer_df):1)){
  for (j in c(j_index)){
    layer_item<-layer_df[j,]
    for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
      i=1
      #dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
      dispersals<-c(1:2)
      df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
      final_df<-NULL
      colors<-rainbow(length(2021:2100))
      #p<-ggplot()
      
      for (i in c(1:nrow(df_list))){
        print(paste(i, nrow(df_list), group, layer_item$LABEL))
        item<-df_list[i,]
        item$sp<-gsub(" ", "_", item$sp)
        if (item$area<=0){
          next()
        }
        target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
        
        if (threshold==5){
          target<-sprintf("%s/dispersal_%d", target_folder, threshold)
        }else{
          target<-sprintf("%s/dispersal", target_folder)
        }
        
        k=1
        
        dispersal<-dispersals[k]
        ttt<-sprintf("../../Objects/dispersal_path_%d/%s/%s_%s_%d.rda", 
                     threshold, group, item$sp, layer_item$LABEL, dispersal)
        
        all_c<-readRDS(ttt)
        if (is.null(all_c)){
          next()
        }
        
        #all_c<-all_c%>%dplyr::filter(YEAR>=2020)
        if (nrow(all_c)<=1){
          next()
        }
        
        
        #sm_path<-data.frame(xspline(all_c$gravity_x, all_c$gravity_y, shape=1, draw=F))
        for (ii in unique(all_c$continent_i)){
          #keyyears<-c(2020, 2040, 2060, 2080, 2100)
          all_c_continent<-all_c%>%filter(continent_i==ii)
          if (nrow(all_c_continent)<=1){
            next()
          }
          keyyears<-unique(all_c_continent$YEAR)
          keyyears<-unique(round(quantile(keyyears, seq(0, 1, 0.2))))
          all_c_few<-all_c_continent%>%dplyr::filter(YEAR %in% keyyears)
          if (nrow(all_c_few)<=1){
            next()
          }
          sm_path<-data.frame(xspline(all_c_few$gravity_x, all_c_few$gravity_y, shape=1, draw=F, repEnds=T))
          sm_path$YEAR<-seq(keyyears[1], keyyears[length(keyyears)], by=(keyyears[length(keyyears)]-keyyears[1])/(nrow(sm_path)-1))
          sm_path$group<-group
          sm_path$sp<-item$sp
          if (max(keyyears)==2100){
            sm_path$survive<-"SURVIVE"
          }else{
            sm_path$survive<-"EXTINCT"
          }
          sm_path$continent_i<-ii
          #geom_path(aes(x=gravity_x, y=gravity_y, color=YEAR))+
          #p<-p+geom_path(data=sm_path, aes(x=x, y=y, color=YEAR))
          smooth_path<-bind(smooth_path, sm_path)
        }
      }
    }
    saveRDS(smooth_path, sprintf("../../Figures/Top_Figure_5/smooth_path_%s.rda", layer_item$LABEL))
  }
}

