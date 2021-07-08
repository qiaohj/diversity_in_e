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


GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
l_i=1
exposure=0
GCM<-GCMs[1]
SSP<-SSPs[1]
mask_p<-data.table(rasterToPoints(raster("../../Raster/mask_100km.tif")))
for (SSP_i in c(1:length(SSPs))){
  SSP<-SSPs[SSP_i]
  for (exposure in c(0, 5)){
    df_all<-list()
    for (GCM_i in c(1:length(GCMs))){
      GCM<-GCMs[GCM_i]
      target_folder<-sprintf("../../Objects/density_based_pathway/%s_%s_exposure_%d.rda", 
                             GCM, SSP, exposure)
      df<-readRDS(target_folder)
      df<-rbindlist(df)
      if (F){
        ggplot(df[["2100"]])+geom_tile(aes(x=x, y=y, fill=count))
        ggplot(df_all_all_year)+geom_tile(aes(x=x, y=y, fill=mean_count))
      }
      df_all[[GCM]]<-df
    }
    df_all<-rbindlist(df_all)
    df_all_all_year<-df_all[, .(mean_count=mean(count)), by=list(x, y, mask_100km, SSP, exposure)]
    df_all_all_year<-df_all_all_year[mask_100km %in% mask_p$mask_100km]
  }
}
