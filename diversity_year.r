library(raster)
#library(rgdal)
#library(rgeos)
#library(MASS)
#library(cluster)
library(dplyr)
library(ggplot2)
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
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, -1), N=c(rep(1,5), c(2:5), 2, 1))
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]


j=1
k=1
for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:nrow(dispersals))){
    layer$M<-dispersals[k, "M"]
    layer$N<-dispersals[k, "N"]
    target_folder<-sprintf("../../Objects/Diversity/%s_%d_%d", layer$LABEL, layer$M, layer$N)
    
    diversity_df<-readRDS(sprintf("%s/diversity_df.rda", target_folder))
  }
  
  
  
}