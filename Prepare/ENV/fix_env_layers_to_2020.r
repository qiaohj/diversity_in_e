setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(data.table)
source("commonFuns/functions.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(1850:2020)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

mask<-data.frame(rasterToPoints(raster("../../Raster/mask_index.tif")))
start_env_layers<-NULL
for (i in c(1:nrow(start_layer_df))){
  print(paste("Init layer list:", i, nrow(start_layer_df)))
  item<-start_layer_df[i,]
  if ((item$Y<2015)&(item$SSP!="SSP119")){
    next()
  }
  pr<-sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum")
  tasmax<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max")
  #tasmin<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min")
  #layers<-stack(c(pr, tasmax, tasmin))
  layers<-stack(c(pr, tasmax))
  #names(layers)<-c("PR", "TEMP_MAX", "TEMP_MIN")
  names(layers)<-c("PR", "TEMP_MAX")
  layers_table<-mask
  #layers_table[,c("PR", "TEMP_MAX", "TEMP_MIN")]<-raster::extract(layers, mask[, c("x", "y")])
  layers_table[,c("PR", "TEMP_MAX")]<-raster::extract(layers, mask[, c("x", "y")])
  layers_table$year<-item$Y
  layers_table$GCM<-item$GCM
  layers_table$SSP<-item$SSP
  
  start_env_layers<-bind(start_env_layers, layers_table)
}

start_env_layers_se<-start_env_layers%>%dplyr::group_by(x, y, mask_index, year, GCM)%>%
  dplyr::summarise(mean_PR=mean(PR),
                   mean_TEMP=mean(TEMP_MAX))
dim(start_env_layers_se)
start_env_layers_se[is.na(start_env_layers_se)]
dim(start_env_layers)
colnames(start_env_layers_se)[c(6,7)]<-c("PR", "TEMP")
start_env_layers_se<-data.table(start_env_layers_se)
saveRDS(start_env_layers_se, "../../Objects/stacked_layers_1850_2020_df.rda")

var_se<-start_env_layers_se%>%dplyr::group_by(year, GCM)%>%
  dplyr::summarise(mean_temp=mean(TEMP),
                   mean_pr=mean(PR))
library(ggplot2)
ggplot(var_se)+geom_line(aes(x=year, y=mean_temp, color=factor(GCM)))


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=predict_range)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, layer_df$Y, sep="_")
mask<-data.frame(rasterToPoints(raster("../../Raster/mask_index.tif")))
env_layers<-list()
for (i in c(1:nrow(layer_df))){
  print(paste("Init layer list:", i, nrow(layer_df)))
  item<-layer_df[i,]
  pr<-sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum")
  tasmax<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max")
  #tasmin<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min")
  #layers<-stack(c(pr, tasmax, tasmin))
  layers<-stack(c(pr, tasmax))
  #names(layers)<-c("PR", "TEMP_MAX", "TEMP_MIN")
  names(layers)<-c("PR", "TEMP")
  layers_table<-mask
  layers_table[,c("PR", "TEMP")]<-raster::extract(layers, mask[, c("x", "y")])
  env_layers[[item$LABEL]]<-layers_table
}
saveRDS(env_layers, "../../Objects/stacked_layers_2021_2100.rda")

