setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(data.table)
source("commonFuns/functions.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")

start_range<-c(1850:2020)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=start_range)

var_tamplate<-"/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d/%s_eck4.tif"
if (F){
  mask<-raster("../../Raster/mask_10km.tif")
  na<-is.na(values(mask))
  values(mask)<-c(1:length(values(mask)))
  values(mask)[na]<-NA
  plot(mask)
  writeRaster(mask, "../../Raster/mask_10km.tif", overwrite=T)
  
  mask<-raster("../../Raster/mask_index.tif")
  na<-is.na(values(mask))
  values(mask)<-c(1:length(values(mask)))
  values(mask)[na]<-NA
  plot(mask)
  writeRaster(mask, "../../Raster/mask_100km.tif", overwrite=T)
  
  mask<-raster("../../Raster/mask_100km.tif")
  continent<-raster("../../Raster/Continent_ect4.tif")
  
  na<-is.na(values(continent))
  
  values(mask)[na]<-NA
  plot(mask)
  writeRaster(mask, "../../Raster/mask_100km_plot.tif", overwrite=T)
  
}
mask<-data.frame(rasterToPoints(raster("../../Raster/mask_100km.tif")))
start_env_layers<-list()
i=1
for (i in c(1:nrow(start_layer_df))){
  print(paste("Init layer list:", i, nrow(start_layer_df)))
  item<-start_layer_df[i,]
  if ((item$Y<2015)&(item$SSP!="SSP119")){
    next()
  }
  bio1<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio1")
  bio5<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio5")
  bio6<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio6")
  bio12<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio12")
  bio13<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio13")
  bio14<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio14")
  
  layers<-stack(c(bio1, bio5, bio6, bio12, bio13, bio14))
  names(layers)<-VARs
  layers_table<-mask
  layers_table[,VARs]<-raster::extract(layers, mask[, c("x", "y")])
  layers_table$year<-item$Y
  layers_table$GCM<-item$GCM
  layers_table$SSP<-item$SSP
  start_env_layers[[length(start_env_layers)+1]]<-layers_table
}
start_env_layers<-rbindlist(start_env_layers)

start_env_layers_se<-start_env_layers[, .(bio1=mean(bio1),
                                          bio5=mean(bio5),
                                          bio6=mean(bio6),
                                          bio12=mean(bio12),
                                          bio13=mean(bio13),
                                          bio14=mean(bio14)), by=list(x, y, mask_100km, year, GCM)]
dim(start_env_layers_se)

dim(start_env_layers)
#
start_env_layers_se<-data.table(start_env_layers_se)
if (F){
  start_env_layers_se<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
  colnames(start_env_layers_se)[c(6,7,8)]<-c("bio1", "bio5", "bio6")
}
setkey(start_env_layers_se, mask_100km)

saveRDS(start_env_layers_se, "../../Objects/stacked_layers_1850_2020_df_100km.rda")

var_se<-start_env_layers_se%>%dplyr::group_by(year, GCM)%>%
  dplyr::summarise(mean_temp_max=mean(TEMP_MAX),
                   mean_temp_min=mean(TEMP_MIN),
                   mean_pr=mean(PR))
library(ggplot2)
#ggplot(var_se)+geom_line(aes(x=year, y=mean_temp_max, color=factor(GCM)))


setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(data.table)
source("commonFuns/functions.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=predict_range)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, layer_df$Y, sep="_")
mask<-data.frame(rasterToPoints(raster("../../Raster/mask_100km.tif")))
env_layers<-list()
var_tamplate<-"/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d/%s_eck4.tif"
i=645
for (i in c(1:nrow(layer_df))){
  print(paste("Init layer list:", i, nrow(layer_df)))
  item<-layer_df[i,]
  bio1<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio1")
  bio5<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio5")
  bio6<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio6")
  bio12<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio12")
  bio13<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio13")
  bio14<-sprintf(var_tamplate, item$GCM, item$SSP, item$Y, "bio14")
  
  layers<-stack(c(bio1, bio5, bio6, bio12, bio13, bio14))
  
  names(layers)<-VARs
  layers_table<-mask
  layers_table[,VARs]<-raster::extract(layers, mask[, c("x", "y")])
  layers_table<-data.table(layers_table)
  setkey(layers_table, mask_100km)
  
  env_layers[[item$LABEL]]<-layers_table
}
#saveRDS(env_layers, "../../Objects/stacked_layers_2021_2100.rda")

for (i in names(env_layers)){
  print(i)
  labels<-strsplit(i, "_")[[1]]
  env_layers[[i]]$GCM<-labels[1]
  env_layers[[i]]$SSP<-labels[2]
  env_layers[[i]]$year<-year<-as.numeric(labels[3])
}
future_env_layers_full<-rbindlist(env_layers)
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
future_env_layers_list<-list()
i=1
for (i in c(1:nrow(layer_df))){
  print(paste(i, paste(nrow(layer_df))))
  item<-layer_df[i,]
  df_item<-future_env_layers_full[(GCM==item$GCM)&(SSP==item$SSP)]
  #setkey(df_item, mask_10km)
  setindex(df_item, mask_100km)
  setindex(df_item, year)
  future_env_layers_list[[item$LABEL]]<-df_item
}

saveRDS(future_env_layers_list, "../../Objects/stacked_layers_2021_2100_list_100km.rda")

if (F){
  setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
  
  library(data.table)
  start_env_layers_se<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
  setkey(start_env_layers_se, mask_100km)
  saveRDS(start_env_layers_se, "../../Objects/stacked_layers_1850_2020_df.rda")
}

if (F){
  setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
  
  library(data.table)
  future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")
  for (i in c(1:length(future_env_layers))){
    print(paste(i, paste(length(future_env_layers))))
    future_env_layers[[i]]<-data.table(future_env_layers[[i]])
    setkey(future_env_layers[[i]], mask_100km)
  }
  
  saveRDS(future_env_layers, "../../Objects/stacked_layers_2021_2100_xx.rda")
}

if (F){
  setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
  
  library(data.table)
  future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")
  i<-names(future_env_layers)[1]
  for (i in names(future_env_layers)){
    print(i)
    labels<-strsplit(i, "_")[[1]]
    future_env_layers[[i]]$GCM<-labels[1]
    future_env_layers[[i]]$SSP<-labels[2]
    future_env_layers[[i]]$year<-year<-as.numeric(labels[3])
  }
  future_env_layers_full<-rbindlist(future_env_layers)
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
  layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
  future_env_layers_list<-list()
  i=1
  for (i in c(1:nrow(layer_df))){
    print(paste(i, paste(nrow(layer_df))))
    item<-layer_df[i,]
    df_item<-future_env_layers_full[(GCM==item$GCM)&(SSP==item$SSP)]
    #setkey(df_item, mask_10km)
    setindex(df_item, mask_10km)
    setindex(df_item, year)
    future_env_layers_list[[item$LABEL]]<-df_item
  }
  
  saveRDS(future_env_layers_list, "../../Objects/stacked_layers_2021_2100_list.rda")
}


if (F){
  setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
  
  library(data.table)
  future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_100km.rda")
  i<-names(future_env_layers)[1]
  for (i in names(future_env_layers)){
    print(i)
    labels<-strsplit(i, "_")[[1]]
    future_env_layers[[i]]$GCM<-labels[1]
    future_env_layers[[i]]$SSP<-labels[2]
    future_env_layers[[i]]$year<-year<-as.numeric(labels[3])
  }
  future_env_layers_full<-rbindlist(future_env_layers)
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
  layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
  future_env_layers_list<-list()
  i=1
  for (i in c(1:nrow(layer_df))){
    print(paste(i, paste(nrow(layer_df))))
    item<-layer_df[i,]
    df_item<-future_env_layers_full[(GCM==item$GCM)&(SSP==item$SSP)]
    #setkey(df_item, mask_10km)
    setindex(df_item, mask_100km)
    setindex(df_item, year)
    future_env_layers_list[[item$LABEL]]<-df_item
  }
  
  saveRDS(future_env_layers_list, "../../Objects/stacked_layers_2021_2100_list_100km.rda")
}


#start_env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
#start_env_layers$TEMP_MAX<-start_env_layers$TEMP_MAX*10-273.16
#start_env_layers$TEMP_MIN<-start_env_layers$TEMP_MIN*10-273.16
#saveRDS(start_env_layers, "../../Objects/stacked_layers_1850_2020_df.rda")
