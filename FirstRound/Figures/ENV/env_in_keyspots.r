library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(scales)
library(rgdal)
library(rgeos)
library(Rmisc)
library(geometry)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")

mask<-raster("../../Raster/mask_index.tif")

mask_p<-data.frame(rasterToPoints(mask))

if (T){
  keyspots<-data.table(name=c("Sahara", "Central Australia"), 
                       x=c(1222400, 12290769), y=c(2026400, -3315500))
  cord.eck4 = SpatialPoints(keyspots[, c("x", "y")], proj4string=crs(mask))
  plot(mask)
  plot(cord.eck4, add=T)
  
  pc100km <- gBuffer( cord.eck4, width=100*10e3, byid=TRUE )
  plot(pc100km, add=T, col="red")
  pc100km <- SpatialPolygonsDataFrame(pc100km, data=keyspots)
  
  mask_p_sp<-SpatialPointsDataFrame(mask_p[, c("x", "y")], 
                                    data=mask_p, proj4string=crs(mask))
  
  mask_p_sp$name<-NA
  for (i in c(1:nrow(pc100km))){
    overed<-!is.na(sp::over(mask_p_sp, pc100km[i,])[,1])
    mask_p_sp[overed, "name"]<-pc100km[i,]$name
  }
  table(mask_p_sp$name)
  
  i=1
  layer_all<-NULL
  for (i in c(1:length(env_layers))){
    print(paste(i, length(env_layers)))
    layer<-env_layers[[i]]
    layer_x<-left_join(layer, mask_p_sp@data[, c("mask_index", "name")], by="mask_index")
    layer_x[is.na(layer_x$name), "name"]<-"Others"
    layer_x_se<-layer_x%>%dplyr::group_by(name)%>%
      dplyr::summarise(mean_temp_max=mean(TEMP_MAX),
                       mean_temp_min=mean(TEMP_MIN),
                       mean_pr=mean(PR))
    labels<-strsplit(names(env_layers)[i], "_")[[1]]
    layer_x_se$GCM<-labels[1]
    layer_x_se$SSP<-labels[2]
    layer_x_se$YEAR<-as.numeric(labels[3])
    layer_all<-bind_dplyr(layer_all, layer_x_se)
    #table(layer_x$name)
  }
  p<-ggplot(layer_all%>%dplyr::filter(SSP=="SSP585"))+geom_line(aes(x=YEAR, y=mean_pr, color=factor(name)))+
    facet_wrap(~GCM, nrow=3)
  
  p
  
  p<-ggplot(layer_all)+geom_line(aes(x=YEAR, y=mean_pr, color=factor(name)))+
    facet_grid(SSP~GCM)
  
  p
  
  
  p<-ggplot(layer_all%>%dplyr::filter(SSP=="SSP585"))+geom_line(aes(x=YEAR, y=mean_temp_max, color=factor(name)))+
    facet_wrap(~GCM, nrow=3)
  
  p
  p<-ggplot(layer_all%>%dplyr::filter(SSP=="SSP585"))+geom_line(aes(x=YEAR, y=mean_temp_min, color=factor(name)))+
    facet_wrap(~GCM, nrow=3)
  
  p
}
