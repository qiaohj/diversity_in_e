library(raster)
library(dplyr)
library(sp)
library(rgdal)
library(Rmisc)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
SSPs<-c("SSP119", "SSP245", "SSP585")
GCMs<-c("EC-Earth3-Veg", "UKESM1", "MRI-ESM2-0")
Groups<-c("Amphbians", "Birds", "Mammals", "Reptiles")
dispersals<-c(0:1)
dis_i<-1
args = commandArgs(trailingOnly=TRUE)
g<-args[1]

SSP_i<-"SSP119"
GCM_i<-GCMs[1]
keyspots<-readOGR("../../Shape/pc100km", "pc100km") 
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
no_na<-!is.na(values(mask))
proj4string(keyspots)<-crs(mask)
y=2020
div_i<-"species.richness"
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-5
}
diectory<-sprintf("Diversity_%d", threshold)
dis_i<-1
#for (g in Groups){
result<-NULL
for (SSP_i in SSPs){
  for (dis_i in c(1:length(dispersals))){
    dis<-dispersals[dis_i]
    for (GCM_i in GCMs){
      df<-readRDS(sprintf("../../Objects/%s/%s/%s_%s_%d/indices_df.rda",
                          diectory, g, GCM_i, SSP_i, dis))
      for (y in c(2020:2100)){
        print(paste(g, SSP_i, dis, GCM_i, y, threshold))
        for (div_i in names(df[[as.character(y)]])){
          df_item<-df[[as.character(y)]][[div_i]]
          dir.create(sprintf("../../Figures/%s/%s/%s", diectory, g, div_i), showWarnings = F, recursive = T)
          dir.create(sprintf("../../Figures/%s/%s/%s/%s_%s_%d", diectory, g, div_i, GCM_i, SSP_i, dis), showWarnings = F, recursive = T)
          layer<-mask
          layer_p<-mask_p
          layer_p<-left_join(layer_p, df_item, by=c("mask_index"="index"))
          values(layer)[no_na]<-layer_p$metric
          if (y==2020){
            layer_base<-layer
          }
          layer_differ<-layer-layer_base
          writeRaster(layer, sprintf("../../Figures/%s/%s/%s/%s_%s_%d/%d.tif", diectory, g, div_i, GCM_i, SSP_i, dis, y), overwrite=T)
          writeRaster(layer_differ, sprintf("../../Figures/%s/%s/%s/%s_%s_%d/%d_differ.tif", diectory, g, div_i, GCM_i, SSP_i, dis, y), overwrite=T)
          over_points_result_se<-df_item%>%dplyr::summarise(mean=mean(metric, na.rm=T), 
                                                            sd=sd(metric),
                                                            CI=CI(metric)[1]-CI(metric)[2])
          over_points_result_se$name<-"Global"
          over_points_result_se$type<-div_i
          over_points_result_se$year<-y
          over_points_result_se$GCM<-GCM_i
          over_points_result_se$M<-dis
          over_points_result_se$SSP<-SSP_i
          over_points_result_se$group<-g
          result<-bind_dplyr(result, over_points_result_se)
          sp_exposure<-SpatialPoints(df_item[, c("x", "y")], proj4string=crs(mask))
          for (i in c(1:nrow(keyspots))){
            keyspot<-keyspots[i,]
            over_points<-sp::over(sp_exposure, keyspots[which(keyspots$name==keyspot$name),])
            over_points_result<-df_item[!is.na(over_points$name),]
            over_points_result_se<-over_points_result%>%dplyr::summarise(mean=mean(metric), 
                                                                         sd=sd(metric),
                                                                         CI=CI(metric)[1]-CI(metric)[2])
            over_points_result_se$name<-keyspot$name
            over_points_result_se$type<-div_i
            over_points_result_se$year<-y
            over_points_result_se$GCM<-GCM_i
            over_points_result_se$M<-dis
            over_points_result_se$SSP<-SSP_i
            over_points_result_se$group<-g
            result<-bind_dplyr(result, over_points_result_se)
          }
        }
      }
    }
  }
}
#}
saveRDS(result, sprintf("../../Figures/%s/%s/%s.rda", diectory, g, g))
