library(raster)
library(dplyr)
library(concaveman)
library(sf)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
alt<-raster("../../Raster/ALT/alt_eck4")


no_na<-!is.na(values(mask))

source("commonFuns/functions.r")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

if (F){
  
  j=6
  threshold=1
  for (j in c(1:nrow(layer_df))){
    for (threshold in c(1, 5)){
      print(paste(j, nrow(layer_df), threshold))
      layer_item<-layer_df[j,]
      raw_path<-readRDS(sprintf("../../Figures/Top_Figure_%d/raw_path_%s.rda", threshold, layer_item$LABEL))
      smooth_path<-data.frame(smooth_path)
      smooth_path$alt<-raster::extract(alt, smooth_path[, c("x", "y")])
      saveRDS(smooth_path, sprintf("../../Figures/Top_Figure_%d/smooth_path_%s_with_alt.rda", threshold, layer_item$LABEL))
    }
  }
}

smooth_path$YEAR<-round(smooth_path$YEAR)

smooth_path_min_max_year<-smooth_path%>%
  dplyr::group_by(group, sp, survive, continent_i)%>%
  dplyr::summarise(MAX_YEAR=max(YEAR),
                   MIN_YEAR=min(YEAR))

smooth_path_2<-inner_join(smooth_path_min_max_year, smooth_path, 
                         by=c("group", "sp", "survive", "continent_i", "MAX_YEAR"="YEAR"))
colnames(smooth_path_2)[7:9]<-c("end_x", "end_y", "end_alt")
smooth_path_2<-inner_join(smooth_path_min_max_year, smooth_path, 
                          by=c("group", "sp", "survive", "continent_i", "MAX_YEAR"="YEAR"))


smooth_path[which(smooth_path$sp=="Abavorana_luctuosa"),]

smooth_path_se<-smooth_path%>%dplyr::group_by(YEAR, group, survive)%>%
  dplyr::summarise(mean_alt=mean(alt, na.rm=T),
                   sd_alt=sd(alt, na.rm=T))

library(ggplot2)
ggplot(smooth_path_se)+
  geom_line(aes(x=YEAR, y=mean_alt, color=factor(group)))+
  facet_wrap(~survive, scale="free", nrow=2)
