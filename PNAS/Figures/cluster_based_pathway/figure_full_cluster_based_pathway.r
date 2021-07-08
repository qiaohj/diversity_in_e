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

mask<-raster("../../Raster/mask_10km.tif")
mask_p<-data.frame(rasterToPoints(mask))

p_bak<-ggplot() + 
  geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
  map_theme

GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1



j=1
mask<-raster("../../Raster/mask_index_100km.tif")
r_continent<-raster("../../Raster/Continent_ect4.tif")
sp_i<-1
exposrue<-5
l_i<-1
group<-"Birds"
width<-13
height<-6

persents<-c(1, 0.5, 0.4, 0.3, 0.2, 0.1)
all_info<-NULL
persent<-0.2
#for (l_i in c(1:nrow(layer_df))){
for (l_i in c(1)){
  layer_item<-layer_df[l_i,]
  for (exposrue in c(5)){
    for (persent in persents){
      full_pathways<-list()
      for (group in c("Birds", "Mammals")){
        print(paste(group, exposrue, layer_item$LABEL, persent))
        target_rda<-sprintf("../../Objects/cluster_based_pathway/merged/%s_%s_exposure_%d_sub_%d.rda",
                            group, layer_item$LABEL, exposrue, persent * 100)
        smooth_path_sub<-readRDS(target_rda)
        
        
        smooth_path_sub$alpha<-((smooth_path_sub$YEAR-2020)/80)^5
        full_pathways[[group]]<-smooth_path_sub
      }
      full_pathways<-rbindlist(full_pathways)
      
      saveRDS(full_pathways, sprintf("../../Objects/cluster_based_pathway/merged/%s_%s_exposure_%d_sub_%d.rda",
                                     "ALL", layer_item$LABEL, exposrue, persent * 100))
      p<-p_bak+geom_path(data=full_pathways, aes(x=x, y=y, alpha=alpha, color=group,
                                                   group=line_group))+
        scale_alpha_continuous()+
        scale_color_manual(values = color_groups)
      ggsave(p, filename=
               sprintf("../../Figures/cluster_based_pathway/%s_%s_exposure_%d_sub_%d.png", 
                       "ALL", layer_item$LABEL, exposrue, persent * 100), 
             width=width, height = height)
      
    }
  }
}
