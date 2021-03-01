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

mask<-raster("../../Raster/mask.tif")
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
mask<-raster("../../Raster/mask_index.tif")
r_continent<-raster("../../Raster/Continent_ect4.tif")
sp_i<-1
threshold<-5
l_i<-1
group<-"Birds"
width<-10
height<-6
persents<-c(1, 0.5, 0.4, 0.3, 0.2, 0.1)
all_info<-NULL
#for (l_i in c(1:nrow(layer_df))){
for (l_i in c(1, 3)){
  layer_item<-layer_df[l_i,]
  for (threshold in c(5)){
    
    for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      print(paste(group, threshold, layer_item$LABEL))
      smooth_path<-readRDS(sprintf("../../Objects_Full_species/cluster_based_pathway/merged/%s_%s_exposure_%d.rda",
                                   group, layer_item$LABEL, threshold))
      smooth_path$line_group<-paste(smooth_path$sp, smooth_path$path_group)
      
      smooth_path_se<-smooth_path[, .(max_year=max(YEAR),
                                      length=.N), by=list(sp, line_group)]
      #smooth_path_se_extant<-smooth_path_se[sp %in% sp_status[STATUS=="extant"]$sp]
      
      sp_N_Pathways<-unique(smooth_path_se[, c("sp", "line_group")])
      sp_N_Pathways_se<-sp_N_Pathways[, .(N_pathways=.N), by=sp]
      info_item<-layer_item
      info_item$mean_N_cluster<-mean(sp_N_Pathways_se$N_pathways)
      info_item$sd_N_cluster<-sd(sp_N_Pathways_se$N_pathways)
      info_item$group<-group
      info_item$exposure<-ifelse(threshold==1, " no exposure", "5-year exposure")
      all_info<-bind(all_info, info_item)
      for (persent in persents){
        print(paste(group, threshold, layer_item$LABEL, persent))
        target_rda<-sprintf("../../Objects_Full_species/cluster_based_pathway/merged/%s_%s_exposure_%d_sub_%d.rda",
                          group, layer_item$LABEL, threshold, persent * 100)
        if (file.exists(target_rda)){
          smooth_path_sub<-readRDS(target_rda)
          
        }else{
          all_line_group<-unique(smooth_path$line_group)
          all_line_group_sample<-all_line_group[sample(length(all_line_group), 
                                                       round(length(all_line_group)*persent))]
          smooth_path_sub<-smooth_path[line_group %in% all_line_group_sample]
          saveRDS(smooth_path_sub, target_rda)
        }
        
        
        smooth_path_sub$alpha<-((smooth_path_sub$YEAR-2020)/80)^5
        p<-p_bak+geom_path(data=smooth_path_sub, aes(x=x, y=y, alpha=alpha, color=group,
                                                     group=line_group))+
          scale_alpha_continuous()+
          scale_color_manual(values = color_groups)
        ggsave(p, filename=
                 sprintf("../../Figures_Full_species/cluster_based_pathway/%s_%s_exposure_%d_sub_%d.png", 
                         group, layer_item$LABEL, threshold, persent * 100), 
               width=width, height = height)
      }
      if (F){
        sp_status<-smooth_path[, .(max_year=max(YEAR)), by=list(sp)]
        sp_status$STATUS<-ifelse(sp_status$max_year==2100, "extant", "extinct")
        table(sp_status$STATUS)
        smooth_path_se<-smooth_path[, .(max_year=max(YEAR),
                                        length=.N), by=list(sp, line_group)]
        smooth_path_se_extant<-smooth_path_se[sp %in% sp_status[STATUS=="extant"]$sp]
        
        sp_N_Pathways<-unique(smooth_path_se[, c("sp", "line_group")])
        sp_N_Pathways_se<-sp_N_Pathways[, .(N_pathways=.N), by=sp]
        mean<-mean(sp_N_Pathways_se$N_pathways)
        sd<-sd(sp_N_Pathways_se$N_pathways)
      }
      
      
    }
  }
}
