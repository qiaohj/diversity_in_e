

library(rgl)
library(ceramic)
library(anglr)
library(ggnewscale)
library(ggplot2)
library(raster)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
alt<-raster("../../Raster/ALT/alt_eck4_high_res.tif")
mask<-raster("../../Raster/mask.tif")
mask_p<-data.frame(rasterToPoints(mask))
source("commonFuns/colors.r")
#p2<-data.frame(rasterToPoints(alt))
#p2[which(p2$x<=-12103059), "alt_eck4_high_res"]<-NA
#p2[which((p2$x>12912000)&(p2$y>5000000)), "alt_eck4_high_res"]<-NA
#values(alt)[!is.na(values(alt))]<-p2$alt_eck4_high_res
if (F){
  slope = terrain(alt, opt='slope')
  aspect = terrain(alt, opt='aspect')
  hill = hillShade(slope, aspect)
  dem_spdf <- as(alt, "SpatialPixelsDataFrame")
  dem_spdf <- as.data.frame(dem_spdf)
  colnames(dem_spdf) <- c("value", "x", "y")
  
  hill_spdf <- as(hill, "SpatialPixelsDataFrame")
  hill_spdf <- as.data.frame(hill_spdf)
  colnames(hill_spdf) <- c("value", "x", "y")
}
p_bak<-ggplot() + 
  geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
  map_theme

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1
df_sp_list<-list()
for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list)

ttt<-2
df_sp_list<-df_sp_list[area>ttt]
df_sp_list$sp<-gsub(" ", "_", df_sp_list$sp)
j=9
euc.dist <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2) ^ 2+(y1-y2)^2)
}
threshold=5
for (j in c(1:nrow(layer_df))){
  for (threshold in c(1, 5)){
    layer_item<-layer_df[j,]
    smooth_path<-readRDS(sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s.rda", 
                                 threshold, layer_item$LABEL))
    print(sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s.rda", threshold, layer_item$LABEL))
    
    smooth_path<-smooth_path[sp %in% df_sp_list$sp]
    smooth_path$line_group<-paste(smooth_path$sp, smooth_path$continent_i)
    smooth_path$index<-raster::extract(mask, smooth_path[, c("x", "y")])
    smooth_path<-smooth_path[!is.na(index)]
    smooth_path$pre_x<--1
    smooth_path[2:nrow(smooth_path), "pre_x"]<-smooth_path[1:(nrow(smooth_path)-1), "x"]
    smooth_path$pre_y<--1
    smooth_path[2:nrow(smooth_path), "pre_y"]<-smooth_path[1:(nrow(smooth_path)-1), "y"]
    smooth_path$dist<-euc.dist(smooth_path$x, smooth_path$y, smooth_path$pre_x, smooth_path$pre_y)
    smooth_path[1, "dist"]<-0
    smooth_path[YEAR==2020]$dist<-0
    smooth_path<-smooth_path[dist<100000]
    
    smooth_path$alpha<-((smooth_path$YEAR-2020)/80)^5
    p<-p_bak+geom_path(data=smooth_path, aes(x=x, y=y, alpha=alpha, color=group,
                                             group=line_group))+
      scale_alpha_continuous()+
      scale_color_manual(values = color_groups)
    
    width<-10
    height<-6
    ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_%d/Top_Figure_ALL_%s_ttt_%d.png", 
                               threshold, layer_item$LABEL, ttt), width=width, height = height)
    ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_%d/Top_Figure_ALL_%s_ttt_%d.pdf", 
                               threshold, layer_item$LABEL, ttt), width=width, height = height)
    
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      print(g)
      p<-p_bak+geom_path(data=smooth_path%>%dplyr::filter(group==g), 
                         aes(x=x, y=y, alpha=alpha, color=group, group=line_group))+
        scale_alpha_continuous()+
        scale_color_manual(values = color_groups)
      
      
      ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_%d/Top_Figure_%s_%s_ttt_%d.png",
                                 threshold, g, layer_item$LABEL, ttt), width=width, height = height)
      ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_%d/Top_Figure_%s_%s_ttt_%d.pdf", 
                                 threshold, g, layer_item$LABEL, ttt), width=width, height = height)
    }
  }
  
}
