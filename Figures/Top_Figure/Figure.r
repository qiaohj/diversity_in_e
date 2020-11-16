

library(rgl)
library(ceramic)
library(anglr)
library(ggnewscale)
library(ggplot2)
library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
alt<-raster("../../Raster/ALT/alt_eck4_high_res.tif")
source("colors.r")
p2<-data.frame(rasterToPoints(alt))
p2[which(p2$x<=-12103059), "alt_eck4_high_res"]<-NA
p2[which((p2$x>12912000)&(p2$y>5000000)), "alt_eck4_high_res"]<-NA
values(alt)[!is.na(values(alt))]<-p2$alt_eck4_high_res

slope = terrain(alt, opt='slope')
aspect = terrain(alt, opt='aspect')
hill = hillShade(slope, aspect)
dem_spdf <- as(alt, "SpatialPixelsDataFrame")
dem_spdf <- as.data.frame(dem_spdf)
colnames(dem_spdf) <- c("value", "x", "y")

hill_spdf <- as(hill, "SpatialPixelsDataFrame")
hill_spdf <- as.data.frame(hill_spdf)
colnames(hill_spdf) <- c("value", "x", "y")

p_bak<-ggplot() + 
  geom_tile(data = hill_spdf, aes(x = x, y = y, fill = value)) + 
  scale_fill_gradient(low = "black", high = "white") + 
  new_scale_fill() + 
  geom_tile(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.4) + 
  scale_fill_gradientn(colours = rev(terrain.colors(10))) + 
  map_theme

args = commandArgs(trailingOnly=TRUE)
j_index<-as.numeric(args[1])


if (is.na(j_index)){
  j_index<-1
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

#for (j in c(1:nrow(layer_df))){
for (j in c(j_index)){
  layer_item<-layer_df[j,]
  smooth_path<-readRDS(sprintf("../../Figures/Top_Figure/smooth_path_%s.rda", layer_item$LABEL))
  print(sprintf("../../Figures/Top_Figure/smooth_path_%s.rda", layer_item$LABEL))

  smooth_path$line_group<-paste(smooth_path$sp, smooth_path$continent_i)
  smooth_path$alpha<-((smooth_path$YEAR-2014)/86)^5
  p<-p_bak+geom_path(data=smooth_path, aes(x=x, y=y, alpha=alpha, color=group,
                                           group=line_group))+
    scale_alpha_continuous()+
    scale_color_manual(values = color_groups)
  
  width<-10
  height<-6
  ggsave(p, filename=sprintf("../../Figures/Top_Figure/Top_Figure_ALL_%s.png", layer_item$LABEL), width=width, height = height)
  ggsave(p, filename=sprintf("../../Figures/Top_Figure/Top_Figure_ALL_%s.pdf", layer_item$LABEL), width=width, height = height)
  
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    print(g)
    p<-p_bak+geom_path(data=smooth_path%>%dplyr::filter(group==g), 
                       aes(x=x, y=y, alpha=alpha, color=group, group=line_group))+
      scale_alpha_continuous()+
      scale_color_manual(values = color_groups)
    
    
    ggsave(p, filename=sprintf("../../Figures/Top_Figure/Top_Figure_%s_%s.png", g, layer_item$LABEL), width=width, height = height)
    ggsave(p, filename=sprintf("../../Figures/Top_Figure/Top_Figure_%s_%s.pdf", g, layer_item$LABEL), width=width, height = height)
  }
  
}
