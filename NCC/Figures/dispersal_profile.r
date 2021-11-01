library(ggplot2)
library(data.table)
library(sf)
library(raster)
source("commonFuns/colors.r")
sp<-"Poicephalus_rufiventris"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
df<-readRDS(sprintf("../../Objects/Birds/%s/UKESM1_SSP585_0.rda", sp))
all_df<-rbindlist(df)
x_range<-range(all_df$x)
x_range[1]<-x_range[1]-100000
x_range[2]<-x_range[2]+100000

y_range<-range(all_df$y)
y_range[1]<-y_range[1]-100000
y_range[2]<-y_range[2]+100000

init_dis<-readRDS(sprintf("../../Objects/Birds/%s/initial_disp_exposure_0.rda", sp))

shape<-readRDS(sprintf("../../Objects/IUCN_Distribution/Birds/RAW/%s.rda", sp))
mask<-raster("../../Raster/mask_100km.tif")
mask_p<-data.frame(rasterToPoints(mask))


p_bak<-ggplot() + 
  geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
  geom_tile(data=init_dis, aes(x=x, y=y), fill=colors_blue[3])+
  geom_sf(data=shape, color=colors_red[7], fill=NA)+
  xlim(x_range)+ylim(y_range)+
  ggtitle(2000)+
  map_theme
ggsave(p_bak, filename=sprintf("../../Figures/Example/dispersal/%d.png", 2020))

YYYY=2022
future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")
fit<-readRDS(sprintf("../../Objects/Birds/%s/fit.rda", sp))
env_item<-future_env_layers[["UKESM1_SSP585"]]

for (YYYY in c(2021:2100)){
  print(YYYY)
  
  suitable_item<-env_item[between(bio1, fit$range_bio1_sd_min, fit$range_bio1_sd_max)&
               between(bio5, fit$range_bio5_sd_min, fit$range_bio5_sd_max)&
               between(bio6, fit$range_bio6_sd_min, fit$range_bio6_sd_max)&
               between(bio12, fit$range_bio12_sd_min, fit$range_bio12_sd_max)&
               between(bio13, fit$range_bio13_sd_min, fit$range_bio13_sd_max)&
               between(bio14, fit$range_bio14_sd_min, fit$range_bio14_sd_max)]
  suitable_item<-suitable_item[year==YYYY]
  df_item<-df[[as.character(YYYY)]]
  dispersaled_points<-df_item[accumulative_disp>0]
  pts     <- sf::st_as_sf(dispersaled_points, coords = c("x", "y"), remove = F, crs=crs(mask))
  pts_buf <- sf::st_buffer(pts, dispersaled_points$accumulative_disp)
  pts_buf_union<-st_union(pts_buf, shape)
  
  p_bak<-ggplot() + 
    geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
    geom_tile(data=suitable_item, aes(x=x, y=y), fill=colors_green[3], alpha=0.5)+
    geom_tile(data=df_item, aes(x=x, y=y), fill=colors_blue[3])+
    geom_tile(data=df_item[accumulative_disp==0], aes(x=x, y=y), fill=colors_red[3])+
    geom_sf(data=pts_buf, color=colors_red[7], fill=NA)+
    xlim(x_range)+ylim(y_range)+
    ggtitle(YYYY)+
    map_theme
  ggsave(p_bak, filename=sprintf("../../Figures/Example/dispersal/%d.png", YYYY))
  
}
