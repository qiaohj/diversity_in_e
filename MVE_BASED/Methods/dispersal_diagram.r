library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(ggpubr)

setwd("Y:/Script/diversity_in_e")

source("colors.R")
source("genCircle.R")
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}

group<-"Amphibians"

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(2000:2014)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_2000_2014_df.rda")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))


var_pair<-start_layer_df[which(start_layer_df$VAR=="pr"),]
var_pair1<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmax"),], by=c("GCM", "SSP", "Y"))
#var_pair2<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmin"),], by=c("GCM", "SSP", "Y"))
var_pair1<-var_pair1%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair1)<-c("GCM", "PR", "Y", "TEMP")
#var_pair2<-var_pair2%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
#colnames(var_pair2)<-c("GCM", "PR", "Y", "TEMP")
#var_pair<-bind_rows(var_pair1, var_pair2)
var_pair<-var_pair1
var_pair$PR_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$PR, var_pair$Y, sep="_")
var_pair$TEMP_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$TEMP, var_pair$Y, sep="_")
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
i=1

if (F){
  j=1
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area<=0){
      next()
    }
    
    
    print(paste(i, nrow(df_list), item$sp))
    occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
    if ((nrow(occ)>20)&(nrow(occ)<25)){
      if (j==1){
        asdf
      }
      j=j+1
    }else{
      next()
    }
  }
}
example_sp<-"Dendropsophus_walfordi"
target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
all_v<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
#plot(mask)
target<-sprintf("%s/dispersal", target_folder)
dispersal_df<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1_1"))

year=2015
for (year in c(2015:2099)){
  print(year)
  dis1<-dispersal_df%>%dplyr::filter(YEAR==year)
  p1<-left_join(mask_p, dis1, by=c("x", "y", "mask_index"))
  
  
  next_dis1<-mask_p%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, dis1)/100000)
  next_dis1<-next_dis1%>%filter(dist<=1)
  p_next_1<-left_join(mask_p, next_dis1, by=c("x", "y", "mask_index"))
  
  
  dis2<-dispersal_df%>%dplyr::filter(YEAR==year+1)
  p2<-left_join(mask_p, dis2, by=c("x", "y", "mask_index"))
  
  
  
  g2<-ggplot()+
    geom_tile(data=mask_p, aes(x=x, y=y), fill=colors_black[3])+
    geom_tile(data=p1 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_black[8])+
    geom_tile(data=p_next_1 %>% dplyr::filter(!is.na(dist)), aes(x=x, y=y), fill=colors_blue[8], alpha=0.5)+
    geom_tile(data=p2 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_red[8], alpha=0.5)+
    xlim(c(-8000000, -3000000))+
    ylim(c(-3300000, 1500000))+
    ggtitle(year)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  
  if (F){
    g_item1<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=colors_black[3])+
      geom_tile(data=p1 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_black[8])+
      xlim(c(-8000000, -3000000))+
      ylim(c(-3300000, 1500000))+
      ggtitle("Before dispersal")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none")
    
    g_item2<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=colors_black[3])+
      geom_tile(data=p1 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_black[8])+
      geom_tile(data=p_next_1 %>% dplyr::filter(!is.na(dist)), aes(x=x, y=y), fill=colors_blue[8], alpha=0.5)+
      xlim(c(-8000000, -3000000))+
      ylim(c(-3300000, 1500000))+
      ggtitle("After dispersal, before detecting the suitable area")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none")
    
    g_item3<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=colors_black[3])+
      geom_tile(data=p2 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_red[8], alpha=0.5)+
      xlim(c(-8000000, -3000000))+
      ylim(c(-3300000, 1500000))+
      ggtitle("After remoing the unsuitable area")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none")
    
    g_item4<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=colors_black[3])+
      geom_tile(data=p1 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_black[8])+
      geom_tile(data=p_next_1 %>% dplyr::filter(!is.na(dist)), aes(x=x, y=y), fill=colors_blue[8], alpha=0.5)+
      geom_tile(data=p2 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_red[8], alpha=0.5)+
      xlim(c(-8000000, -3000000))+
      ylim(c(-3300000, 1500000))+
      ggtitle("Combined")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none")
    
    ppp<-ggarrange(g_item1, g_item2, g_item3, g_item4)
    ggsave(ppp, filename=sprintf("../../Figures/Methods/example_%d.png", year), width = 12, height=10)
  }
  #g2
  ggsave(g2, filename=sprintf("../../Figures/Methods/dispersal_diagram/%d.png", year), width = 6, height=5)
}

#ffmpeg -r 1 -start_number 2015 -i %04d.png -y ../dispersal_diagram.gif
#ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../dispersal_diagram.mp4
