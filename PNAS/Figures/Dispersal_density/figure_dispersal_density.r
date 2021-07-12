library(ceramic)
library(raster)
library(data.table)
library(igraph)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(tidyr)
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
rm(list=ls())
source("commonFuns/colors.r")
source("commonFuns/functions.r")

setDTthreads(threads=1)
print(sprintf("%d CPUs are using", getDTthreads()))


GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
l_i=1
exposure=0
GCM<-GCMs[1]
SSP<-SSPs[1]
mask_p<-data.table(rasterToPoints(raster("../../Raster/mask_100km_plot.tif")))
colnames(mask_p)[3]<-"mask_100km"
df_final<-NULL
for (SSP_i in c(1:length(SSPs))){
  SSP<-SSPs[SSP_i]
  for (exposure in c(0, 5)){
    df_all<-list()
    for (GCM_i in c(1:length(GCMs))){
      GCM<-GCMs[GCM_i]
      target_folder<-sprintf("../../Objects/density_based_pathway/%s_%s_exposure_%d.rda", 
                             GCM, SSP, exposure)
      df<-readRDS(target_folder)
      df<-rbindlist(df)
      if (F){
        ggplot(df[["2100"]])+geom_tile(aes(x=x, y=y, fill=count))
        ggplot(df_all_all_year)+geom_tile(aes(x=x, y=y, fill=mean_count))
      }
      df_all[[GCM]]<-df
    }
    df_all<-rbindlist(df_all)
    df_all_all_year<-df_all[, .(mean_count=mean(count)), by=list(x, y, mask_100km, SSP, exposure)]
    df_all_all_year<-df_all_all_year[mask_100km %in% mask_p$mask_100km]
    df_all_all_year$SSP<-SSP
    df_final<-bind(df_final, df_all_all_year)
  }
}
df_final$exposure_label<-ifelse(df_final$exposure==0, " no exposure", "5-year exposure")
df_final$SSP<-as.character(df_final$SSP)
p<-ggplot()+
  geom_tile(data=mask_p, aes(x=x, y=y), fill="#e5f5f9")+
  #geom_tile(data=df_final%>%dplyr::filter((exposure==5)&(SSP=="SSP119")), aes(x=x, y=y, fill=mean_count))+
  geom_tile(data=df_final, aes(x=x, y=y, fill=mean_count))+
  facet_grid(exposure_label~SSP, scale="free")+
  scale_fill_gradient(low="#e5f5f9", high="#2ca25f")+
  
  
  labs(fill = "Number of uses for dispersal")+
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.background = element_rect(fill = map_background, color = NA), 
    panel.background = element_blank(), 
    #legend.background = element_rect(fill = map_background, color = NA),
    #panel.border = element_blank(),
    #panel.border = element_rect(colour = "black", fill=NA),
    #strip.background = element_blank()
  )
p
ggsave(p, filename="../../Figures/density_based_pathway/density_based_pathway.png", width=13, height=6)

alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
xy<-c("x", "y")
df_final$alt<-raster::extract(alt, df_final[,..xy])
df_final$slope<-raster::extract(slope, df_final[,..xy])
df_final<-df_final[!is.na(alt)]
df_final<-df_final[!is.na(slope)]

alt_p<-data.frame(rasterToPoints(alt))
df_final_alt<-full_join(df_final, alt_p, by=c("x", "y"))
df_final_alt[is.na(mean_count)]$mean_count<-0
df_final_alt<-df_final_alt[!is.na(alt)]

cor(df_final_alt$mean_count, df_final_alt$alt)

p<-ggplot(df_final_alt, aes(x=alt, y=mean_count))+geom_point(size=1, color="grey")+
  geom_smooth(method="lm")+
  facet_grid(exposure_label~SSP, scale="free")+
  xlab("Elevation (m)")+
  ylab("Number of uses for dispersal")+
  theme_bw()
p
ggsave(p, filename="../../Figures/density_based_pathway/density_vs_elevation.png", width=13, height=6)
p<-ggplot(df_final, aes(x=slope, y=mean_count))+geom_point(size=1, color="grey")+
  geom_smooth(method="lm")+
  facet_grid(exposure_label~SSP, scale="free")+
  xlab("Slope")+
  ylab("Number of uses for dispersal")+
  theme_bw()
p
ggsave(p, filename="../../Figures/density_based_pathway/density_vs_slope.png", width=13, height=6)

for (SSP_i in c(1:length(SSPs))){
  SSP_l<-SSPs[SSP_i]
  for (exposure_l in c(0, 5)){
    df_final_item<-df_final[((SSP==SSP_l)&(exposure==exposure_l))]
    p_item<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill="#e5f5f9")+
      geom_tile(data=df_final_item, aes(x=x, y=y, fill=mean_count))+
      scale_fill_gradient(low="#e5f5f9", high="#2ca25f")+
      labs(fill = "Number of uses for dispersal")+
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.background = element_rect(fill = map_background, color = NA), 
        panel.background = element_blank(), 
        #legend.background = element_rect(fill = map_background, color = NA),
        #panel.border = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA),
        #strip.background = element_blank()
      )
    p_item
    ggsave(p_item, filename=
             sprintf("../../Figures/density_based_pathway/density_based_pathway_%s_exposure_%d.png", SSP_l, exposure_l),
           width=12, height=6)
  }
}
if (F){
  SSP_i=1
  df_all<-list()
  for (SSP_i in c(1:length(SSPs))){
    SSP<-SSPs[SSP_i]
    for (exposure in c(0, 5)){
      for (GCM_i in c(1:length(GCMs))){
        GCM<-GCMs[GCM_i]
        target_folder<-sprintf("../../Objects/density_based_pathway/%s_%s_exposure_%d.rda", 
                               GCM, SSP, exposure)
        df<-readRDS(target_folder)
        df<-rbindlist(df)
        if (F){
          ggplot(df[["2100"]])+geom_tile(aes(x=x, y=y, fill=count))
          ggplot(df_all_all_year)+geom_tile(aes(x=x, y=y, fill=mean_count))
        }
        df_all[[paste(GCM, SSP, exposure)]]<-df
      }
    }
  }
  df_all<-rbindlist(df_all)
  df_all_by_year<-df_all[, .(mean_count=mean(count)), by=list(x, y, mask_100km, SSP, exposure, year)]
  df_all_by_year<-df_all_by_year[mask_100km %in% mask_p$mask_100km]
  saveRDS(df_all_by_year, "../../Figures/density_based_pathway/df_all_by_year.rda")
}

df_all_by_year<-readRDS("../../Figures/density_based_pathway/df_all_by_year.rda")
df_all_by_year<-df_all_by_year[mask_100km %in% mask_p$mask_100km]
df_all_by_year$exposure<-ifelse(df_all_by_year$exposure==0, " no exposure", "5-year exposure")
df_all_by_year$SSP<-as.character(df_all_by_year$SSP)
yyy=2021
max_count<-max(df_all_by_year$mean_count)
for (yyy in c(2021:2100)){
  print(yyy)
  if (F){
    ggplot(df_all_by_year)+geom_histogram(aes(x=mean_count))+
      facet_grid(SSP~exposure)
  }
  p<-ggplot()+
    geom_tile(data=mask_p, aes(x=x, y=y), fill="#e5f5f9")+
    geom_tile(data=df_all_by_year[year==yyy], aes(x=x, y=y, fill=mean_count))+
    scale_fill_gradient(low="#e5f5f9", high="#2ca25f",
                        limits=c(0, 100), oob=squish,
                        breaks=seq(0, 100, by=20),
                        labels=c(as.character(seq(0, 100, by=20)[1:5]), 
                                 sprintf(">100, up to %d", max_count)))+
    labs(fill = "Number of uses for dispersal")+
    ggtitle(yyy)+
    facet_grid(exposure~SSP)+
    theme(
      legend.position = "bottom",
      legend.key.width=unit(0.8, "in"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #plot.background = element_rect(fill = map_background, color = NA), 
      panel.background = element_blank(), 
      #legend.background = element_rect(fill = map_background, color = NA),
      #panel.border = element_blank(),
      #panel.border = element_rect(colour = "black", fill=NA),
      #strip.background = element_blank()
    )
  ggsave(p, filename = sprintf("../../Figures/density_based_pathway/Movies/%d.png", yyy),
         width=13, height=6)
}

#cd /media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/density_based_pathway/Movies
#ffmpeg -r 5 -start_number 2021 -i %04d.png -y ../density_based_pathway.mp4

