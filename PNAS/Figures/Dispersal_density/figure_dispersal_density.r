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
  geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  #geom_tile(data=df_final%>%dplyr::filter((exposure==5)&(SSP=="SSP119")), aes(x=x, y=y, fill=mean_count))+
  geom_tile(data=df_final, aes(x=x, y=y, fill=mean_count))+
  facet_grid(exposure_label~SSP, scale="free")+
  scale_fill_gradient(low=color_two_map[1], high=color_two_map[2])+
  
  
  labs(fill = "Number of usage")+
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

for (SSP_i in c(1:length(SSPs))){
  SSP_l<-SSPs[SSP_i]
  for (exposure_l in c(0, 5)){
    df_final_item<-df_final[((SSP==SSP_l)&(exposure==exposure_l))]
    p_item<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(data=df_final_item, aes(x=x, y=y, fill=mean_count))+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2])+
      labs(fill = "Number of usage")+
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
