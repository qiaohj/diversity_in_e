library(raster)
library(ggplot2)
library(dplyr)
library(dplyr)
library(Rmisc)
library(rayshader)
library(Hmisc)
library(plot3D)
library(rgl)
library(magick)
library(data.table)
library(ggpubr)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/functions.r")
source("commonFuns/colors.r")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Amphibians"
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

i=1
j=1
k=1
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
#dispersals<-data.frame(M=c(0:5), N=1)
dispersals<-c(0:1)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
SSP_i<-SSPs[1]
year_i<-2025
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-1
}
threshold_N<-0.1
threshold_left<-0.3

if (F){
  df_threshold_1<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/threshold_all.%d.rda", 1, threshold_N*100))
  df_threshold_5<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/threshold_all.%d.rda", 5, threshold_N*100))
  df_threshold_1$threshold<-1
  df_threshold_5$threshold<-5
  df_threshold<-rbind(df_threshold_1, df_threshold_5)
  
  df_threshold$type_str<-""
  
  df_threshold[which(df_threshold$type=="HIGH HIGH HIGH"), "type_str"]<-"stable hotspot"
  df_threshold[which(df_threshold$type=="HIGH HIGH LOW"), "type_str"]<-"turnover hotspot"
  
  df_threshold[which((df_threshold$low_high_2020=="HIGH")&
                       (df_threshold$low_high_year=="LOW")), "type_str"]<-"vulnerable hotspot"
  df_threshold[which((df_threshold$low_high_2020=="LOW")&
                       (df_threshold$low_high_year=="HIGH")), "type_str"]<-"novel hotspot"
  
  table(df_threshold$type_str)
  
  df_threshold_se<-df_threshold%>%dplyr::group_by(YEAR, mask_index, x, y, 
                                                  type_str, M, SSP,
                                                  threshold)%>%
    dplyr::summarise(sum_n_overlap=sum(mean_n_overlap),
                     sum_n_loss=sum(mean_n_loss),
                     sum_n_gain=sum(mean_n_gain),
                     sum_n_2020=sum(mean_n_2020),
                     sum_n_year=sum(mean_n_year))
  df_threshold_se$sum_left_2020<- df_threshold_se$sum_n_overlap/df_threshold_se$sum_n_2020
  
  saveRDS(df_threshold_se, "../../Figures/Species_gain_loss_all/Curves/data.rda")
}


df_threshold_se<-readRDS("../../Figures/Species_gain_loss_all/Curves/data.rda")
df_threshold_se<-df_threshold_se%>%dplyr::filter(M!=2)



df_sm<-df_threshold_se%>%dplyr::group_by(YEAR, type_str, M, SSP, threshold)%>%
  dplyr::summarise(N=n())
df_sm$dispersal<-ifelse(df_sm$M==0, "no dispersal", "with dispersal")
df_sm$exposure<-ifelse(df_sm$threshold==1, " no exposure", "5-year exposure")
df_sm$label<-paste(df_sm$exposure, "&", df_sm$dispersal)


colors_types<-c("HIGH HIGH HIGH"=colors_red[9],
                "HIGH LOW HIGH"=colors_green[9],
                "LOW HIGH HIGH"=colors_red[5],
                "LOW LOW HIGH"=colors_green[5],
                "HIGH HIGH LOW"=colors_purple[5],
                "HIGH LOW LOW"=colors_blue[5],
                "LOW HIGH LOW"=colors_purple[9],
                "LOW LOW LOW"=colors_blue[9])

colors_types_str<-c("stable hotspot"=colors_blue[9],
                "vulnerable hotspot"=colors_red[9],
                "novel hotspot"=colors_blue[5],
                "turnover hotspot"=colors_red[5])
p<-ggplot(df_sm%>%dplyr::filter((type_str!="")&
                                  (type_str!="turnover hotspot")&
                                  (type_str!="novel hotspot")))+
  geom_line(aes(x=YEAR, y=N, color=type_str, linetype=SSP))+
  scale_linetype_manual(values=linetype_ssp)+
  scale_color_manual(values=colors_types_str)+
  theme_bw()+
  labs(color="Type", x="Year", y="Number of pixels")+
  scale_y_log10()+
  facet_grid(exposure~dispersal, scale="free")

ggsave(p, filename="../../Figures/Species_gain_loss_all/Curves/gain_loss_all.png",
       width=16, height=8)   
  
ggsave(p, filename="../../Figures/Species_gain_loss_all/Curves/gain_loss_all.pdf",
       width=16, height=8)   

df_threshold_se<-data.table(df_threshold_se)
year_i<-2021
for (year_i in c(2100:2021)){
  print(year_i)
  df_threshold_item<-df_threshold_se[(YEAR==year_i)]
  df_threshold_item<-df_threshold_item[(type_str!="")&(type_str!="turnover hotspot")&
                                         (type_str!="novel hotspot")]
  df_threshold_item$dispersal<-ifelse(df_threshold_item$M==0, "no dispersal", "with dispersal")
  df_threshold_item$exposure<-ifelse(df_threshold_item$threshold==1, " no exposure", "5-year exposure")
  
  df_threshold_item$label<-paste(df_threshold_item$exposure, "&", df_threshold_item$dispersal)
  
  p<-ggplot()+geom_tile(data=points, aes(x=x, y=y), fill=colors_black[3])+
    geom_tile(data=df_threshold_item, aes(x=x, y=y, fill=type_str))+
    scale_fill_manual(values=colors_types_str)+
    theme_bw()+
    ggtitle(year_i)+
    labs(fill="Type")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    facet_grid(SSP~label)
  ggsave(p, filename=sprintf("../../Figures/Species_gain_loss_all/Maps/gain_loss_%d.png", year_i),
         width=16, height=8)   
  
  ggsave(p, filename=sprintf("../../Figures/Species_gain_loss_all/Maps/gain_loss_%d.pdf", year_i),
         width=16, height=8)    
}
#ffmpeg -r 2 -start_number 2021 -i gain_loss_%04d.png -y ../loss_gain.mp4
