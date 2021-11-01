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
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/functions.r")
source("commonFuns/colors.r")
args = commandArgs(trailingOnly=TRUE)
threshold<-as.numeric(args[1])
if (is.na(threshold)){
  threshold<-1
}

dispersal<-as.numeric(args[2])
if (is.na(dispersal)){
  dispersal<-1
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
mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
SSP_i<-SSPs[3]
year_i<-2100
GCM_i<-GCMs[1]
output_figures<-T

for (SSP_i in SSPs){
  df_end_full_all<-NULL
  for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    print(paste(SSP_i, group))
    df_end_full_list<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d.rda", 
                                      threshold, group, SSP_i, dispersal))
    df_end_full<-rbindlist(df_end_full_list)
    df_end_full_all<-bind_dplyr(df_end_full_all, df_end_full)
  }
  df_end_full_all_final<-df_end_full_all%>%dplyr::group_by(YEAR, mask_index, x, y)%>%
    dplyr::summarise(sum_n_overlap=sum(mean_n_overlap),
                     sum_n_loss=sum(mean_n_loss),
                     sum_n_gain=sum(mean_n_gain),
                     sum_n_2020=sum(mean_n_2020),
                     sum_n_year=sum(mean_n_year),
                     sum_left_2020=sum(mean_left_2020))
  
  year_str<-2100
  min_max<-data.frame(max_2020=max(df_end_full_all_final$sum_n_2020),
                      max_year=max(df_end_full_all_final$sum_n_year))
  for (year_str in c(2021:2100)){
    df_end_se<-df_end_full_all_final%>%dplyr::filter(YEAR==year_str)
    g1<-ggplot()+
      geom_point(data=df_end_se, aes(x=sum_n_2020, y=sum_n_year, color = sum_left_2020))+
      xlab("2020")+
      ylab(year_str)+
      xlim(c(0, min_max$max_2020))+
      ylim(c(0, min_max$max_year))+
      scale_color_gradient(low=colors_black[4], high=colors_red[8])+
      #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersal, sep=""))+
      theme(legend.position = "none")
    t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/2D/fine/ALL/%s/%d", threshold, SSP_i, dispersal)
    dir.create(t_dir, recursive = T, showWarnings = F)
    ggsave(g1, filename=sprintf("%s/%s.png", t_dir, year_str))
    
    plot_gg(g1, multicore = F, raytrace = F, width = 5, height = 5,
            scale = 100, windowsize = c(1400,866), zoom = 0.5, phi = 30)
    t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/3D/fine/ALL/%s/%d", threshold, SSP_i, dispersal)
    dir.create(t_dir, recursive = T, showWarnings = F)
    render_snapshot(filename=sprintf("%s/%s.png", t_dir, year_str), clear = TRUE)
    
    cuts_x<-seq(0, ceiling(min_max$max_2020/10)*10, by=ceiling(min_max$max_2020/10)/3)
    cuts_y<-seq(0, ceiling(min_max$max_year/10)*10, by=ceiling(min_max$max_year/10)/3)
    
    df_end_se$p_2020_cut<-cut(df_end_se$sum_n_2020, breaks=cuts_x, labels=cuts_x[2:length(cuts_x)], include.lowest=T)
    df_end_se$p_year_cut<-cut(df_end_se$sum_n_year, breaks=cuts_y, labels=cuts_y[2:length(cuts_y)], include.lowest=T)
    df_end_full_se<-df_end_se%>%dplyr::group_by(p_2020_cut, p_year_cut)%>%dplyr::summarise(sum_left=mean(sum_left_2020, na.rm=T))
    cuts_factor_x<-unique(cut(c(1:max(cuts_x)), breaks=cuts_x, labels=cuts_x[2:length(cuts_x)], include.lowest=T))
    cuts_factor_x<-cuts_factor_x[!is.na(cuts_factor_x)]
    cuts_factor_y<-cut(c(1:max(cuts_y)), breaks=cuts_y, labels=cuts_y[2:length(cuts_y)], include.lowest=T)
    cuts_factor_y<-cuts_factor_y[!is.na(cuts_factor_y)]
    
    expand_df<-expand.grid(p_2020_cut=cuts_factor_x, p_year_cut=cuts_factor_y)
    df_end_full_se<-full_join(expand_df, df_end_full_se, by=c("p_2020_cut", "p_year_cut"))
    df_end_full_se[is.na(df_end_full_se$sum_left), "sum_left"]<-0
    df_end_full_se$p_2020_cut<-as.numeric(as.character(df_end_full_se$p_2020_cut))
    df_end_full_se$p_year_cut<-as.numeric(as.character(df_end_full_se$p_year_cut))
    df_end_full_se[which(is.na(df_end_full_se$p_year_cut)),]
    gg<-ggplot()+
      geom_tile(data=df_end_full_se, aes(x=p_2020_cut, y=p_year_cut, fill = sum_left), 
                color=colors_black[7])+
      xlab("2020")+
      ylab(year_str)+
      scale_fill_gradient(low=colors_black[4], high=colors_red[8])+
      #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersal, sep=""))+
      theme(legend.position = "none",
            axis.text.x = element_text( angle = 90 ))
    t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/2D/rough/ALL/%s/%d", threshold, SSP_i, dispersal)
    dir.create(t_dir, recursive = T, showWarnings = F)
    ggsave(gg, filename=sprintf("%s/%s.png", t_dir, year_str))
    
    plot_gg(gg, multicore = F, raytrace = F, width = 5, height = 5,
            scale = 100, windowsize = c(1400,866), zoom = 0.5, phi = 30)
    t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/3D/rough/ALL/%s/%d", threshold, SSP_i, dispersal)
    dir.create(t_dir, recursive = T, showWarnings = F)
    render_snapshot(filename=sprintf("%s/%s.png", t_dir, year_str), clear = TRUE)
  }
  

}

#/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Species_gain_loss_5/Movies/RawValue/3D/rough/ALL
#ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../../../gain_loss_3D.mp4