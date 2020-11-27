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
dispersals<-c(0)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
SSP_i<-SSPs[1]
year_i<-2025
GCM_i<-GCMs[1]
output_figures<-T
threshold<-5
for (SSP_i in SSPs){
  for (k in c(1:length(dispersals))){
    if (T){
      df_end_full_list<-list()
      max_2020<-0
      max_year<-0
      for (year_i in c(2021:2100)){
        df_end_full<-NULL
        for (GCM_i in GCMs){
          print(paste(group, SSP_i, dispersals[k], year_i, GCM_i, threshold))
          target_folder<-sprintf("../../Objects/Diversity_%d/%s/%s_%s_%d", threshold, group, GCM_i, SSP_i, dispersals[k])
          target<-sprintf("%s/loss_gain.rda", target_folder)
          df<-readRDS(target)
          
          df$n_2020<-df$n_loss+df$n_overlap
          df$n_year<-df$n_gain+df$n_overlap
          df$YEAR<-as.numeric(as.character(df$YEAR))
          #print(range(df$n_year))
          
          df_end<-df%>%dplyr::filter(YEAR==year_i)
          df_end[is.na(df_end)]<-0
          
          df_end$left_2020<-df_end$n_overlap/df_end$n_2020
          df_end[is.nan(df_end$left_2020), "left_2020"]<-0
          
          df_end_full<-bind_dplyr(df_end_full, df_end)
        }
        df_end_se<-df_end_full%>%dplyr::group_by(YEAR, mask_index, x, y)%>%
          dplyr::summarise(mean_n_overlap=mean(n_overlap),
                           mean_n_loss=mean(n_loss),
                           mean_n_gain=mean(n_gain),
                           mean_n_2020=mean(n_2020),
                           mean_n_year=mean(n_year),
                           mean_left_2020=mean(left_2020)
          )
        range(df_end_se$mean_n_year)
        
        max_temp<-ceiling(max(df_end_se$mean_n_2020))
        max_2020<-ifelse(max_2020>max_temp, max_2020, max_temp) 
        
        max_temp<-ceiling(max(df_end_se$mean_n_year))
        max_year<-ifelse(max_year>max_temp, max_year, max_temp) 
        
        df_end_full_list[[as.character(year_i)]]<-df_end_se
        
      }
      saveRDS(data.frame(max_2020=max_2020, max_year=max_year), 
              sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d_min_max.rda", threshold, group, SSP_i, dispersals[k]))
      saveRDS(df_end_full_list, 
              sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d.rda", threshold, group, SSP_i, dispersals[k]))
    }
    if (!output_figures){
      next()
    }
    df_end_full_list<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d.rda", threshold, group, SSP_i, dispersals[k]))
    min_max<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d_min_max.rda", threshold, group, SSP_i, dispersals[k]))
    year_str<-"2021"
    for (year_str in names(df_end_full_list)){
      df_end_se<-df_end_full_list[[year_str]]
      g1<-ggplot()+
        geom_point(data=df_end_se, aes(x=mean_n_2020, y=mean_n_year, color = mean_left_2020))+
        xlab("2020")+
        ylab(year_str)+
        xlim(c(0, min_max$max_2020))+
        ylim(c(0, min_max$max_year))+
        scale_color_gradient(low=colors_black[4], high=colors_red[8])+
        #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersals[k], sep=""))+
        theme(legend.position = "none")
      t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/2D/fine/%s/%s/%d", threshold, group, SSP_i, dispersals[k])
      dir.create(t_dir, recursive = T, showWarnings = F)
      ggsave(g1, filename=sprintf("%s/%s.png", t_dir, year_str))
      
      plot_gg(g1, multicore = F, raytrace = F, width = 5, height = 5,
              scale = 100, windowsize = c(1400,866), zoom = 0.6, phi = 30)
      t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/3D/fine/%s/%s/%d", threshold, group, SSP_i, dispersals[k])
      dir.create(t_dir, recursive = T, showWarnings = F)
      render_snapshot(filename=sprintf("%s/%s.png", t_dir, year_str), clear = TRUE)
      
      cuts_x<-seq(0, ceiling(min_max$max_2020/10)*10, by=10)
      cuts_y<-seq(0, ceiling(min_max$max_year/10)*10, by=10)
      
      df_end_se$p_2020_cut<-cut2(df_end_se$mean_n_2020+0.001, cuts=cuts_x, levels.mean=F)
      df_end_se$p_year_cut<-cut2(df_end_se$mean_n_year+0.001, cuts=cuts_y, levels.mean=F)
      df_end_full_se<-df_end_se%>%dplyr::group_by(p_2020_cut, p_year_cut)%>%dplyr::summarise(mean_left=mean(mean_left_2020))
      cuts_factor_x<-cut2(c(1:max(cuts_x)), cuts=cuts_x)
      cuts_factor_y<-cut2(c(1:max(cuts_y)), cuts=cuts_y)
      
      expand_df<-expand.grid(p_2020_cut=cuts_factor_x, p_year_cut=cuts_factor_y)
      df_end_full_se<-full_join(expand_df, df_end_full_se, by=c("p_2020_cut", "p_year_cut"))
      df_end_full_se[is.na(df_end_full_se$mean_left), "mean_left"]<-0
      
      gg<-ggplot()+
        geom_tile(data=df_end_full_se, aes(x=p_2020_cut, y=p_year_cut, fill = mean_left), color=colors_black[7])+
        xlab("2020")+
        ylab(year_str)+
        scale_fill_gradient(low=colors_black[4], high=colors_red[8])+
        #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersals[k], sep=""))+
        theme(legend.position = "none",
              axis.text.x = element_text( angle = 90 ))
      t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/2D/rough/%s/%s/%d", threshold, group, SSP_i, dispersals[k])
      dir.create(t_dir, recursive = T, showWarnings = F)
      ggsave(gg, filename=sprintf("%s/%s.png", t_dir, year_str))
      
      plot_gg(gg, multicore = F, raytrace = F, width = 5, height = 5,
              scale = 100, windowsize = c(1400,866), zoom = 0.6, phi = 30)
      t_dir<-sprintf("../../Figures/Species_gain_loss_%d/Movies/RawValue/3D/rough/%s/%s/%d", threshold, group, SSP_i, dispersals[k])
      dir.create(t_dir, recursive = T, showWarnings = F)
      render_snapshot(filename=sprintf("%s/%s.png", t_dir, year_str), clear = TRUE)
    }
    
  }
}

#/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Species_gain_loss_5/Movies/RawValue/3D/rough/ALL
#ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../../../gain_loss_3D.mp4