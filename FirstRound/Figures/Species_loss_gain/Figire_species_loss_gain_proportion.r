library(raster)
library(ggplot2)
library(dplyr)
library(dplyr)
library(Rmisc)
library(rayshader)
library(Hmisc)
library(plot3D)
library( rgl )
library(magick)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("functions.r")
source("colors.r")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Reptiles"
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

i=1
j=1
k=1
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
#dispersals<-data.frame(M=c(0:5), N=1)
dispersals<-data.frame(M=c(1:2), N=1)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
SSP_i<-SSPs[1]
year_i<-2025
GCM_i<-GCMs[1]
for (SSP_i in SSPs){
  for (k in c(1:nrow(dispersals))){
    for (year_i in c(2015:2100)){
      df_end_full<-NULL
      for (GCM_i in GCMs){
        print(paste(SSP_i, k, year_i, GCM_i))
        target_folder<-sprintf("../../Objects/Diversity/%s/%s_%s_%d_%d", group, GCM_i, SSP_i, dispersals[k, "M"], dispersals[k, "N"])
        target<-sprintf("%s/loss_gain.rda", target_folder)
        df<-readRDS(target)
        df$n_2014<-df$n_loss+df$n_overlap
        df$n_year<-df$n_gain+df$n_overlap
        df$YEAR<-as.numeric(as.character(df$YEAR))
        df<-df%>%dplyr::group_by(x, y, mask_index)%>%dplyr::mutate(max_n_year=max(n_year))
        
        df_end<-df%>%dplyr::filter(YEAR==year_i)
        df_end[is.na(df_end)]<-0
        
        df_end$left_2014<-df_end$n_overlap/df_end$n_2014
        df_end[is.nan(df_end$left_2014), "left_2014"]<-0
        df_end$proportion_2014<-df_end$n_overlap/df_end$n_2014
        df_end[which(is.nan(df_end$proportion_2014)), "proportion_2014"]<-0
        df_end$proportion_year<-df_end$n_year/df_end$max_n_year
        df_end_full<-bind(df_end_full, df_end)
      }
      cuts<-seq(0, 1, by=0.05)
      df_end_se<-df_end_full%>%dplyr::group_by(YEAR, mask_index, x, y)%>%
        dplyr::summarise(mean_n_overlap=mean(n_overlap),
                         mean_n_loss=mean(n_loss),
                         mean_n_gain=mean(n_gain),
                         mean_n_2014=mean(n_2014),
                         mean_n_year=mean(n_year),
                         mean_max_n_year=mean(max_n_year),
                         mean_left_2014=mean(left_2014),
                         mean_proportion_2014=mean(proportion_2014),
                         mean_proportion_year=mean(proportion_year)
                                )
      range(df_end_se$mean_left_2014)
      
      if (F){
        table(df_end_full_se$proportion_2014_cut)
        df_end_full_se%>%dplyr::filter(between(proportion_2014_cut, 0.14, 0.16))
      }
      g1<-ggplot()+
        geom_point(data=df_end_se, aes(x=mean_proportion_2014, y=mean_proportion_year, color = mean_left_2014))+
        xlab("2014")+
        ylab(year_i)+
        scale_color_gradient(low=colors_black[4], high=colors_red[8])+
        #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersals[k, "M"], sep=""))+
        theme(legend.position = "none")
      t_dir<-sprintf("../../Figures/Species_gain_loss/Movies/2D/fine/%s/%s/%d", group, SSP_i, dispersals[k, "M"])
      dir.create(t_dir, recursive = T, showWarnings = F)
      ggsave(g1, filename=sprintf("%s/%d.png", t_dir, year_i))
      
      plot_gg(g1, multicore = F, raytrace = F, width = 5, height = 5,
              scale = 100, windowsize = c(1400,866), zoom = 0.6, phi = 30)
      t_dir<-sprintf("../../Figures/Species_gain_loss/Movies/3D/fine/%s/%s/%d", group, SSP_i, dispersals[k, "M"])
      dir.create(t_dir, recursive = T, showWarnings = F)
      render_snapshot(filename=sprintf("%s/%d.png", t_dir, year_i), clear = TRUE)
      
      df_end_se[which(df_end_se$mean_proportion_2014>1), "mean_proportion_2014"]<-1
      df_end_se[which(df_end_se$mean_proportion_year>1), "mean_proportion_year"]<-1
      
      cuts_factor<-unique(cut2(seq(0, 1, by=0.01), cuts=cuts, levels.mean=F))
      df_end_se[which(df_end_se$mean_proportion_2014>=1), "mean_proportion_2014"]<-0.99999
      df_end_se[which(df_end_se$mean_proportion_year>=1), "mean_proportion_year"]<-0.99999
      
      df_end_se$proportion_2014_cut<-cut2(df_end_se$mean_proportion_2014, cuts=cuts, levels.mean=F)
      
      df_end_se$proportion_year_cut<-cut2(df_end_se$mean_proportion_year, cuts=cuts, levels.mean=F)
      df_end_full_se<-df_end_se%>%dplyr::group_by(proportion_2014_cut, proportion_year_cut)%>%dplyr::summarise(mean_left=mean(mean_left_2014))
      
      expand_df<-expand.grid(proportion_2014_cut=cuts_factor, proportion_year_cut=cuts_factor)
      df_end_full_se<-full_join(expand_df, df_end_full_se, by=c("proportion_2014_cut", "proportion_year_cut"))
      df_end_full_se[is.na(df_end_full_se$mean_left), "mean_left"]<-0
      
      gg<-ggplot()+
        geom_tile(data=df_end_full_se, aes(x=proportion_2014_cut, y=proportion_year_cut, fill = mean_left), color=colors_black[7])+
        xlab("2014")+
        ylab(year_i)+
        scale_fill_gradient(low=colors_black[4], high=colors_red[8])+
        #ggtitle(paste(group, ", ", SSP_i, ", Dispersal distance=", dispersals[k, "M"], sep=""))+
        theme(legend.position = "none",
              axis.text.x = element_text( angle = 90 ))
      t_dir<-sprintf("../../Figures/Species_gain_loss/Movies/2D/rough/%s/%s/%d", group, SSP_i, dispersals[k, "M"])
      dir.create(t_dir, recursive = T, showWarnings = F)
      ggsave(gg, filename=sprintf("%s/%d.png", t_dir, year_i))
      
      plot_gg(gg, multicore = F, raytrace = F, width = 5, height = 5,
              scale = 100, windowsize = c(1400,866), zoom = 0.6, phi = 30)
      t_dir<-sprintf("../../Figures/Species_gain_loss/Movies/3D/rough/%s/%s/%d", group, SSP_i, dispersals[k, "M"])
      dir.create(t_dir, recursive = T, showWarnings = F)
      render_snapshot(filename=sprintf("%s/%d.png", t_dir, year_i), clear = TRUE)
    }
  }
}