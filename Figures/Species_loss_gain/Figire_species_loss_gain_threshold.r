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
dispersals<-c(0:2)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")
SSP_i<-SSPs[1]
year_i<-2025
threshold<-1
threshold_N<-0.1
threshold_left<-0.3
for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df_end_full_all<-NULL
  for (SSP_i in SSPs){
    for (k in c(1:length(dispersals))){
      print(paste(group, SSP_i, dispersals[k]))
      df_end_full_list<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/%s_%s_%d.rda", threshold, group, SSP_i, dispersals[k]))
      df_end_full<-rbindlist(df_end_full_list)
      threshold_N_2020<-quantile(df_end_full[(mean_n_2020>0)&(YEAR==2100),]$mean_n_2020, c(threshold_N, 1-threshold_N))
      threshold_N_year<-quantile(df_end_full[(mean_n_year>0),]$mean_n_year, c(threshold_N, 1-threshold_N))
      df_end_full$low_high_2020<-"MIDDLE"
      df_end_full[mean_n_2020<=threshold_N_2020[1], "low_high_2020"]<-"LOW"
      df_end_full[mean_n_2020>=threshold_N_2020[2], "low_high_2020"]<-"HIGH"
      df_end_full$low_high_year<-"MIDDLE"
      df_end_full[mean_n_year<=threshold_N_year[1], "low_high_year"]<-"LOW"
      df_end_full[mean_n_year>=threshold_N_year[2], "low_high_year"]<-"HIGH"
      table(df_end_full$low_high_year)
      
      df_end_full$low_high_left<-"MIDDLE"
      df_end_full[mean_left_2020<=threshold_left, "low_high_left"]<-"LOW"
      df_end_full[mean_left_2020>=(1-threshold_left), "low_high_left"]<-"HIGH"
      table(df_end_full$low_high_left)
      
      df_end_full$type<-paste(df_end_full$low_high_2020, df_end_full$low_high_year, df_end_full$low_high_left)
      table(df_end_full$type)
      
      df_end_full$M<-dispersals[k]
      df_end_full$SSP<-SSP_i
      df_end_full_all<-bind(df_end_full_all, df_end_full)
    }
  }
  saveRDS(df_end_full_all, file = sprintf("../../Figures/Species_gain_loss_%d/%s_threshold.%d.rda", threshold, group, threshold_N*100))
}
types<-c("LOW", "HIGH")
types_list<-expand.grid(a=types, b=types, c=types)
types_list$type<-paste(types_list$a, types_list$b, types_list$c)
df_threshold<-NULL

for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(group)
  df_end_full_all<-readRDS(file = sprintf("../../Figures/Species_gain_loss_%d/%s_threshold.%d.rda", threshold, group, threshold_N*100))
  df_end_full_all_sub<-df_end_full_all[type %in% types_list$type]
  df_end_full_all_sub$group<-group
  df_threshold<-bind(df_threshold, df_end_full_all_sub)
  
}
df_threshold$type_index<-as.factor(df_threshold$type)
df_threshold$type_index_number<-as.numeric(df_threshold$type_index)


saveRDS(df_threshold, file = sprintf("../../Figures/Species_gain_loss_%d/threshold_all.%d.rda", threshold, threshold_N*100))


df_threshold<-readRDS(sprintf("../../Figures/Species_gain_loss_%d/threshold_all.%d.rda", threshold, threshold_N*100))
colors_types<-c("HIGH HIGH HIGH"=colors_red[9],
                "HIGH LOW HIGH"=colors_red[7],
                "LOW HIGH HIGH"=colors_red[5],
                "LOW LOW HIGH"=colors_red[3],
                "HIGH HIGH LOW"=colors_blue[3],
                "HIGH LOW LOW"=colors_blue[5],
                "LOW HIGH LOW"=colors_blue[7],
                "LOW LOW LOW"=colors_blue[9])



#unique(df_sm[, c("type", "type_index_number")])

df_sm<-df_threshold%>%dplyr::group_by(YEAR, type, M, SSP, type_index, type_index_number, group)%>%
  dplyr::summarise(N=n())
p<-ggplot(df_sm%>%dplyr::filter(M==1))+geom_line(aes(x=YEAR, y=N, color=type, linetype=SSP))+
  scale_linetype_manual(values=linetype_ssp)+
  scale_color_manual(values=colors_types)+
  theme_bw()+
  facet_wrap(~group, scale="free", ncol=1)

ggsave(p, filename=sprintf("../../Figures/Species_gain_loss_%d/threshold.%d.png", threshold, threshold_N*100),
                           width=8, height=8)  
ggsave(p, filename=sprintf("../../Figures/Species_gain_loss_%d/threshold.%d.pdf", threshold, threshold_N*100),
       width=8, height=8)  
group_i<-"Amphibians"
p<-ggplot()+geom_tile(data=points, aes(x=x, y=y), fill=colors_black[3])+
  geom_tile(data=df_threshold, aes(x=x, y=y, fill=type))+
  scale_fill_manual(values=colors_types)
legend_g<-get_legend(p)


for (year_i in c(2021:2100)){
  print(year_i)
  p_list<-list()
  for (group_i in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    for (SSP_i in SSPs){
      df_threshold_item<-df_threshold[(YEAR==year_i)&(M==1)&(SSP==SSP_i)&(group==group_i)]
      p<-ggplot()+geom_tile(data=points, aes(x=x, y=y), fill=colors_black[3])+
        geom_tile(data=df_threshold_item, aes(x=x, y=y, fill=type))+
        scale_fill_manual(values=colors_types)+
        ggtitle(paste(year_i, group_i, SSP_i))+
        theme_bw()+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position = "none")
      p_list[[paste(group_i, SSP_i)]]<-p
    }
    
  }
  ppp<-ggarrange(plotlist=p_list, nrow=4, ncol=3, common.legend = T, legend = "right", legend.grob = legend_g)
  #ggsave(ppp, 
  #       filename=sprintf("../../Figures/Species_gain_loss_%d/threshold_map_year.%d/%d.png", threshold, threshold_N*100, year_i),
  #       width=12, height = 9)
  ggsave(ppp, 
         filename=sprintf("../../Figures/Species_gain_loss_%d/threshold_map_year.%d/%d.pdf", threshold, threshold_N*100, year_i),
         width=12, height = 9)
}

  #ffmpeg -r 2 -start_number 2021 -i %04d.png -y ../threshold_map_year.10.mp4
