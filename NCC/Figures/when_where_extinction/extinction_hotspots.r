library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(scales)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
upper_threshold<-0.95

if (T){
  g<-"Amphibians"
  mask<-raster("../../Raster/mask_index.tif")
  mask_p<-data.frame(rasterToPoints(mask))
  threshold<-1
  SSP_i<-"SSP585"
  da=0
  
  result<-NULL
  upper_percentiles<-list()
  for (threshold in c(1, 5)){
    for (ttt in c(0, 1, 2)){
      for (SSP_i in c("SSP119", "SSP245", "SSP585")){
        for (da in c(0:1)){
          n_ext_all<-NULL
          for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
            
            print(paste(threshold, SSP_i, g, da, ttt))
            
            n_ext_r<-raster(sprintf("../../Figures_Full_species/when_where_extinction_%d/where/%s_%s_n_extinct_dispersal_%d_ttt_%d.tif", 
                                    threshold, g, SSP_i, da, ttt))
            
            n_ext<-data.frame(rasterToPoints(n_ext_r))
            colnames(n_ext)[3]<-"V"
            n_ext<-left_join(mask_p, n_ext, by=c("x", "y"))
            n_ext$group<-g
            n_ext$ttt<-ttt
            n_ext[which(is.na(n_ext$V)), "V"]<-0
            n_ext_all<-bind_dplyr(n_ext_all, n_ext)
          }
          n_ext_all<-n_ext_all%>%dplyr::group_by(x, y, mask_index, ttt)%>%
            dplyr::summarise(sum_V=sum(V))
          if (da==0){
            upper_percentile<-quantile(pull(n_ext_all[which(n_ext_all$sum_V>0), "sum_V"]), upper_threshold)
            upper_percentiles[[paste(threshold, ttt, SSP_i)]]<-upper_percentile
          }else{
            upper_percentile<-upper_percentiles[[paste(threshold, ttt, SSP_i)]]
          }
          hotspots<-n_ext_all%>%dplyr::filter(sum_V>upper_percentile)
          if (F){
            ggplot(hotspots)+geom_point(aes(x=x, y=y, color=sum_V))
          }
          hotspots$upper_threshold<-upper_threshold
          hotspots$upper_percentile<-upper_percentile
          hotspots$SSP<-SSP_i
          hotspots$dispersal<-ifelse(da==0, "no dispersal", "with dispersal")
          hotspots$ttt<-ttt
          hotspots$exposure<-ifelse(threshold==1, " no exposure", "5-year exposure")
          result<-bind_dplyr(result, hotspots)
        }
      }
    }
  }
  saveRDS(result, sprintf("../../Figures_Full_species/Extinction_hotspots/Extinction_hotspots_%.2f.rda", upper_threshold))
}
hotspots<-readRDS(sprintf("../../Figures_Full_species/Extinction_hotspots/Extinction_hotspots_%.2f.rda", upper_threshold))

exposure_i<-"no exposure"
ttt_i=1
SSP_i="SSP119"
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
for (ttt_i in c(0, 1, 2)){
  df<-NULL
  for (exposure_i in c(" no exposure", "5-year exposure")){
    for (SSP_i in c("SSP119", "SSP245", "SSP585")){
      item<-hotspots%>%dplyr::filter((exposure==exposure_i)&(ttt==ttt_i)&(SSP==SSP_i))
      item1<-item%>%dplyr::filter(dispersal=="no dispersal")
      item2<-item%>%dplyr::filter(dispersal=="with dispersal")
      merged_item<-left_join(item1, item2, by=c("x", "y", "mask_index", "ttt",
                                                "upper_threshold", "upper_percentile",
                                                "SSP", "exposure"))
      merged_item_map<-merged_item%>%dplyr::filter(is.na(sum_V.y))
      df<-bind_dplyr(df, merged_item_map)
    }
  }
  title1<-sprintf("Distribution>%d", ttt_i)
  p<-ggplot(df)+
    geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
    geom_tile(aes(x=x, y=y, fill=factor(round(upper_percentile))))+
    facet_grid(exposure~SSP)+
    ggtitle(title1)+
    labs(fill = sprintf("%.2f percentile", upper_threshold))+
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = map_background, color = NA), 
      panel.background = element_blank(), 
      legend.background = element_rect(fill = map_background, color = NA),
      panel.border = element_blank()
    )
  p
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_hotspots/ttt_%d_%.2f.png", ttt_i, upper_threshold),
         width=9, height=4)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_hotspots/ttt_%d_%.2f.pdf", ttt_i, upper_threshold),
         width=9, height=4)
}

