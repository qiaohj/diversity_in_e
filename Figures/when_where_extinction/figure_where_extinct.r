library(dplyr)
library(raster)
library(gglpot2)
library(Rmisc)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

if (F){
  g<-"Birds"
  for (threshold in c(1, 5)){
    for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
      print(paste(g, threshold))
      df<-readRDS(sprintf("../../Objects/when_where_extinction_%d/%s.rda", threshold, g))
      mask<-raster("../../Raster/mask_index.tif")
      df$mask_index<-raster::extract(mask, df[, c("x", "y")])
      where_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, mask_index, extinct_year)
      where_extinct<-where_extinct%>%dplyr::group_by(group, GCM, SSP, mask_index, extinct_year)%>%
        dplyr::summarise(n_sp=n())
      
      mask_p<-data.frame(rasterToPoints(mask))
      
      where_extinct<-left_join(mask_p, where_extinct, by="mask_index")
      sp_richness<-readRDS(sprintf("../../Objects/Diversity_%d/%s/EC-Earth3-Veg_SSP119_1/indices_df.rda", 
                                   threshold, g))
      
      sp_richness<-sp_richness[["2020"]][["species.richness"]]
      
      where_extinct<-left_join(where_extinct, sp_richness, by=c("mask_index"="index"))
      
      where_extinct$extinct_ratio<-where_extinct$n_sp/where_extinct$metric
      where_extinct<-where_extinct%>%dplyr::filter(extinct_ratio<=1)
      if (F){
        df%>%dplyr::filter(mask_index==1335)
        xx<-readRDS("../../Objects/IUCN_Distribution/Amphibians/Rana_amurensis.rda")
        xx$index<-extract(mask, xx[, c("x", "y")])
        xx%>%dplyr::filter(index==1335)
        
        xx<-readRDS("../../Objects/IUCN_Distribution/Amphibians/Salamandrella_keyserlingii.rda")
        xx$index<-extract(mask, xx[, c("x", "y")])
        xx%>%dplyr::filter(index==1335)
      }
      where_extinct_se<-where_extinct%>%dplyr::group_by(mask_index, group, SSP)%>%
        dplyr::summarise(mean_extinct_ratio=mean(extinct_ratio),
                         mean_n_extinct=mean(n_sp),
                         mean_species_richness=mean(metric))
      where_extinct_se<-where_extinct_se%>%dplyr::filter(!is.na(group))
      saveRDS(where_extinct_se, sprintf("../../Objects/when_where_extinction_%d/where_extinct_%s.rda", 
                                        threshold, g))
    }
  }
}
g<-"Amphibians"
mask<-raster("../../Raster/mask_index.tif")

SSP_i<-"SSP585"
for (threshold in c(1, 5)){
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    print(paste(g, threshold))
    where_extinct<-readRDS(sprintf("../../Objects/when_where_extinction_%d/where_extinct_%s.rda", threshold, g))
    png(sprintf("../../Figures/when_where_extinction_%d/where/%s_extinct_ratio.png", threshold, g), width=1000, height=200)
    par(mfrow=c(1, 3))
    for (SSP_i in c("SSP119", "SSP245", "SSP585")){
      item<-where_extinct%>%dplyr::filter(SSP==SSP_i)
      p_mask<-data.frame(rasterToPoints(mask))
      p_mask<-left_join(p_mask, item, by="mask_index")
      r<-mask
      values(r)[!is.na(values(mask))]<-p_mask$mean_extinct_ratio
      
      plot(r, main=paste(g, SSP_i))
      
      writeRaster(r, sprintf("../../Figures/when_where_extinction_%d/where/%s_%s_extinct_ratio.tif", threshold, g, SSP_i), overwrite=T)
    }
    dev.off()
    
    png(sprintf("../../Figures/when_where_extinction_%d/where/%s_n_extinct.png", threshold, g), width=1000, height=200)
    par(mfrow=c(1, 3))
    for (SSP_i in c("SSP119", "SSP245", "SSP585")){
      item<-where_extinct%>%dplyr::filter(SSP==SSP_i)
      p_mask<-data.frame(rasterToPoints(mask))
      p_mask<-left_join(p_mask, item, by="mask_index")
      r<-mask
      values(r)[!is.na(values(mask))]<-p_mask$mean_n_extinct
      
      plot(r, main=paste(g, SSP_i))
      
      writeRaster(r, sprintf("../../Figures/when_where_extinction_%d/where/%s_%s_n_extinct.tif", threshold, g, SSP_i), overwrite=T)
    }
    dev.off()
  }
}

if (F){
  mask<-raster("../../Raster/mask_index.tif")
  mask_p<-data.frame(rasterToPoints(mask))
  
  ratio_final<-NULL
  n_ext_final<-NULL
  for (threshold in c(1, 5)){
    for (SSP_i in c("SSP119", "SSP245", "SSP585")){
      ratio_all<-NULL
      n_ext_all<-NULL
      for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
        print(paste(threshold, SSP_i, g))
        ratio<-raster(sprintf("../../Figures/when_where_extinction_%d/where/%s_%s_extinct_ratio.tif", threshold, g, SSP_i))
        ratio<-data.frame(rasterToPoints(ratio))
        colnames(ratio)[3]<-"V"
        ratio<-left_join(mask_p, ratio, by=c("x", "y"))
        ratio[which(is.na(ratio$V)), "V"]<-0
        ratio_all<-bind_dplyr(ratio_all, ratio)
        n_ext<-raster(sprintf("../../Figures/when_where_extinction_%d/where/%s_%s_n_extinct.tif", threshold, g, SSP_i))
        n_ext<-data.frame(rasterToPoints(n_ext))
        colnames(n_ext)[3]<-"V"
        n_ext<-left_join(mask_p, n_ext, by=c("x", "y"))
        n_ext[which(is.na(n_ext$V)), "V"]<-0
        n_ext_all<-bind_dplyr(n_ext_all, n_ext)
      }
      ratio_all_se<-ratio_all%>%dplyr::group_by(x, y, mask_index)%>%
        dplyr::summarise(mean_V=mean(V, na.rm=T))
      ratio_all_se$SSP<-SSP_i
      ratio_all_se$threshold<-threshold
      
      n_ext_all_se<-n_ext_all%>%dplyr::group_by(x, y, mask_index)%>%
        dplyr::summarise(sum_V=sum(V, na.rm=T))
      n_ext_all_se$SSP<-SSP_i
      n_ext_all_se$threshold<-threshold
      
      ratio_final<-bind_dplyr(ratio_final, ratio_all_se)
      n_ext_final<-bind_dplyr(n_ext_final, n_ext_all_se)
      
    }
  }
  saveRDS(ratio_final, "../../Figures/when_where_extinction_all/ratio_final.rda")
  saveRDS(n_ext_final, "../../Figures/when_where_extinction_all/n_ext_final.rda")
  
}
myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))
ratio_final<-readRDS("../../Figures/when_where_extinction_all/ratio_final.rda")
ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$threshold)
colors<-myPalette(max(p_full_se$MEAN_V, na.rm=T))

p<-ggplot(ratio_final)+geom_tile(aes(x=x, y=y, fill=mean_V))+
  facet_wrap(~label, ncol=2)+
  scale_fill_gradient(low=color_two_map[2], high=color_two_map[1])+
  ggtitle("Extinction proportion")+
  labs(fill = "Extinction Proportion")+
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
ggsave(p, filename="../../Figures/when_where_extinction_all/ratio_final.png", width=8, height=8)
ggsave(p, filename="../../Figures/when_where_extinction_all/ratio_final.pdf", width=8, height=8)

n_ext_final<-readRDS("../../Figures/when_where_extinction_all/n_ext_final.rda")
n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$threshold)

p<-ggplot(n_ext_final)+geom_tile(aes(x=x, y=y, fill=sum_V))+
  facet_wrap(~label, ncol=2)+
  scale_fill_gradient(low=color_two_map[2], high=color_two_map[1])+
  ggtitle("Number of extinct species")+
  labs(fill = "Number of extinct species")+
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
ggsave(p, filename="../../Figures/when_where_extinction_all/n_ext_final.png", width=8, height=8)
ggsave(p, filename="../../Figures/when_where_extinction_all/n_ext_final.pdf", width=8, height=8)
