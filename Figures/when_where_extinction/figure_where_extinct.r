library(dplyr)
library(raster)
library(gglpot2)
library(Rmisc)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")

if (F){
  g<-"Birds"
  for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
    print(g)
    df<-readRDS(sprintf("../../Objects/when_where_extinction_5/%s.rda", g))
    mask<-raster("../../Raster/mask_index.tif")
    df$mask_index<-raster::extract(mask, df[, c("x", "y")])
    where_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, mask_index, extinct_year)
    where_extinct<-where_extinct%>%dplyr::group_by(group, GCM, SSP, mask_index, extinct_year)%>%
      dplyr::summarise(n_sp=n())
    
    mask_p<-data.frame(rasterToPoints(mask))
    
    where_extinct<-left_join(mask_p, where_extinct, by="mask_index")
    sp_richness<-readRDS(sprintf("../../Objects/Diversity_5/%s/EC-Earth3-Veg_SSP119_1/indices_df.rda", group))
    
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
    saveRDS(where_extinct_se, sprintf("../../Objects/when_where_extinction_5/where_extinct_%s.rda", g))
  }
}
g<-"Amphibians"
mask<-raster("../../Raster/mask_index.tif")

SSP_i<-"SSP585"
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  where_extinct<-readRDS(sprintf("../../Objects/when_where_extinction_5/where_extinct_%s.rda", g))
  png(sprintf("../../Figures/when_where_extinction_5/where/%s.png", g), width=1000, height=200)
  par(mfrow=c(1, 3))
  for (SSP_i in c("SSP119", "SSP245", "SSP585")){
    item<-where_extinct%>%dplyr::filter(SSP==SSP_i)
    p_mask<-data.frame(rasterToPoints(mask))
    p_mask<-left_join(p_mask, item, by="mask_index")
    r<-mask
    values(r)[!is.na(values(mask))]<-p_mask$mean_extinct_ratio
    
    plot(r, main=paste(g, SSP_i))
    
    writeRaster(r, sprintf("../../Figures/when_where_extinction_5/where/%s_%s.tif", g, SSP_i), overwrite=T)
  }
  dev.off()
  
  png(sprintf("../../Figures/when_where_extinction_5/where/%s_n_extinct.png", g), width=1000, height=200)
  par(mfrow=c(1, 3))
  for (SSP_i in c("SSP119", "SSP245", "SSP585")){
    item<-where_extinct%>%dplyr::filter(SSP==SSP_i)
    p_mask<-data.frame(rasterToPoints(mask))
    p_mask<-left_join(p_mask, item, by="mask_index")
    r<-mask
    values(r)[!is.na(values(mask))]<-p_mask$mean_n_extinct
    
    plot(r, main=paste(g, SSP_i))
    
    writeRaster(r, sprintf("../../Figures/when_where_extinction_5/where/%s_%s_n_extinct.tif", g, SSP_i), overwrite=T)
  }
  dev.off()
}

