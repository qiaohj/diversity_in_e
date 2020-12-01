library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
g<-"Reptiles"
threshold<-as.numeric(args[1])
if (is.na(threshold)){
  threshold<-1
}
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))

groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")


SSPs<-c("SSP119", "SSP245", "SSP585")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")


SSP_i<-SSPs[1]
GCM_i<-GCMs[1]
M_i<-0
y<-2020
if (F){
  for (threshold in c(1, 5)){
    diectory<-sprintf("Diversity_%d", threshold)
    for (g in groups){
      p_full<-NULL
      for (SSP_i in SSPs){
        for (y in c(2020, 2100)){
          for (M_i in c(0, 1, 2)){
            for (GCM_i in GCMs){
              raster_f<-sprintf("../../Figures/%s/%s/%s/%s_%s_%d/%d.tif", 
                                diectory, g, "species.richness", GCM_i, SSP_i, M_i, y)
              print(raster_f)
              r<-raster(raster_f)
              p<-data.frame(rasterToPoints(r))
              colnames(p)[3]<-"V"
              p<-left_join(mask_p, p, by=c("x", "y"))
              p$GCM<-GCM_i
              p$M<-M_i
              p$YEAR<-y
              p$SSP<-SSP_i
              
              p_full<-bind_dplyr(p_full, p)
            }
            
          }
        }
      }
      saveRDS(p_full, file=sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", threshold, g))
    }
  }
}
myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))
for (threshold in c(1, 5)){
  p_full_all<-NULL
  for (g in groups){
    print(paste(g, threshold))
    p_full<-readRDS(sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", threshold, g))
    p_full$group<-g
    p_full_all<-bind_dplyr(p_full_all, p_full)
    p_full_se<-p_full%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP)%>%
      dplyr::summarise(MEAN_V=mean(V, na.rm=T))
    colors<-myPalette(max(p_full_se$MEAN_V, na.rm=T))
    d1<-p_full_se%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119"))
    p1<-ggplot(p_full_sthresholde%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")))+
      geom_tile(aes(x=x, y=y, fill=MEAN_V))+
      scale_fill_gradientn(colors=colors[c(1:max(d1$MEAN_V, na.rm=T))])+
      ggtitle(paste(g, "Exposure year:", threshold))+
      map_theme
    p1
    
    d2<-p_full_se%>%dplyr::filter((YEAR==2100)&(M==0))
    p2<-ggplot(d2)+
      geom_tile(aes(x=x, y=y, fill=MEAN_V))+
      scale_fill_gradientn(colors=colors[c(1,max(d2$MEAN_V, na.rm=T))])+
      ggtitle("Without dispersal")+
      facet_wrap(~SSP)+
      map_theme
    p2
    
    d3<-p_full_se%>%dplyr::filter((YEAR==2100)&(M==1))
    p3<-ggplot(d3)+
      geom_tile(aes(x=x, y=y, fill=MEAN_V))+
      scale_fill_gradientn(colors=colors[c(1,max(d3$MEAN_V, na.rm=T))])+
      labs(fill = "Richness")+
      ggtitle("With 1 dispersal")+
      facet_wrap(~SSP)
    g_legend<-get_legend(p3)
    p3<-p3+map_theme
    
    d4<-p_full_se%>%dplyr::filter((YEAR==2100)&(M==2))
    p4<-ggplot(d4)+
      geom_tile(aes(x=x, y=y, fill=MEAN_V))+
      scale_fill_gradientn(colors=colors[c(1,max(d3$MEAN_V, na.rm=T))])+
      labs(fill = "Richness")+
      ggtitle("With 2 dispersal")+
      facet_wrap(~SSP)
    g_legend<-get_legend(p4)
    p4<-p4+map_theme
    
    
    p<-ggarrange(p1, ggarrange(p2, p3, p4, ncol=1, nrow=3), 
                 common.legend = T, legend = "right", legend.grob = g_legend, widths = c(8, 7))  
    p
    ggsave(p, filename=sprintf("../../Figures/Diversity_%d/MAP.%s.species.richness.png", threshold, g),
           width=14, height=5)
    
    ggsave(p, filename=sprintf("../../Figures/Diversity_%d/MAP.%s.species.richness.pdf", threshold, g),
           width=14, height=5)
  }
  
}

if (F){
  
  threshold<-5
  p_full_all<-NULL
  for (threshold in c(1, 5)){
    for (g in groups){
      print(paste(g, threshold))
      p_full<-readRDS(sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", threshold, g))
      p_full_se<-p_full%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP)%>%
        dplyr::summarise(MEAN_V=mean(V, na.rm=T))
      p_full_se$group<-g
      p_full_se$threshold<-threshold
      p_full_all<-bind_dplyr(p_full_all, p_full_se)
      
    }
  }
  p_full_all_se<-p_full_all%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP, threshold)%>%
    dplyr::summarise(SUM_V=sum(MEAN_V, na.rm=T))
  saveRDS(p_full_all_se, "../../Figures/Diversity_all/species.richness.rda")
}

p_full_all_se<-readRDS("../../Figures/Diversity_all/species.richness.rda")


for (tttt in c(1, 5)){
  print(tttt)
  p_full_all_se_1<-p_full_all_se%>%dplyr::filter(threshold==tttt)
  colors<-myPalette(max(p_full_all_se_1$SUM_V, na.rm=T))
  d1<-p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119"))
  p1<-ggplot(p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")))+
    geom_tile(aes(x=x, y=y, fill=SUM_V))+
    scale_fill_gradientn(colors=colors[c(1:max(d1$SUM_V, na.rm=T))])+
    ggtitle(paste("Species richness, Exposure year:", tttt))+
    map_theme
  p1
  
  d2<-p_full_all_se_1%>%dplyr::filter((YEAR==2100)&(M==0))
  p2<-ggplot(d2)+
    geom_tile(aes(x=x, y=y, fill=SUM_V))+
    scale_fill_gradientn(colors=colors[c(1,max(d2$SUM_V, na.rm=T))])+
    ggtitle("Without dispersal")+
    facet_wrap(~SSP)+
    map_theme
  p2
  
  d3<-p_full_all_se_1%>%dplyr::filter((YEAR==2100)&(M==1))
  p3<-ggplot(d3)+
    geom_tile(aes(x=x, y=y, fill=SUM_V))+
    scale_fill_gradientn(colors=colors[c(1,max(d3$SUM_V, na.rm=T))])+
    labs(fill = "Richness")+
    ggtitle("With 1 dispersal")+
    facet_wrap(~SSP)
  g_legend<-get_legend(p3)
  p3<-p3+map_theme
  
  d4<-p_full_all_se_1%>%dplyr::filter((YEAR==2100)&(M==2))
  p4<-ggplot(d4)+
    geom_tile(aes(x=x, y=y, fill=SUM_V))+
    scale_fill_gradientn(colors=colors[c(1,max(d3$SUM_V, na.rm=T))])+
    labs(fill = "Richness")+
    ggtitle("With 2 dispersal")+
    facet_wrap(~SSP)
  g_legend<-get_legend(p4)
  p4<-p4+map_theme
  
  
  p<-ggarrange(p1, ggarrange(p2, p3, p4, ncol=1, nrow=3), 
               common.legend = T, legend = "right", legend.grob = g_legend, widths = c(8, 7))  
  p
  ggsave(p, filename=sprintf("../../Figures/Diversity_all/MAP.species.richness.threshold.%d.png", tttt),
         width=14, height=5)
  
  ggsave(p, filename=sprintf("../../Figures/Diversity_all/MAP.species.richness.threshold.%d.pdf", tttt),
         width=14, height=5)
}
