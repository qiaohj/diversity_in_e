library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
g<-"Reptiles"
threshold<-5
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))

groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")


SSPs<-c("SSP119", "SSP245", "SSP585")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")

diectory<-sprintf("Diversity_%d", threshold)
SSP_i<-SSPs[1]
GCM_i<-GCMs[1]
M_i<-0
y<-2020
if (F){
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
    saveRDS(p_full, file=sprintf("../../Figures/Diversity_5/species.richness.%s.rda", g))
  }
}
myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))

for (g in groups){
  print(g)
  p_full<-readRDS(sprintf("../../Figures/Diversity_5/species.richness.%s.rda", g))
  p_full_se<-p_full%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP)%>%
    dplyr::summarise(MEAN_V=mean(V, na.rm=T))
  colors<-myPalette(max(p_full_se$MEAN_V, na.rm=T))
  d1<-p_full_se%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119"))
  p1<-ggplot(p_full_se%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")))+
    geom_tile(aes(x=x, y=y, fill=MEAN_V))+
    scale_fill_gradientn(colors=colors[c(1:max(d1$MEAN_V, na.rm=T))])+
    ggtitle(g)+
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
    ggtitle("With dispersal")+
    facet_wrap(~SSP)
  g_legend<-get_legend(p3)
  p3<-p3+map_theme
  p<-ggarrange(p1, ggarrange(p2, p3, ncol=1, nrow=2), 
               common.legend = T, legend = "right", legend.grob = g_legend, widths = c(6, 8))  
  p
  ggsave(p, filename=sprintf("../../Figures/Diversity_5/MAP.%s.species.richness.png",g),
         width=14, height=5)
}
