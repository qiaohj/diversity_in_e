library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Rmisc)
library(rgdal)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
g<-"Mammals"
exposure<-as.numeric(args[1])
if (is.na(exposure)){
  exposure<-0
}
mask<-raster("../../Raster/mask_100km.tif")
mask_p<-data.frame(rasterToPoints(mask))

groups<-c("Birds", "Mammals")


SSPs<-c("SSP119", "SSP245", "SSP585")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")


SSP_i<-SSPs[1]
GCM_i<-GCMs[1]
M_i<-0
y<-2020
if (F){
  for (exposure in c(0, 5)){
    for (M_i in c(0, 1)){
      diectory<-sprintf("Diversity_exposure_%d_dispersal_%d", exposure, M_i)
      for (g in groups){
        p_full<-NULL
        for (SSP_i in SSPs){
          for (y in c(2020:2100)){
            
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
        saveRDS(p_full, file=
                  sprintf("../../Figures/Diversity_exposure_%d_dispersal_%d/species.richness.%s.rda", 
                          exposure, M_i, g))
      }
      
    }
  }
}
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))

for (exposure in c(0, 5)){
  for (M_i in c(0, 1)){
    p_full_all<-NULL
    for (g in groups){
      print(paste(g, exposure))
      p_full<-readRDS(sprintf("../../Figures/Diversity_exposure_%d_dispersal_%d/species.richness.%s.rda", 
                              exposure, M_i, g))
      p_full$group<-g
      p_full_all<-bind_dplyr(p_full_all, p_full)
      p_full_se<-p_full%>%dplyr::group_by(x, y, mask_100km, M, YEAR, SSP)%>%
        dplyr::summarise(MEAN_V=mean(V, na.rm=T))
      
      dd_df<-p_full_se%>%dplyr::filter((YEAR %in% c(2020, 2100))&(M %in% c(0, 1)))
      max_v<-max(dd_df$MEAN_V, na.rm=T)
      colors<-myPalette(max_v)
      d1<-p_full_se%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")&(!is.na(MEAN_V)))
      d1$log_mean_v<-log2(d1$MEAN_V)
      exposure_label<-ifelse(exposure==0, "no exposure", "5-year exposure")
      dispersal_label<-ifelse(M_i==0, "without dispersal", "with dispersal")
      p1<-ggplot(d1)+
        geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
        geom_tile(aes(x=x, y=y, fill=MEAN_V))+
        scale_fill_gradientn(colors=colors[c(1,max(d1$MEAN_V, na.rm=T))])+
        #ggtitle(paste("2020", g, exposure_label, dispersal_label, sep=", "))+
        ggtitle(g)+
        map_theme
      
      p1
      
      
      d2<-p_full_se%>%dplyr::filter((YEAR==2100)&(!is.na(MEAN_V)))
      p2<-ggplot(d2)+
        geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
        geom_tile(aes(x=x, y=y, fill=MEAN_V))+
        scale_fill_gradientn(colors=colors[c(1,max(d2$MEAN_V, na.rm=T))])+
        ggtitle("Without dispersal")+
        facet_wrap(~SSP)+
        map_theme
      p2
      
      d3<-p_full_se%>%dplyr::filter((YEAR==2100)&(M==1)&(!is.na(MEAN_V)))
      p3<-ggplot(d3)+
        geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
        geom_tile(aes(x=x, y=y, fill=MEAN_V))+
        scale_fill_gradientn(colors=colors[c(1,max(d3$MEAN_V, na.rm=T))])+
        labs(fill = "Richness")+
        ggtitle("With dispersal")+
        facet_wrap(~SSP)
      g_legend<-get_legend(p3)
      p3<-p3+map_theme
      
      if (F){
        d4<-p_full_se%>%dplyr::filter((YEAR==2100)&(M==2)&(!is.na(MEAN_V)))
        p4<-ggplot(d4)+
          geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
          geom_tile(aes(x=x, y=y, fill=MEAN_V))+
          scale_fill_gradientn(colors=colors[c(1,max(d3$MEAN_V, na.rm=T))])+
          labs(fill = "Richness")+
          ggtitle("With 2 dispersal")+
          facet_wrap(~SSP)
        g_legend<-get_legend(p4)
        p4<-p4+map_theme
      }
      
      p<-ggarrange(p1, ggarrange(p2, p3, ncol=1, nrow=2), 
                   common.legend = T, legend = "right", legend.grob = g_legend, widths = c(6, 7))  
      p
      ggsave(p, filename=sprintf("../../Figures/Diversity_%d/MAP.%s.species.richness.png", exposure, g),
             width=14, height=5)
      
      ggsave(p, filename=sprintf("../../Figures/Diversity_%d/MAP.%s.species.richness.pdf", exposure, g),
             width=14, height=5)
    }
  }
}

if (F){
  
  exposure<-5
  p_full_all<-NULL
  for (exposure in c(1, 5)){
    for (g in groups){
      print(paste(g, exposure))
      p_full<-readRDS(sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", exposure, g))
      p_full_se<-p_full%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP)%>%
        dplyr::summarise(MEAN_V=mean(V, na.rm=T))
      p_full_se$group<-g
      p_full_se$exposure<-exposure
      p_full_all<-bind_dplyr(p_full_all, p_full_se)
      
    }
  }
  p_full_all_se<-p_full_all%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP, exposure)%>%
    dplyr::summarise(SUM_V=sum(MEAN_V, na.rm=T))
  saveRDS(p_full_all_se, "../../Figures/Diversity_all/species.richness.rda")
  p_full_all_se_by_group<-p_full_all%>%dplyr::group_by(x, y, mask_index, M, YEAR, SSP, exposure, group)%>%
    dplyr::summarise(SUM_V=sum(MEAN_V, na.rm=T))
  saveRDS(p_full_all_se_by_group, "../../Figures/Diversity_all/species.richness_by_group.rda")
}

p_full_all_se<-readRDS("../../Figures/Diversity_all/species.richness.rda")
p_full_all_se$exposure<-" no exposure"
p_full_all_se[which(p_full_all_se$exposure==5), "exposure"]<-"5-year exposure"

p_full_all_summary<-p_full_all_se%>%dplyr::filter(YEAR %in% c(2020, 2100))%>%
  dplyr::group_by(YEAR, M, SSP, exposure, exposure)%>%
  dplyr::summarise(max_V=round(max(SUM_V)),
                   min_V=round(min(SUM_V)))
p_full_all_summary<-p_full_all_summary%>%dplyr::filter(M!=2)
write.csv(p_full_all_summary, "../../Figures/Diversity_all/richness.csv")
da=0
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
for (da in c(0, 1)){
  print(da)
  p_full_all_se_1<-p_full_all_se%>%dplyr::filter(M==da)
  colors<-myPalette(max(p_full_all_se_1$SUM_V, na.rm=T))
  d1<-p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")&(!is.na(SUM_V)))
  p1<-ggplot(p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")))+
    geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
    geom_tile(aes(x=x, y=y, fill=SUM_V), color=NA)+
    scale_fill_gradientn(colors=colors[c(1:max(d1$SUM_V, na.rm=T))])+
    scale_color_gradientn(colors=colors[c(1:max(d1$SUM_V, na.rm=T))])+
    ggtitle(ifelse(da==0, "Species richness, no dispersal", "Species richness, with dispersal"))+
    map_theme
  p1
  
  d2<-p_full_all_se_1%>%dplyr::filter((YEAR==2100)&(!is.na(SUM_V)))
  p2<-ggplot(d2)+
    geom_tile(aes(x=x, y=y, fill=SUM_V), color=NA)+
    scale_fill_gradientn(colors=colors[c(1,max(d2$SUM_V, na.rm=T))])+
    labs(fill = "Richness")+
    facet_grid(exposure~SSP)
  g_legend<-get_legend(p2)
  p2<-p2+  map_theme
  p2
  
  # p<-ggarrange(p1, ggarrange(p2, p3, p4, ncol=1, nrow=3), 
  #              common.legend = T, legend = "right", legend.grob = g_legend, widths = c(8, 7))  
  p<-ggarrange(p1, p2, 
               common.legend = T, legend = "right", legend.grob = g_legend, widths = c(6, 8))  
  p
  ggsave(p, filename=sprintf("../../Figures/Diversity_all/MAP.species.richness.dispersal.%d.png", da),
         width=14, height=4)
  
  ggsave(p, filename=sprintf("../../Figures/Diversity_all/MAP.species.richness.dispersal.%d.pdf", da),
         width=14, height=4)
}


p_full_all_se_by_group<-readRDS("../../Figures/Diversity_all/species.richness_by_group.rda")
p_full_all_se_by_group$exposure<-" no exposure"
p_full_all_se_by_group[which(p_full_all_se_by_group$exposure==5), "exposure"]<-"5-year exposure"

p_full_all_summary<-p_full_all_se_by_group%>%dplyr::filter(YEAR %in% c(2020, 2100))%>%
  dplyr::group_by(YEAR, M, SSP, exposure, exposure, group)%>%
  dplyr::summarise(max_V=round(max(SUM_V)),
                   min_V=round(min(SUM_V)))
p_full_all_summary<-p_full_all_summary%>%dplyr::filter(M!=2)
write.csv(p_full_all_summary, "../../Figures/Diversity_all/richness_by_group.csv")
da=0
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))

for (da in c(0, 1)){
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    print(paste(da, g))
    p_full_all_se_1<-p_full_all_se_by_group%>%dplyr::filter((M==da)&(group==g))
    colors<-myPalette(max(p_full_all_se_1$SUM_V, na.rm=T))
    d1<-p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")&(!is.na(SUM_V)))
    p1<-ggplot(p_full_all_se_1%>%dplyr::filter((YEAR==2020)&(SSP=="SSP119")))+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=SUM_V), color=NA)+
      scale_fill_gradientn(colors=colors[c(1:max(d1$SUM_V, na.rm=T))])+
      scale_color_gradientn(colors=colors[c(1:max(d1$SUM_V, na.rm=T))])+
      ggtitle(paste(g, "-",  
                    ifelse(da==0, "Species richness, no dispersal", "Species richness, with dispersal")))+
      map_theme
    p1
    
    d2<-p_full_all_se_1%>%dplyr::filter((YEAR==2100)&(!is.na(SUM_V)))
    p2<-ggplot(d2)+
      geom_tile(aes(x=x, y=y, fill=SUM_V), color=NA)+
      scale_fill_gradientn(colors=colors[c(1,max(d2$SUM_V, na.rm=T))])+
      labs(fill = "Richness")+
      facet_grid(exposure~SSP)
    g_legend<-get_legend(p2)
    p2<-p2+  map_theme
    p2
    
    # p<-ggarrange(p1, ggarrange(p2, p3, p4, ncol=1, nrow=3), 
    #              common.legend = T, legend = "right", legend.grob = g_legend, widths = c(8, 7))  
    p<-ggarrange(p1, p2, 
                 common.legend = T, legend = "right", legend.grob = g_legend, widths = c(6, 8))  
    p
    ggsave(p, filename=
             sprintf("../../Figures/Diversity_all/MAP.species.richness.dispersal.%d_%s.png", da, g),
           width=14, height=4)
    
    ggsave(p, filename=
             sprintf("../../Figures/Diversity_all/MAP.species.richness.dispersal.%d_%s.pdf", da, g),
           width=14, height=4)
  }
}

if (F){
  df1<-readRDS(sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", 1, "Birds"))
  df2<-readRDS(sprintf("../../Figures/Diversity_%d/species.richness.%s.rda", 5, "Birds"))
  df1<-df1%>%dplyr::filter((YEAR==2100)&(!is.na(V)))
  df2<-df2%>%dplyr::filter((YEAR==2100)&(!is.na(V)))
  df_xx<-inner_join(df1, df2, by=c("x", "y", "mask_index", "GCM", "M", "YEAR", "SSP"))
  df_xx$differ<-df_xx$V.x-df_xx$V.y
  range(df_xx$differ)
}
