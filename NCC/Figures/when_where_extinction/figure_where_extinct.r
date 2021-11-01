library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(scales)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

if (F){
  exposure=0
  g<-"Mammals"
  for (exposure in c(0, 5)){
    for (g in c("Birds", "Mammals")){
      print(paste(g, exposure))
      df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
      SP_N_CELL<-df%>%dplyr::group_by(group, sp, GCM, SSP, extinct_year, dispersal)%>%
        dplyr::summarise(N_CELL=n_distinct(x, y))
      SP_N_CELL<-SP_N_CELL%>%ungroup()%>%dplyr::distinct(group, sp, N_CELL)
      df_bak<-inner_join(df, SP_N_CELL, by=c("group", "sp"))
      
      mask<-raster("../../Raster/mask_100km.tif")
      df$mask_100km<-raster::extract(mask, df[, c("x", "y")])
      where_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, mask_100km, dispersal)
      where_extinct<-where_extinct%>%dplyr::group_by(group, GCM, SSP, mask_100km, dispersal)%>%
        dplyr::summarise(n_sp=n())
      
      mask_p<-data.frame(rasterToPoints(mask))
      
      where_extinct<-left_join(mask_p, where_extinct, by="mask_100km")
      
      dis=1
      where_extinct_all<-NULL
      for (dis in c(0, 1)){
        where_extinct_item<-where_extinct%>%dplyr::filter(dispersal==dis)
        sp_richness<-readRDS(sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d/%s/EC-Earth3-Veg_SSP119/indices_df.rda", 
                                     exposure, dis, g))
        
        sp_richness<-sp_richness[["2020"]][["species.richness"]]
        
        where_extinct_item<-left_join(where_extinct_item, sp_richness, by=c("mask_100km"="index"))
        
        where_extinct_item$extinct_ratio<-where_extinct_item$n_sp/where_extinct_item$metric
        #where_extinct<-where_extinct%>%dplyr::filter(extinct_ratio<=1)
        if (F){
          df%>%dplyr::filter(mask_index==1335)
          xx<-readRDS("../../Objects/IUCN_Distribution/Amphibians/Rana_amurensis.rda")
          xx$index<-extract(mask, xx[, c("x", "y")])
          xx%>%dplyr::filter(index==1335)
          
          xx<-readRDS("../../Objects/IUCN_Distribution/Amphibians/Salamandrella_keyserlingii.rda")
          xx$index<-extract(mask, xx[, c("x", "y")])
          xx%>%dplyr::filter(index==1335)
        }
        where_extinct_se<-where_extinct_item%>%dplyr::group_by(mask_100km, group, SSP, dispersal)%>%
          dplyr::summarise(mean_extinct_ratio=mean(extinct_ratio),
                           mean_n_extinct=mean(n_sp),
                           mean_species_richness=mean(metric))
        where_extinct_se<-where_extinct_se%>%dplyr::filter(!is.na(group))
        where_extinct_all<-bind_dplyr(where_extinct_all, where_extinct_se)
        
      }
      saveRDS(where_extinct_all, sprintf("../../Objects/when_where_extinction_exposure_%d/where_extinct_%s.rda", 
                                        exposure, g))
    }
  }
}
g<-"Mammals"
mask<-raster("../../Raster/mask_100km_plot.tif")
exposure<-0
SSP_i<-"SSP585"
da=1

for (exposure in c(0, 5)){
  for (g in c("Birds", "Mammals")){
    for (ttt in c(2)){
      where_extinct<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/where_extinct_%s.rda",
                                     exposure, g))
      for (da in c(0:1)){
        print(paste(g, exposure, da))
        png(sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_extinct_ratio_dispersal_%d.png", 
                    exposure, g, da), width=1000, height=200)
        par(mfrow=c(1, 3))
        for (SSP_i in c("SSP119", "SSP245", "SSP585")){
          item<-where_extinct%>%dplyr::filter((SSP==SSP_i)&(dispersal==da))
          p_mask<-data.frame(rasterToPoints(mask))
          colnames(p_mask)[3]<-"mask_100km"
          p_mask<-left_join(p_mask, item, by="mask_100km")
          r<-mask
          values(r)[!is.na(values(mask))]<-p_mask$mean_extinct_ratio
          
          plot(r, main=paste(g, SSP_i, da))
          
          writeRaster(r, sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_%s_extinct_ratio_dispersal_%d.tif",
                                 exposure, g, SSP_i, da), overwrite=T)
        }
        dev.off()
        
        png(sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_n_extinct_dispersal_%d.png", 
                    exposure, g, da), width=1000, height=200)
        par(mfrow=c(1, 3))
        for (SSP_i in c("SSP119", "SSP245", "SSP585")){
          item<-where_extinct%>%dplyr::filter((SSP==SSP_i)&(dispersal==da))
          p_mask<-data.frame(rasterToPoints(mask))
          colnames(p_mask)[3]<-"mask_100km"
          p_mask<-left_join(p_mask, item, by="mask_100km")
          r<-mask
          values(r)[!is.na(values(mask))]<-p_mask$mean_n_extinct
          
          plot(r, main=paste(g, SSP_i, da))
          
          writeRaster(r, sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_%s_n_extinct_dispersal_%d.tif", 
                                 exposure, g, SSP_i, da), overwrite=T)
        }
        dev.off()
      }
    }
  }
}

if (F){
  mask<-raster("../../Raster/mask_100km.tif")
  mask_p<-data.frame(rasterToPoints(mask))
  
  ratio_final<-NULL
  n_ext_final<-NULL
  ratio_final_group<-NULL
  n_ext_final_group<-NULL
  
  for (exposure in c(0, 5)){
    for (SSP_i in c("SSP119", "SSP245", "SSP585")){
      for (da in c(0:1)){
        ratio_all<-NULL
        n_ext_all<-NULL
        for (g in c("Birds", "Mammals")){
          
          print(paste(exposure, SSP_i, g, da))
          ratio<-raster(sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_%s_extinct_ratio_dispersal_%d.tif", 
                                exposure, g, SSP_i, da))
          ratio<-data.frame(rasterToPoints(ratio))
          colnames(ratio)[3]<-"V"
          ratio<-left_join(mask_p, ratio, by=c("x", "y"))
          ratio[which(is.na(ratio$V)), "V"]<-0
          ratio$group<-g
          ratio_all<-bind_dplyr(ratio_all, ratio)
          n_ext<-raster(sprintf("../../Figures/when_where_extinction_exposure_%d/where/%s_%s_n_extinct_dispersal_%d.tif", 
                                exposure, g, SSP_i, da))
          n_ext<-data.frame(rasterToPoints(n_ext))
          colnames(n_ext)[3]<-"V"
          n_ext<-left_join(mask_p, n_ext, by=c("x", "y"))
          n_ext$group<-g
          
          n_ext[which(is.na(n_ext$V)), "V"]<-0
          n_ext_all<-bind_dplyr(n_ext_all, n_ext)
          
        }
        ratio_all_se<-ratio_all%>%ungroup()%>%dplyr::group_by(x, y, mask_100km)%>%
          dplyr::summarise(mean_V=mean(V, na.rm=T))
        ratio_all_se$SSP<-SSP_i
        ratio_all_se$exposure<-exposure
        ratio_all_se$dispersal<-da
        
        n_ext_all_se<-n_ext_all%>%ungroup()%>%dplyr::group_by(x, y, mask_100km)%>%
          dplyr::summarise(sum_V=sum(V, na.rm=T))
        n_ext_all_se$SSP<-SSP_i
        n_ext_all_se$exposure<-exposure
        n_ext_all_se$dispersal<-da
        
        ratio_final<-bind_dplyr(ratio_final, ratio_all_se)
        n_ext_final<-bind_dplyr(n_ext_final, n_ext_all_se)
        
        ratio_all_group_se<-ratio_all%>%ungroup()%>%dplyr::group_by(x, y, mask_100km, group)%>%
          dplyr::summarise(mean_V=mean(V, na.rm=T))
        ratio_all_group_se$SSP<-SSP_i
        ratio_all_group_se$exposure<-exposure
        ratio_all_group_se$dispersal<-da
        
        n_ext_all_group_se<-n_ext_all%>%ungroup()%>%dplyr::group_by(x, y, mask_100km, group)%>%
          dplyr::summarise(sum_V=sum(V, na.rm=T))
        n_ext_all_group_se$SSP<-SSP_i
        n_ext_all_group_se$exposure<-exposure
        n_ext_all_group_se$dispersal<-da
        
        ratio_final_group<-bind_dplyr(ratio_final_group, ratio_all_group_se)
        n_ext_final_group<-bind_dplyr(n_ext_final_group, n_ext_all_group_se)
      }
      
    }
    
  }
  saveRDS(ratio_final, "../../Figures/when_where_extinction_all/ratio_final.rda")
  saveRDS(n_ext_final, "../../Figures/when_where_extinction_all/n_ext_final.rda")
  saveRDS(ratio_final_group, "../../Figures/when_where_extinction_all/ratio_final_group.rda")
  saveRDS(n_ext_final_group, "../../Figures/when_where_extinction_all/n_ext_final_group.rda")
  
}
da=1
mask<-raster("../../Raster/mask_100km_plot.tif")
mask_p<-data.frame(rasterToPoints(mask))
colnames(mask_p)[3]<-"mask_100km"
for (da in c(0:1)){
  print(paste(da))
  if (da==0){
    title1<-"Extinction proportion (no dispersal)"
    title2<-"Number of extinct species (no dispersal)"
  }else{
    title1<-"Extinction proportion (with dispersal)"
    title2<-"Number of extinct species (with dispersal)"
  }
  
  myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))
  ratio_final<-readRDS("../../Figures/when_where_extinction_all/ratio_final.rda")
  
  
  ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$exposure)
  ratio_final<-ratio_final%>%dplyr::filter(mean_V>0)
  ratio_final<-data.frame(ratio_final)
  ratio_final$exposure<-ifelse(ratio_final$exposure==0, " no exposure", "5-year exposure")
  ratio_final[which(ratio_final$exposure==0), "label"]<-
    paste(as.vector(ratio_final[which(ratio_final$exposure==0), "SSP"]), " no exposure", sep=", ")
  ratio_final[which(ratio_final$exposure==5), "label"]<-
    paste(ratio_final[which(ratio_final$exposure==5), "SSP"], "5-year exposure", sep=", ")
  
  hist(ratio_final$mean_V)
  p<-ggplot(ratio_final)+
    geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
    geom_tile(aes(x=x, y=y, fill=mean_V))+
    facet_grid(exposure~SSP)+
    scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                        limits=c(0, 0.6), oob=squish,
                        breaks=seq(0, 0.6, by=0.1),
                        labels=c("0%", "10%", "20%", "30%", "40%", "50%",
                                 sprintf(">60%%, up to %.0f%%", max(ratio_final$mean_V)*100)))+
    ggtitle(title1)+
    labs(fill = "Extinction proportion")+
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
  ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/ratio_final_dispersal_%d.png", da), 
         width=9, height=4)
  ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/ratio_final_dispersal_%d.pdf", da),
         width=9, height=4)
  
  n_ext_final<-readRDS("../../Figures/when_where_extinction_all/n_ext_final.rda")
  n_ext_final<-n_ext_final%>%dplyr::filter((dispersal==da))
  n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$exposure)
  n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
  n_ext_final$exposure<-ifelse(n_ext_final$exposure==0, " no exposure", "5-year exposure")
  n_ext_final<-data.frame(n_ext_final)
  n_ext_final[which(n_ext_final$exposure==0), "label"]<-
    paste(as.vector(n_ext_final[which(n_ext_final$exposure==0), "SSP"]), " no exposure", sep=", ")
  n_ext_final[which(n_ext_final$exposure==5), "label"]<-
    paste(n_ext_final[which(n_ext_final$exposure==5), "SSP"], "5-year exposure", sep=", ")
  hist(n_ext_final$sum_V)
  p<-ggplot(n_ext_final)+
    geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
    geom_tile(aes(x=x, y=y, fill=sum_V))+
    facet_grid(exposure~SSP)+
    scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                        limits=c(0, 100), oob=squish,
                        breaks=c(0, 20, 40, 60, 80, 100),
                        labels=c("0", "20", "40", "60", "80", 
                                 sprintf(">100, up to %d", round(max(n_ext_final$sum_V)))))+
    ggtitle(title2)+
    labs(fill = "Number of \nextinct species")+
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
  ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/n_ext_final_dispersal_%d.png", da),
         width=9, height=4)
  ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/n_ext_final_dispersal_%d.pdf", da),
         width=9, height=4)
  
  
  
  
  
  
  
  ratio_final<-readRDS("../../Figures/when_where_extinction_all/ratio_final_group.rda")
  ratio_final<-ratio_final%>%dplyr::filter((dispersal==da))
  ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$exposure)
  ratio_final<-ratio_final%>%dplyr::filter(mean_V>0)
  ratio_final<-data.frame(ratio_final)
  ratio_final$exposure<-ifelse(ratio_final$exposure==0, " no exposure", "5-year exposure")
  ratio_final[which(ratio_final$exposure==0), "label"]<-
    paste(as.vector(ratio_final[which(ratio_final$exposure==0), "SSP"]), " no exposure", sep=", ")
  ratio_final[which(ratio_final$exposure==5), "label"]<-
    paste(ratio_final[which(ratio_final$exposure==5), "SSP"], "5-year exposure", sep=", ")
  
  hist(ratio_final$mean_V)
  
  n_ext_final<-readRDS("../../Figures/when_where_extinction_all/n_ext_final_group.rda")
  n_ext_final<-n_ext_final%>%dplyr::filter((dispersal==da))
  n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$exposure)
  n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
  n_ext_final$exposure<-ifelse(n_ext_final$exposure==0, " no exposure", "5-year exposure")
  n_ext_final<-data.frame(n_ext_final)
  n_ext_final[which(n_ext_final$exposure==0), "label"]<-
    paste(as.vector(n_ext_final[which(n_ext_final$exposure==0), "SSP"]), " no exposure", sep=", ")
  n_ext_final[which(n_ext_final$exposure==5), "label"]<-
    paste(n_ext_final[which(n_ext_final$exposure==5), "SSP"], "5-year exposure", sep=", ")
  hist(n_ext_final$sum_V)
  
  for (g in c("Birds", "Mammals")){
    ratio_final_item<-ratio_final%>%dplyr::filter(group==g)
    #hist(ratio_final_item$mean_V)
    ggg<-scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                             limits=c(0, 0.6), oob=squish,
                             breaks=seq(0, 0.6, by=0.1),
                             labels=c("0%", "10%", "20%", "30%", "40%", "50%",
                                      sprintf(">60%%, up to %.0f%%", max(ratio_final_item$mean_V)*100)))
    #if (g=="Birds"){
    #  ggg<-scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
    #                           limits=c(0, 0.4), oob=squish,
    #                           breaks=seq(0, 0.4, by=0.1),
    #                           labels=c("0%", "10%", "20%", "30%", 
    #                                    sprintf(">40%%, up to %.0f%%", max(ratio_final_item$mean_V)*100)))
    #}
    p<-ggplot(ratio_final_item)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=mean_V))+
      facet_grid(exposure~SSP)+
      ggg+
      ggtitle(paste(g, title1, sep=" - "))+
      labs(fill = "Extinction proportion")+
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
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/ratio_final_dispersal_%d_%s.png", da, g), 
           width=9, height=4)
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/ratio_final_dispersal_%d_%s.pdf", da, g),
           width=9, height=4)
    
    n_ext_final_item<-n_ext_final%>%dplyr::filter(group==g)
    
    p<-ggplot(n_ext_final_item)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=sum_V))+
      facet_grid(exposure~SSP)+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                          limits=c(0, 40), oob=squish,
                          breaks=c(0, 10, 20, 30, 40),
                          labels=c("0", "10", "20", "30",
                                   sprintf(">40, up to %d", round(max(n_ext_final_item$sum_V)))))+
      ggtitle(paste(title2, g, sep=" - "))+
      labs(fill = "Number of \nextinct species")+
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
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/n_ext_final_dispersal_%d_%s.png", da, g),
           width=9, height=4)
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_all/n_ext_final_dispersal_%d_%s.pdf", da, g),
           width=9, height=4)
  }
  
}
