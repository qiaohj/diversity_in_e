library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(scales)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

g<-"Amphibians"
mask<-raster("../../Raster/mask_index.tif")
threshold<-1
SSP_i<-"SSP585"
da=1
tttt=2
mask_p<-data.frame(rasterToPoints(mask))
if (T){
  #for (da in c(0:1)){
    da=0
    print(paste(da, tttt))
    if (da==0){
      title1<-"Extinct proportion (no dispersal)"
      title2<-"Number of extinct species (no dispersal)"
    }else{
      title1<-"Extinct proportion (with dispersal)"
      title2<-"Number of extinct species (with dispersal)"
    }
    
    myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))
    ratio_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/ratio_final_with_ttt.rda")
    
    ratio_final<-ratio_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$threshold)
    ratio_final<-ratio_final%>%dplyr::filter(mean_V>0)
    ratio_final<-data.frame(ratio_final)
    ratio_final$exposure<-ifelse(ratio_final$threshold==1, " no exposure", "5-year exposure")
    ratio_final[which(ratio_final$threshold==1), "label"]<-
      paste(as.vector(ratio_final[which(ratio_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    ratio_final[which(ratio_final$threshold==5), "label"]<-
      paste(ratio_final[which(ratio_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
    
    n_ext_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/n_ext_final_with_ttt.rda")
    n_ext_final<-n_ext_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$threshold)
    n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
    n_ext_final$exposure<-ifelse(n_ext_final$threshold==1, " no exposure", "5-year exposure")
    n_ext_final<-data.frame(n_ext_final)
    n_ext_final[which(n_ext_final$threshold==1), "label"]<-
      paste(as.vector(n_ext_final[which(n_ext_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    n_ext_final[which(n_ext_final$threshold==5), "label"]<-
      paste(n_ext_final[which(n_ext_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
    
    p_ratio_no_da<-ggplot(ratio_final)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=mean_V))+
      facet_grid(exposure~SSP, scale="free")+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                          limits=c(0, 0.6), oob=squish,
                          breaks=seq(0, 0.6, by=0.1),
                          labels=c("0%", "10%", "20%", "30%", "40%", "50%",
                                   sprintf(">60%%, up to %.0f%%", max(ratio_final$mean_V)*100)))+
      #ggtitle(title1)+
      labs(fill = "Extinct proportion")+
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
        panel.border = element_blank(),
        legend.position="bottom",
        legend.key.width=unit(0.8, "in")
      )
    p_ratio_no_da
    
    p_n_ext_no_da<-ggplot(n_ext_final)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=sum_V))+
      facet_grid(exposure~SSP)+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                          limits=c(0, 100), oob=squish,
                          breaks=c(0, 20, 40, 60, 80, 100),
                          labels=c("0", "20", "40", "60", "80", 
                                   sprintf(">100, up to %d", round(max(n_ext_final$sum_V)))))+
      #ggtitle(title2)+
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
        panel.border = element_blank(),
        legend.position="bottom",
        legend.key.width=unit(0.8, "in")
      )
    p_n_ext_no_da

    da=1
    print(paste(da, tttt))
    if (da==0){
      title1<-"Extinct proportion (no dispersal)"
      title2<-"Number of extinct species (no dispersal)"
    }else{
      title1<-"Extinct proportion (with dispersal)"
      title2<-"Number of extinct species (with dispersal)"
    }
    
    myPalette <- colorRampPalette(c(color_two_map[2], color_two_map[1]))
    ratio_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/ratio_final_with_ttt.rda")
    
    ratio_final<-ratio_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$threshold)
    ratio_final<-ratio_final%>%dplyr::filter(mean_V>0)
    ratio_final<-data.frame(ratio_final)
    ratio_final$exposure<-ifelse(ratio_final$threshold==1, " no exposure", "5-year exposure")
    ratio_final[which(ratio_final$threshold==1), "label"]<-
      paste(as.vector(ratio_final[which(ratio_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    ratio_final[which(ratio_final$threshold==5), "label"]<-
      paste(ratio_final[which(ratio_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
    
    n_ext_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/n_ext_final_with_ttt.rda")
    n_ext_final<-n_ext_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$threshold)
    n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
    n_ext_final$exposure<-ifelse(n_ext_final$threshold==1, " no exposure", "5-year exposure")
    n_ext_final<-data.frame(n_ext_final)
    n_ext_final[which(n_ext_final$threshold==1), "label"]<-
      paste(as.vector(n_ext_final[which(n_ext_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    n_ext_final[which(n_ext_final$threshold==5), "label"]<-
      paste(n_ext_final[which(n_ext_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
    
    p_ratio_with_da<-ggplot(ratio_final)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=mean_V))+
      facet_grid(exposure~SSP, scale="free")+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                          limits=c(0, 0.6), oob=squish,
                          breaks=seq(0, 0.6, by=0.1),
                          labels=c("0%", "10%", "20%", "30%", "40%", "50%",
                                   sprintf(">60%%, up to %.0f%%", max(ratio_final$mean_V)*100)))+
      #ggtitle(title1)+
      labs(fill = "Extinct proportion")+
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
        panel.border = element_blank(),
        legend.position="bottom",
        legend.key.width=unit(0.8, "in")
      )
    p_ratio_with_da
    
    p_n_ext_with_da<-ggplot(n_ext_final)+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(aes(x=x, y=y, fill=sum_V))+
      facet_grid(exposure~SSP)+
      scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                          limits=c(0, 100), oob=squish,
                          breaks=c(0, 20, 40, 60, 80, 100),
                          labels=c("0", "20", "40", "60", "80", 
                                   sprintf(">100, up to %d", round(max(n_ext_final$sum_V)))))+
      #ggtitle(title2)+
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
        panel.border = element_blank(),
        legend.position="bottom",
        legend.key.width=unit(0.8, "in")
      )
    p_n_ext_no_da_formatted<-p_n_ext_no_da
    
    p_ratio_no_da_formatted<-p_ratio_no_da+theme(strip.background.x = element_blank(),
                                                 strip.text.x = element_blank())
    p_ratio_with_da_formatted<-p_ratio_with_da+theme(strip.background.x = element_blank(),
                                                 strip.text.x = element_blank())
    p_n_ext_with_da_formatted<-p_n_ext_with_da+theme(strip.background.x = element_blank(),
                                                 strip.text.x = element_blank())
    
    pp<-ggarrange(p_n_ext_no_da_formatted,
                  p_ratio_no_da_formatted, 
                  p_n_ext_with_da_formatted, 
                  p_ratio_with_da_formatted, ncol=1, nrow=4, labels=c("(a)", "(b)", "(c)", "(d)"))
    
    ggsave(pp, 
           filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_all_ttt_%d.png", tttt), 
           width=7, height=13)
    
    ggsave(pp, 
           filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_all_ttt_%d.pdf", tttt), 
           width=7, height=13)
}   

all_p_list<-list()
for (da in c(0, 1)){
    ratio_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/ratio_final_group_with_ttt.rda")
    ratio_final<-ratio_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    ratio_final$label<-paste(ratio_final$SSP, "Exposure year:", ratio_final$threshold)
    ratio_final<-ratio_final%>%dplyr::filter(mean_V>0)
    ratio_final<-data.frame(ratio_final)
    ratio_final$exposure<-ifelse(ratio_final$threshold==1, " no exposure", "5-year exposure")
    ratio_final[which(ratio_final$threshold==1), "label"]<-
      paste(as.vector(ratio_final[which(ratio_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    ratio_final[which(ratio_final$threshold==5), "label"]<-
      paste(ratio_final[which(ratio_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
    
    
    n_ext_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/n_ext_final_group_with_ttt.rda")
    n_ext_final<-n_ext_final%>%dplyr::filter((dispersal==da)&(ttt==tttt))
    n_ext_final$label<-paste(n_ext_final$SSP, "Exposure year:", n_ext_final$threshold)
    n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
    n_ext_final$exposure<-ifelse(n_ext_final$threshold==1, " no exposure", "5-year exposure")
    n_ext_final<-data.frame(n_ext_final)
    n_ext_final[which(n_ext_final$threshold==1), "label"]<-
      paste(as.vector(n_ext_final[which(n_ext_final$threshold==1), "SSP"]), " no exposure", sep=", ")
    n_ext_final[which(n_ext_final$threshold==5), "label"]<-
      paste(n_ext_final[which(n_ext_final$threshold==5), "SSP"], "5-year exposure", sep=", ")
   
    
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      print(paste(g, da))
      ratio_final_item<-ratio_final%>%dplyr::filter(group==g)
      #hist(ratio_final_item$mean_V)
      ggg<-scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                               limits=c(0, 0.6), oob=squish,
                               breaks=seq(0, 0.6, by=0.1),
                               labels=c("0%", "10%", "20%", "30%", "40%", "50%",
                                        sprintf(">60%%, up to %.0f%%", max(ratio_final_item$mean_V)*100)))
      if (g=="Birds"){
        ggg<-scale_fill_gradient(low=color_two_map[1], high=color_two_map[2], 
                                 limits=c(0, 0.4), oob=squish,
                                 breaks=seq(0, 0.4, by=0.1),
                                 labels=c("0%", "10%", "20%", "30%", 
                                          sprintf(">40%%, up to %.0f%%", max(ratio_final_item$mean_V)*100)))
      }
      p1<-ggplot(ratio_final_item)+
        geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
        geom_tile(aes(x=x, y=y, fill=mean_V))+
        facet_grid(exposure~SSP)+
        ggg+
        #ggtitle(paste(g, title1, sep=" - "))+
        labs(fill = "Extinct proportion")+
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
          panel.border = element_blank(),
          legend.position="bottom",
          legend.key.width=unit(0.8,"in"),
          #strip.background.x = element_blank(),
          #strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )
      p1
      all_p_list[[sprintf("%s_%d_ratio", g, da)]]<-p1
      n_ext_final_item<-n_ext_final%>%dplyr::filter(group==g)
      
      p2<-ggplot(n_ext_final_item)+
        geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
        geom_tile(aes(x=x, y=y, fill=sum_V))+
        facet_grid(exposure~SSP)+
        scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                            limits=c(0, 40), oob=squish,
                            breaks=c(0, 10, 20, 30, 40),
                            labels=c("0", "10", "20", "30",
                                     sprintf(">40, up to %d", round(max(n_ext_final_item$sum_V)))))+
        #ggtitle(paste(title2, g, sep=" - "))+
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
          panel.border = element_blank(),
          legend.position="bottom",
          legend.key.width=unit(0.8,"in"),
          plot.title = element_text(hjust = 0.5)
        )
      p2
      all_p_list[[sprintf("%s_%d_n_ext", g, da)]]<-p2
      
      #pp<-ggarrange(p2, p1, ncol=1, nrow=2, labels=c("(a)", "(b)"))
      
      #ggsave(pp, 
      #       filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_%d_ttt_%d_%s.png", 
      #                        da, tttt, g), 
      #       width=9, height=8)
      
      #ggsave(pp, 
      #       filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_%d_ttt_%d_%s.pdf", 
      #                        da, tttt, g), 
      #       width=9, height=8)
    #}
  }
}

for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  p_n_ext_no_da_formatted<-all_p_list[[sprintf("%s_%d_n_ext", g, 0)]]
  
  p_ratio_no_da_formatted<-all_p_list[[sprintf("%s_%d_ratio", g, 0)]]+theme(strip.background.x = element_blank(),
                                               strip.text.x = element_blank())
  p_ratio_with_da_formatted<-all_p_list[[sprintf("%s_%d_ratio", g, 1)]]+theme(strip.background.x = element_blank(),
                                                   strip.text.x = element_blank())
  p_n_ext_with_da_formatted<-all_p_list[[sprintf("%s_%d_n_ext", g, 0)]]+theme(strip.background.x = element_blank(),
                                                   strip.text.x = element_blank())
  
  pp<-ggarrange(p_n_ext_no_da_formatted,
                p_ratio_no_da_formatted, 
                p_n_ext_with_da_formatted, 
                p_ratio_with_da_formatted, ncol=1, nrow=4, labels=c("(a)", "(b)", "(c)", "(d)")
                )
  
  ggsave(pp, 
         filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_all_ttt_%d_%s.png", tttt, g), 
         width=7, height=13)
  
  ggsave(pp, 
         filename=sprintf("../../Figures_Full_species/when_where_extinction_all/combined_final_da_all_ttt_%d_%s.pdf", tttt, g), 
         width=7, height=13)
}
