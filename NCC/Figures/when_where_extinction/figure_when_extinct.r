library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
exposure<-0

if (F){
    for (exposure in c(0,5)){
      g<-"Mammals"
      for (g in c("Birds", "Mammals")){
        df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
        sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
        sp_list$sp2<-gsub(" ", "_", sp_list$SP)
        sp_list$sp<-sp_list$SP
        df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
        when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
        
        when_extinct<-when_extinct%>%dplyr::group_by(group, GCM, SSP, extinct_year, dispersal)%>%
          dplyr::summarise(n_sp=n())
        coms<-when_extinct%>%ungroup()%>%dplyr::distinct(group, GCM, SSP, dispersal)
        i=1
        when_extinct_df<-NULL
        for (i in c(1:nrow(coms))){
          print(paste(i, nrow(coms)))
          com<-coms[i,]
          item<-when_extinct%>%dplyr::filter((group==com$group)&(GCM==com$GCM)&
                                               (SSP==com$SSP)&(dispersal==com$dispersal))
          y=2021
          for (y in c(2021:2100)){
            print(paste(g, i, nrow(coms), y))
            item2<-item%>%dplyr::filter(extinct_year<=y)
            if (nrow(item2)==0){
              com$n_sp<-0
            }else{
              com$n_sp<-sum(item2$n_sp)
            }
            com$extinct_year<-y
            when_extinct_df<-bind_dplyr(when_extinct_df, com)
          }
        }
        saveRDS(when_extinct_df, sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_%s.rda", exposure, g))
      }
    }
  
}
if (F){
  exposure<-0
  for (exposure in c(0,5)){
    when_extinct<-NULL
    
    for (g in c("Birds", "Mammals")){
      sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
      when_extinct_df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_%s.rda", exposure, g))
      when_extinct_df$extinct_ratio<-when_extinct_df$n_sp/nrow(sp_list)
      
      when_extinct_df_se<-when_extinct_df%>%dplyr::group_by(group, SSP, extinct_year, dispersal)%>%
        dplyr::summarise(mean_n_sp=mean(n_sp),
                         sd_n_sp=sd(n_sp),
                         CI_n_sp=CI(n_sp)[2]-CI(n_sp)[3],
                         mean_extinct_ratio=mean(extinct_ratio),
                         sd_extinct_ratio=sd(extinct_ratio),
                         CI_extinct_ratio=CI(extinct_ratio)[2]-CI(extinct_ratio)[3])
      when_extinct_df_se$all_sp<-nrow(sp_list)
      when_extinct<-bind_dplyr(when_extinct, when_extinct_df_se)
    }
    saveRDS(when_extinct, sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", exposure))
  }
  
}

if (F){
  
  for (exposure in c(0,5)){
    when_extinct<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", exposure))
    names(when_extinct)[1]<-"Group"
    when_extinct$dispersal_label<-ifelse(when_extinct$dispersal==0, "no dispersal", "with dispersal")
    p1<-ggplot(when_extinct, aes(x=extinct_year, y=mean_n_sp, color=Group))+
      #geom_errorbar(aes(ymin=mean_n_sp-CI_n_sp, ymax=mean_n_sp+CI_n_sp, color=Group), alpha=0.7, width=0.25)+
      geom_line(aes(linetype=SSP))+
      scale_color_manual(values=color_groups)+
      scale_linetype_manual(values=linetype_ssp)+
      xlab("Year")+
      ylab("Average number of extinctions")+
      theme_bw()+
      facet_wrap(~dispersal_label)
    p1
    
    p2<-ggplot(when_extinct, aes(x=extinct_year, y=mean_extinct_ratio, color=Group))+
      geom_line(aes(linetype=SSP))+
      scale_color_manual(values=color_groups)+
      scale_linetype_manual(values=linetype_ssp)+
      xlab("Year")+
      ylab("Extinction proportion")+
      theme_bw()+
      facet_wrap(~dispersal_label)
    p2
    legend_g<-get_legend(p1)
    p<-ggarrange(p1, p2, ncol=1, common.legend = T, legend="right", legend.grob=legend_g)
    p
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when.pdf", exposure), width=8, height=6)
    ggsave(p, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when.png", exposure), width=8, height=6)
    
    ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when_number.pdf", exposure), width=8, height=4)
    ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when_number.png", exposure), width=8, height=4)
    
    ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when_proportion.pdf", exposure), width=8, height=4)
    ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_exposure_%d/when_proportion.png", exposure), width=8, height=4)
  }
  
}


  when_extinct_1<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", 0))
  names(when_extinct_1)[1]<-"Group"
  when_extinct_1$label<-ifelse((when_extinct_1$dispersal==0), 
                               "no dispersal, no exposure",
                               "with dispersal, no exposure")
  when_extinct_1$exposure<-" no exposure"
  when_extinct_1$da<-ifelse((when_extinct_1$dispersal==0), 
                            "no dispersal",
                            "with dispersal")
  
  when_extinct_5<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", 5))
  names(when_extinct_5)[1]<-"Group"
  when_extinct_5$label<-ifelse((when_extinct_5$dispersal==0), 
                               "no dispersal, 5-year exposure",
                               "with dispersal, 5-year exposure")
  when_extinct_5$exposure<-"5-year exposure"
  when_extinct_5$da<-ifelse((when_extinct_5$dispersal==0), 
                            "no dispersal",
                            "with dispersal")
  
  
  
  when_extinct<-bind_rows(when_extinct_1, when_extinct_5)
  
  df_item<-when_extinct
  p1<-ggplot()+
    geom_ribbon(data=df_item, 
                aes(x=extinct_year,
                    ymin=mean_n_sp-CI_n_sp, 
                    ymax=mean_n_sp+CI_n_sp, 
                    fill=Group, linetype=SSP), alpha=0.2)+
    geom_line(data=df_item, aes(x=extinct_year, y=mean_n_sp, 
                                color=Group, linetype=SSP))+
    scale_color_manual(values=color_groups)+
    scale_fill_manual(values=color_groups)+
    scale_linetype_manual(values=linetype_ssp)+
    xlab("Year")+
    ylab("Average number of extinctions")+
    labs(linetype = "SSP scenario")+
    theme_bw()+
    #xlim(2020, 2110)+
    facet_grid(exposure~da, scale="free")
  
  p1
  ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_all/when_number.pdf"),
         width=11, height=6)
  ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_all/when_number.png"),
         width=11, height=6)
  
  df_item<-when_extinct
  p2<-ggplot(when_extinct)+
    geom_ribbon(data=df_item, 
                aes(x=extinct_year,
                    ymin=mean_extinct_ratio-CI_extinct_ratio, 
                    ymax=mean_extinct_ratio+CI_extinct_ratio, 
                    fill=Group, linetype=SSP), alpha=0.2)+
    
    geom_line(data=df_item, 
              aes(x=extinct_year, y=mean_extinct_ratio, 
                  color=Group, linetype=SSP))+
    scale_color_manual(values=color_groups)+
    scale_fill_manual(values=color_groups)+
    scale_linetype_manual(values=linetype_ssp)+
    xlab("Year")+
    ylab("Extinction proportion")+
    labs(linetype = "SSP scenario")+
    theme_bw()+
    #xlim(2020, 2110)+
    facet_grid(exposure~da)
  p2
  ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_all/when_proportion.pdf"),
         width=11, height=6)
  ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_all/when_proportion.png"),
         width=11, height=6)
  
  
  write.csv(when_extinct,
            sprintf("../../Figures/when_where_extinction_all/when_extinct.csv"))
