library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
threshold<-1
ttt=0
if (F){
  for (ttt in c(0, 1, 2)){
    for (threshold in c(1,5)){
      g<-"Mammals"
      for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
        df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
        sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
        sp_list<-sp_list[which(sp_list$area>ttt),]
        sp_list$sp2<-gsub(" ", "_", sp_list$sp)
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
            print(paste(g, i, nrow(coms), y, ttt))
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
        saveRDS(when_extinct_df, sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_%s_%d.rda", threshold, g, ttt))
      }
    }
  }
}
if (F){
  threshold<-1
  for (ttt in c(0, 1, 2)){
    for (threshold in c(1,5)){
      when_extinct<-NULL
      
      for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
        sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
        sp_list<-sp_list[which(sp_list$area>ttt),]
        
        when_extinct_df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_%s_%d.rda", threshold, g, ttt))
        when_extinct_df$extinct_ratio<-when_extinct_df$n_sp/nrow(sp_list)
        when_extinct_df%>%dplyr::filter((SSP=="SSP119")&(extinct_year==2021))
             + )
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
      saveRDS(when_extinct, sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_final_%d.rda", threshold, ttt))
    }
  }
}

if (F){
  for (ttt in c(0, 1, 2)){
    for (threshold in c(1,5)){
      when_extinct<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_final_%d.rda", threshold, ttt))
      names(when_extinct)[1]<-"Group"
      p1<-ggplot(when_extinct%>%dplyr::filter(dispersal!=2), aes(x=extinct_year, y=mean_n_sp, color=Group))+
        #geom_errorbar(aes(ymin=mean_n_sp-CI_n_sp, ymax=mean_n_sp+CI_n_sp, color=Group), alpha=0.7, width=0.25)+
        geom_line(aes(linetype=SSP))+
        scale_color_manual(values=color_groups)+
        scale_linetype_manual(values=linetype_ssp)+
        xlab("Year")+
        ylab("Average number of extinctions")+
        ggtitle(sprintf("Distribution>%d", ttt))+
        theme_bw()+
        facet_wrap(~dispersal)
      p1
      
      p2<-ggplot(when_extinct%>%dplyr::filter(dispersal!=2), aes(x=extinct_year, y=mean_extinct_ratio, color=Group))+
        geom_line(aes(linetype=SSP))+
        scale_color_manual(values=color_groups)+
        scale_linetype_manual(values=linetype_ssp)+
        xlab("Year")+
        ylab("Extinction proportion")+
        ggtitle(sprintf("Distribution>%d", ttt))+
        theme_bw()
      p2
      legend_g<-get_legend(p1)
      p<-ggarrange(p1, p2, ncol=1, common.legend = T, legend="right", legend.grob=legend_g)
      p
      ggsave(p, filename=sprintf("../../Figures/when_where_extinction_%d/when_%d.pdf", threshold, ttt), width=8, height=6)
      ggsave(p, filename=sprintf("../../Figures/when_where_extinction_%d/when_%d.png", threshold, ttt), width=8, height=6)
      
      ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_%d/when_number_%d.pdf", threshold, ttt), width=8, height=4)
      ggsave(p1, filename=sprintf("../../Figures/when_where_extinction_%d/when_number_%d.png", threshold, ttt), width=8, height=4)
      
      ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_%d/when_proportion_%d.pdf", threshold, ttt), width=8, height=4)
      ggsave(p2, filename=sprintf("../../Figures/when_where_extinction_%d/when_proportion_%d.png", threshold, ttt), width=8, height=4)
    }
  }
}

for (ttt in c(0, 1, 2)){
  when_extinct_1<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_final_%d.rda", 1, ttt))
  names(when_extinct_1)[1]<-"Group"
  when_extinct_1$label<-ifelse((when_extinct_1$dispersal==0), 
                               "no dispersal, no exposure",
                               "with dispersal, no exposure")
  when_extinct_1$exposure<-" no exposure"
  when_extinct_1$da<-ifelse((when_extinct_1$dispersal==0), 
                            "no dispersal",
                            "with dispersal")
  
  when_extinct_5<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/when_extinct_final_%d.rda", 5, ttt))
  names(when_extinct_5)[1]<-"Group"
  when_extinct_5$label<-ifelse((when_extinct_5$dispersal==0), 
                               "no dispersal, 5-year exposure",
                               "with dispersal, 5-year exposure")
  when_extinct_5$exposure<-"5-year exposure"
  when_extinct_5$da<-ifelse((when_extinct_5$dispersal==0), 
                            "no dispersal",
                            "with dispersal")
  
  
  
  when_extinct<-bind_rows(when_extinct_1, when_extinct_5)
  
  df_item<-when_extinct%>%dplyr::filter(dispersal!=2)
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
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(linetype = "SSP scenario")+
    theme_bw()+
    #xlim(2020, 2110)+
    facet_grid(exposure~da, scale="free")
  
  p1
  ggsave(p1, filename=sprintf("../../Figures_Full_species/when_where_extinction_all/when_number_%d.pdf", ttt),
         width=11, height=6)
  ggsave(p1, filename=sprintf("../../Figures_Full_species/when_where_extinction_all/when_number_%d.png", ttt),
         width=11, height=6)
  
  df_item<-when_extinct%>%dplyr::filter(dispersal!=2)
  p2<-ggplot(when_extinct%>%dplyr::filter(dispersal!=2),
  )+
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
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(linetype = "SSP scenario")+
    theme_bw()+
    #xlim(2020, 2110)+
    facet_grid(exposure~da)
  p2
  ggsave(p2, filename=sprintf("../../Figures_Full_species/when_where_extinction_all/when_proportion_%d.pdf", ttt),
         width=11, height=6)
  ggsave(p2, filename=sprintf("../../Figures_Full_species/when_where_extinction_all/when_proportion_%d.png", ttt),
         width=11, height=6)
  
  
  write.csv(when_extinct%>%dplyr::filter((dispersal!=2)),
            sprintf("../../Figures_Full_species/when_where_extinction_all/when_extinct_%d.csv", ttt))
}