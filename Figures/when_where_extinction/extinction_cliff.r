library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
threshold<-1

if (F){
  result<-NULL
  for (threshold in c(1,5)){
    g<-"Mammals"
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
      df<-df%>%dplyr::filter(!is.infinite(extinct_year))
      
      sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
      sp_list$sp2<-gsub(" ", "_", sp_list$sp)
      df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
      extinct_sp_list<-df%>%dplyr::distinct(group, sp, GCM, SSP, dispersal, extinct_year)
      i=10
      for (i in c(1:nrow(extinct_sp_list))){
        print(paste(threshold, g, i, nrow(extinct_sp_list)))
        item<-extinct_sp_list[i,]
        df_item<-df%>%dplyr::filter((sp==item$sp)&(GCM==item$GCM)&
                                      (SSP==item$SSP)&(dispersal==item$dispersal))
        N_st_Cell<-nrow(df_item)
        if (item$extinct_year==2021){
          N_extinct_Cell<-N_st_Cell
          Max_N_Cell<-N_st_Cell
        }else{
          future_dis<-readRDS(sprintf("../../Objects_Full_species/Niche_Models/%s/%s/dispersal_%d/%s_%s_%d.rda", 
                                      item$group, item$sp, threshold, item$GCM, item$SSP, item$dispersal))
          future_dis_sum<-future_dis%>%dplyr::group_by(YEAR)%>%
            dplyr::summarise(N_CELL=n())
          Max_N_Cell<-max(future_dis_sum$N_CELL)
          N_extinct_Cell<-pull(future_dis_sum[which(future_dis_sum$YEAR==(item$extinct_year-1)), "N_CELL"])
        }
        item$N_st_Cell<-N_st_Cell
        item$N_extinct_Cell<-N_extinct_Cell
        item$Max_N_Cell<-Max_N_Cell
        item$exposure<-ifelse(threshold==1, " no exposure", "5-year exposure")
        result<-bind_dplyr(result, item)
      }
    }
  }
  saveRDS(result, "../../Figures_Full_species/Extinction_cliff/Extinction_cliff.rda")
}

cliff<-readRDS("../../Figures_Full_species/Extinction_cliff/Extinction_cliff.rda")



range(cliff$N_extinct_Cell)
hist(cliff$N_extinct_Cell)



for (ttt in c(0, 1, 2)){
  
  cliff_item<-cliff%>%dplyr::filter(N_st_Cell>ttt)
  cliff_se<-cliff_item%>%dplyr::group_by(SSP, dispersal, exposure, group)%>%
    dplyr::summarise(mean_extinct_Cell=mean(N_extinct_Cell),
                     sd_extinct_Cell=sd(N_extinct_Cell),
                     CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])
  
  cliff_se_all<-cliff_item%>%dplyr::group_by(SSP, dispersal, exposure)%>%
    dplyr::summarise(group="ALL",
                     mean_extinct_Cell=mean(N_extinct_Cell),
                     sd_extinct_Cell=sd(N_extinct_Cell),
                     CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])
  cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")
  cliff_item$da=ifelse(cliff_item$dispersal==0, "no dispersal", "with dispersal")
  cliff_item[which(cliff_item$exposure=="no exposure"), "exposure"]<-" no exposure"
  cliff_item[which(cliff_item$exposure==" 5-year exposure"), "exposure"]<-"5-year exposure"
  
  p<-ggplot(cliff_item)+geom_density(aes(x=N_extinct_Cell, color=da))+
    facet_grid(exposure~SSP)+
    xlab("N Cells")+
    ylab("Density")+
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(color="")+
    scale_x_log10()+
    scale_color_manual(values=color_da)+
    theme_bw()
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_ttt_%d.pdf", ttt), 
         width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_ttt_%d.png", ttt), 
         width=12, height=6)
  
  
  cliff_se_all<-cliff_item%>%dplyr::group_by(SSP, dispersal, exposure, extinct_year)%>%
    dplyr::summarise(group="ALL",
                     mean_extinct_Cell=mean(N_extinct_Cell),
                     sd_extinct_Cell=sd(N_extinct_Cell),
                     CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])
  cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")
  
  p<-ggplot(cliff_se_all)+
    #geom_ribbon(aes(x=extinct_year,
    #                ymin=mean_extinct_Cell-CI_extinct_Cell, 
    #                ymax=mean_extinct_Cell+CI_extinct_Cell, 
    #                fill=da), alpha=0.2)+
    geom_line(aes(x=extinct_year, y=mean_extinct_Cell, color=da))+
    facet_grid(exposure~SSP, scale="free")+
    xlab("Year")+
    ylab("N Cells")+
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(color="", fill="")+
    scale_color_manual(values=color_da)+
    scale_fill_manual(values=color_da)+
    theme_bw()
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff__by_year_ttt_%d.pdf", ttt), 
         width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff__by_year_ttt_%d.png", ttt), 
         width=12, height=6)
  
  cliff_item$extinct_proportion_st<-cliff_item$N_extinct_Cell/cliff_item$N_st_Cell
  cliff_item$Max_N_Cell<-ifelse(cliff_item$N_extinct_Cell>cliff_item$Max_N_Cell, cliff_item$N_extinct_Cell, cliff_item$Max_N_Cell)
  cliff_item$extinct_proportion_max<-cliff_item$N_extinct_Cell/cliff_item$Max_N_Cell
  
  cliff_se_all<-cliff_item%>%dplyr::group_by(SSP, dispersal, exposure, extinct_year)%>%
    dplyr::summarise(group="ALL",
                     mean_extinct_st=mean(extinct_proportion_st),
                     sd_extinct_st=sd(extinct_proportion_st),
                     CI_extinct_st=CI(extinct_proportion_st)[1]-CI(extinct_proportion_st)[2],
                     mean_extinct_max=mean(extinct_proportion_max),
                     sd_extinct_max=sd(extinct_proportion_max),
                     CI_extinct_max=CI(extinct_proportion_max)[1]-CI(extinct_proportion_max)[2])
  cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")
  
  p<-ggplot(cliff_se_all)+
    #geom_ribbon(aes(x=extinct_year,
    #                ymin=mean_extinct_Cell-CI_extinct_Cell, 
    #                ymax=mean_extinct_Cell+CI_extinct_Cell, 
    #                fill=da), alpha=0.2)+
    geom_line(aes(x=extinct_year, y=mean_extinct_st, color=da))+
    facet_grid(exposure~SSP, scale="free")+
    xlab("Year")+
    ylab("N Cells/Initial N Cells")+
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(color="", fill="")+
    scale_color_manual(values=color_da)+
    scale_fill_manual(values=color_da)+
    theme_bw()
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_by_year_inital_cell_ttt_%d.pdf", ttt), 
         width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_by_year_inital_cell_ttt_%d.png", ttt), 
         width=12, height=6)
  
  p<-ggplot(cliff_se_all)+
    #geom_ribbon(aes(x=extinct_year,
    #                ymin=mean_extinct_Cell-CI_extinct_Cell, 
    #                ymax=mean_extinct_Cell+CI_extinct_Cell, 
    #                fill=da), alpha=0.2)+
    geom_line(aes(x=extinct_year, y=mean_extinct_max, color=da))+
    facet_grid(exposure~SSP, scale="free")+
    xlab("Year")+
    ylab("N Cells/Max N Cells")+
    ggtitle(sprintf("Distribution>%d", ttt))+
    labs(color="", fill="")+
    scale_color_manual(values=color_da)+
    scale_fill_manual(values=color_da)+
    theme_bw()
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_by_year_max_cell_ttt_%d.pdf", ttt), 
         width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Extinction_cliff/Extinction_cliff_by_year_max_cell_ttt_%d.png", ttt), 
         width=12, height=6)
}
