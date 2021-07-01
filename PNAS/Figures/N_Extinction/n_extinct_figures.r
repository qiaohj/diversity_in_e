library(raster)
#library(rgdal)
#library(rgeos)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
exposure=0

if (F){
  for (exposure in c(0, 5)){
    GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
    SSPs<-c("SSP119", "SSP245", "SSP585")
    
    predict_range<-c(2021:2100)
    layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
    layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
    
    folders<-c(sprintf("Diversity_exposure_%d", exposure))
    #df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
    i=1
    j=4
    k=1
    #dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
    sp_dis_all<-NULL
    sp_dis_all_1<-NULL
    sp_dis_all_2<-NULL
    folder<-folders[1]
    group<-"Mammals"
    for (group in c("Birds", "Mammals")){
      for (j in c(1:nrow(layer_df))){
        layer<-layer_df[j,]
        
        for (folder in folders){
          layer$TYPE<-folder
          target_folder<-sprintf("../../Objects/%s/%s/%s", folder, group, layer$LABEL)
          target<-sprintf("%s/indices_df.rda", target_folder)
          
          print(paste("READING DATA", target_folder))
          
          indices_df<-readRDS(target)
          sp_dis<-readRDS(sprintf("%s/sp_dis.rda", target_folder))
          colnames(sp_dis)[which(colnames(sp_dis)=="N")]<-"N_CELL"
          
          keys<-seq(2040, 2100, by=20)
          sp_dis_key<-NULL
          start_dis<-sp_dis%>%filter(YEAR==2020)
          colnames(start_dis)[c(1, 3:ncol(start_dis))]<-
            paste("st", colnames(start_dis)[c(1, 3:ncol(start_dis))], sep="_")
          colnames(start_dis)[which(colnames(start_dis)=="st_N")]<-"st_N_CELL"
          for (key in keys){
            sub<-sp_dis%>%filter(YEAR==key)
            sub<-full_join(sub, start_dis, by=c("sp"))
            sub$YEAR<-key
            sp_dis_key<-bind(sp_dis_key, sub)
          }
          
          sp_dis_key[which(is.na(sp_dis_key$N_CELL)), "N_CELL"]<-0
          sp_dis_key$group<-group
          sp_dis_key$GCM<-layer$GCM
          sp_dis_key$SSP<-layer$SSP
          
          sp_dis_key$TYPE<-folder
          sp_dis_key$N_type<-""
          sp_dis_key[which(sp_dis_key$N_CELL>sp_dis_key$st_N_CELL), "N_type"]<-"INCREASE"
          sp_dis_key[which(sp_dis_key$N_CELL<sp_dis_key$st_N_CELL), "N_type"]<-"DECREASE"
          sp_dis_key[which(sp_dis_key$N_CELL==sp_dis_key$st_N_CELL), "N_type"]<-"STABLE"
          sp_dis_key[which(sp_dis_key$N_CELL==0), "N_type"]<-"EXTINCT"
          sp_dis_all<-bind(sp_dis_all, sp_dis_key)
          
          #ggplot(sp_dis_key, aes(x=YEAR, fill=factor(N_type)))+geom_bar()
        }
        
      }
    }
    N_SP<-sp_dis_all%>%dplyr::group_by(group)%>%dplyr::summarise(N_SP=n_distinct(sp))
    
    sp_dis_all<-inner_join(sp_dis_all, N_SP, by=c("group"))
    sp_dis_all$Label1<-paste(sp_dis_all$GCM, sp_dis_all$SSP)
    saveRDS(sp_dis_all, sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", exposure))
    
  }
}
ttt<-0
exposure<-0
for (ttt in c(0)){
  sp_dis_all_sub_N_all<-NULL
  sp_dis_extinct<-NULL
  sp_dis_all_se_all<-NULL
  for (exposure in c(0, 5)){
    rda<-sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", exposure)
    print(paste("Reading", rda))
    sp_dis_all<-readRDS(rda)
    sp_dis_all_sub_1<-sp_dis_all%>%dplyr::filter(YEAR==2100)
    sp_dis_all_sub<-sp_dis_all_sub_1%>%dplyr::filter(N_type=="EXTINCT")
    sp_dis_all_sub_N<-sp_dis_all_sub%>%dplyr::group_by(group, Label1, GCM, SSP, N_type, N_SP, TYPE)%>%
      dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
    sp_dis_all_sub_N$persentile<-sp_dis_all_sub_N$N_SP_EXTINCT/sp_dis_all_sub_N$N_SP
    
    sp_dis_all_se<-sp_dis_all_sub%>%dplyr::group_by(group, SSP, N_SP, TYPE)%>%
      dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
    sp_dis_all_se$persentile<-sp_dis_all_se$N_SP_EXTINCT/sp_dis_all_se$N_SP
    
    if (exposure==0){
      sp_dis_all_sub_N$Label<-paste(sp_dis_all_sub_N$group, "  (no exposure)", sep="")
      sp_dis_all_sub_N$exposure<-" no exposure"
      sp_dis_all_sub_1$Label<-paste(sp_dis_all_sub_1$group, "  (no exposure)", sep="")
      sp_dis_all_sub_1$exposure<-" no exposure"
      sp_dis_all_se$Label<-paste(sp_dis_all_se$group, "  (no exposure)", sep="")
      sp_dis_all_se$exposure<-" no exposure"
      
    }else{
      sp_dis_all_sub_N$Label<-paste(sp_dis_all_sub_N$group, " (5-year exposure)", sep="")
      sp_dis_all_sub_N$exposure<-"5-year exposure"
      sp_dis_all_sub_1$Label<-paste(sp_dis_all_sub_1$group, "  (5-year exposure)", sep="")
      sp_dis_all_sub_1$exposure<-"5-year exposure"
      sp_dis_all_se$Label<-paste(sp_dis_all_se$group, "  (5-year exposure)", sep="")
      sp_dis_all_se$exposure<-"5-year exposure"
    }
    sp_dis_extinct<-bind_dplyr(sp_dis_extinct, sp_dis_all_sub_1)
    sp_dis_all_sub_N_all<-bind_dplyr(sp_dis_all_sub_N_all, sp_dis_all_sub_N)
    sp_dis_all_se_all<-bind_dplyr(sp_dis_all_se_all, sp_dis_all_se)
  }
  
  sp_mean_gcm<-sp_dis_all_sub_N_all%>%dplyr::filter(M!=2)%>%
    dplyr::group_by(GCM, SSP, M, N_type)%>%
    dplyr::summarise(persentile_MEAN=sum(persentile*N_SP)/sum(N_SP))
  
  sp_mean4<-sp_mean_gcm%>%
    dplyr::group_by(SSP, M, N_type)%>%
    dplyr::summarise(persentile=mean(persentile_MEAN),
                     persentile_SD=sd(persentile_MEAN))
  
  #sp_mean<-sp_dis_all_sub_N_all%>%dplyr::filter(M!=2)%>%
  #  dplyr::group_by(SSP, M, N_type)%>%
  #  dplyr::summarise(persentile_MEAN=mean(persentile),
  #                   persentile_SD=sd(persentile))
  
  write.csv(sp_mean4, sprintf("../../Figures/N_Extinction/Extinction_%d_all_exposure.csv", ttt))
  
  sp_mean<-sp_dis_all_sub_N_all%>%dplyr::filter(M!=2)%>%
    dplyr::group_by(group, SSP, M, N_type, N_SP, TYPE, Label, exposure)%>%
    dplyr::summarise(persentile_MEAN=mean(persentile),
                     persentile_SD=sd(persentile))
  
  
  sp_mean$exposure<-gsub("\\(", "", sp_mean$exposure)
  sp_mean$exposure<-gsub("\\)", "", sp_mean$exposure)
  write.csv(sp_mean, sprintf("../../Figures/N_Extinction/Extinction_%d.csv", ttt))
  
  sp_mean_gcm<-sp_dis_all_sub_N_all%>%dplyr::filter(M!=2)%>%
    dplyr::group_by(GCM, SSP, M, N_type, exposure)%>%
    dplyr::summarise(persentile_MEAN=sum(persentile*N_SP)/sum(N_SP))
  
  sp_mean_3<-sp_mean_gcm%>%
    dplyr::group_by(SSP, M, N_type, exposure)%>%
    dplyr::summarise(persentile=mean(persentile_MEAN),
                     persentile_SD=sd(persentile_MEAN))
  
  
  write.csv(sp_mean_3, sprintf("../../Figures/N_Extinction/Extinction_%d_by_exposure.csv", ttt))
  
  p<-ggplot(sp_mean, aes(y=persentile_MEAN, x=SSP))+
    geom_bar(stat="identity", position=position_dodge(), aes(fill=factor(M)))+
    geom_errorbar(position=position_dodge(.9), width=0.2,
                  aes(ymin=persentile_MEAN-persentile_SD, 
                      ymax=persentile_MEAN+persentile_SD,
                      group=factor(M))) +
    xlab("SSP scenario")+
    #ggtitle(sprintf("Distribution>%d", ttt))+
    #ylim(c(0, 1))+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
    facet_grid(exposure~group)+
    scale_fill_manual(values=color_dispersal, breaks=c(0:1), 
                      labels = c("no dispersal", "with dispersal"))+
    labs(fill = "Dispersal")+
    ylab("Mean extinction proportion")
  p
  
  
  
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_%d.pdf", ttt), width=10, height=6)
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_%d.png", ttt), width=10, height=6)
  
  sp_dis_extinct<-sp_dis_extinct%>%dplyr::filter(M!=2)
  sp_dis_extinct<-data.frame(sp_dis_extinct)
  sp_dis_extinct[which(sp_dis_extinct$M==0), "Label"]<-
    paste(sp_dis_extinct[which(sp_dis_extinct$M==0), "SSP"], "(no dispersal)")
  sp_dis_extinct[which(sp_dis_extinct$M==1), "Label"]<-
    paste(sp_dis_extinct[which(sp_dis_extinct$M==0), "SSP"], "(with dispersal)")
  
  p<-ggplot(sp_dis_extinct)+
    geom_histogram(aes(x=st_N_CELL), fill=colors_black[4], bins=20)+
    geom_histogram(data=sp_dis_extinct%>%dplyr::filter(N_type=="EXTINCT"), 
                   aes(x=st_N_CELL), fill=colors_red[9], bins=20)+
    scale_x_log10()+
    theme_bw()+
    xlab("Range size (log transformed)")+
    ylab("Number of species")+
    #ggtitle(sprintf("Distribution>%d", ttt))+
    facet_grid(exposure~Label)
  saveRDS(p, "../../Figures/NB_hist_combined/nb_range_size.rda")
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%d.pdf", ttt), width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%d.png", ttt), width=12, height=6)
  
  
  
  for (g in c("Birds", "Mammals")){
    print(g)
    sp_dis_extinct_item<-sp_dis_extinct%>%dplyr::filter(group==g)
    p<-ggplot(sp_dis_extinct_item)+
      geom_histogram(aes(x=st_N_CELL), fill=colors_black[4], bins=20)+
      geom_histogram(data=sp_dis_extinct_item%>%dplyr::filter(N_type=="EXTINCT"), 
                     aes(x=st_N_CELL), fill=colors_red[9], bins=20)+
      #ggtitle(sprintf("Distribution>%d (%s)", ttt, g))+
      scale_x_log10()+
      theme_bw()+
      xlab("Range size")+
      ylab("Number of species")+
      facet_grid(exposure~Label)
    ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%s_%d.pdf", g, ttt), 
           width=12, height=6)
    ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%s_%d.png", g, ttt),
           width=12, height=6)
    
  }
}
