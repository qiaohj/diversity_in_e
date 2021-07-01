library(raster)
#library(rgdal)
library(Rmisc)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")


if (F){
  for (threshold in c(1, 5)){
    GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
    SSPs<-c("SSP119", "SSP245", "SSP585")
    VARs<-c("pr", "tasmax")
    
    predict_range<-c(2021:2100)
    layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
    layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
    
    folders<-c(sprintf("Diversity_%d", threshold))
    #df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
    i=1
    j=4
    k=1
    #dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
    dispersals<-c(0:1)
    sp_dis_all<-NULL
    folder<-folders[1]
    group<-"Amphibians"
    for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      for (j in c(1:nrow(layer_df))){
        layer<-layer_df[j,]
        for (k in c(1:length(dispersals))){
          for (folder in folders){
            layer$M<-dispersals[k]
            layer$TYPE<-folder
            target_folder<-sprintf("../../Objects_Full_species/%s/%s/%s_%d", folder, group, layer$LABEL, layer$M)
            target<-sprintf("%s/indices_df.rda", target_folder)
            
            print(paste("READING DATA", target_folder))
            
            indices_df<-readRDS(target)
            sp_dis<-readRDS(sprintf("%s/sp_dis.rda", target_folder))
            colnames(sp_dis)[which(colnames(sp_dis)=="N")]<-"N_CELL"
            
            keys<-seq(2021, 2100, by=1)
            sp_dis_key<-NULL
            start_dis<-sp_dis%>%filter(year==2020)
            colnames(start_dis)[c(1, 3:ncol(start_dis))]<-paste("st", colnames(start_dis)[c(1, 3:ncol(start_dis))], sep="_")
            colnames(start_dis)[which(colnames(start_dis)=="st_N")]<-"st_N_CELL"
            for (key in keys){
              sub<-sp_dis%>%filter(year==key)
              sub<-full_join(sub, start_dis, by=c("sp"))
              sub$year<-key
              sp_dis_key<-bind(sp_dis_key, sub)
            }
            
            sp_dis_key[which(is.na(sp_dis_key$N_CELL)), "N_CELL"]<-0
            sp_dis_key$group<-group
            sp_dis_key$GCM<-layer$GCM
            sp_dis_key$SSP<-layer$SSP
            sp_dis_key$M<-layer$M
            sp_dis_key$TYPE<-folder
            sp_dis_key$N_type<-""
            sp_dis_key[which(sp_dis_key$N_CELL>sp_dis_key$st_N_CELL), "N_type"]<-"INCREASE"
            sp_dis_key[which(sp_dis_key$N_CELL<sp_dis_key$st_N_CELL), "N_type"]<-"DECREASE"
            sp_dis_key[which(sp_dis_key$N_CELL==sp_dis_key$st_N_CELL), "N_type"]<-"STABLE"
            sp_dis_key[which(sp_dis_key$N_CELL==0), "N_type"]<-"EXTINCT"
            sp_dis_all<-bind(sp_dis_all, sp_dis_key)
          }
        }
      }
    }
    
    saveRDS(sp_dis_all, sprintf("../../Objects_Full_species/N_Cell_year/sp_dis_all_%d.rda", threshold))
    
  }
}
ttt<-1
if (F){
  all_df<-NULL
  for (threshold in c(1, 5)){
    rda<-sprintf("../../Objects_Full_species/N_Cell_year/sp_dis_all_%d.rda", threshold)
    print(paste("Reading", rda))
    sp_dis_all<-readRDS(rda)
    sp_dis_all$N_CELL_Change<-sp_dis_all$N_CELL/sp_dis_all$st_N_CELL
    max_N<-sp_dis_all%>%dplyr::group_by(sp, group, GCM, SSP, M, TYPE)%>%
      dplyr::summarise(MAX_N_CELL=max(N_CELL, st_N_CELL))
    sp_dis_all<-left_join(sp_dis_all, max_N, by=c("sp", "group", "GCM", "SSP", "M", "TYPE"))
    sp_dis_all$N_CELL_Change_Max<-sp_dis_all$N_CELL/sp_dis_all$MAX_N_CELL
    
    sp_dis_all$Label<-paste(sp_dis_all$sp, sp_dis_all$group, 
                            sp_dis_all$GCM, sp_dis_all$SSP, 
                            sp_dis_all$M, sp_dis_all$TYPE)
    for (ttt in c(2)){
      sp_dis_all_filtered<-sp_dis_all%>%dplyr::filter(st_N_CELL>ttt)
      sp_dis_all_sub_1<-sp_dis_all_filtered%>%dplyr::filter(year==2100)
      extincted_species<-sp_dis_all_sub_1%>%dplyr::filter(N_type=="EXTINCT")
      sp_dis_all_extincted<-sp_dis_all_filtered%>%dplyr::filter(Label %in% extincted_species$Label)
      
      
      
      sp_dis_all_extincted_se<-sp_dis_all_extincted%>%
        dplyr::group_by(year, group, SSP, M, TYPE)%>%
        dplyr::summarise(mean_N_CELL=mean(N_CELL),
                         sd_N_CELL=sd(N_CELL),
                         CI_N_CELL=CI(N_CELL)[1]-CI(N_CELL)[2],
                         mean_N_CELL_Change=mean(N_CELL_Change),
                         sd_N_CELL_Change=sd(N_CELL_Change),
                         CI_N_CELL_Change=CI(N_CELL_Change)[1]-CI(N_CELL_Change)[2],
                         mean_N_CELL_Change_Max=mean(N_CELL_Change_Max),
                         sd_N_CELL_Change_Max=sd(N_CELL_Change_Max),
                         CI_N_CELL_Change_Max=CI(N_CELL_Change_Max)[1]- CI(N_CELL_Change_Max)[2])
      
      sp_dis_all_extincted_se_all<-sp_dis_all_extincted%>%
        dplyr::group_by(year, SSP, M, TYPE)%>%
        dplyr::summarise(mean_N_CELL=mean(N_CELL),
                          sd_N_CELL=sd(N_CELL),
                          CI_N_CELL=CI(N_CELL)[1]-CI(N_CELL)[2],
                          mean_N_CELL_Change=mean(N_CELL_Change),
                          sd_N_CELL_Change=sd(N_CELL_Change),
                          CI_N_CELL_Change=CI(N_CELL_Change)[1]-CI(N_CELL_Change)[2],
                          mean_N_CELL_Change_Max=mean(N_CELL_Change_Max),
                          sd_N_CELL_Change_Max=sd(N_CELL_Change_Max),
                          CI_N_CELL_Change_Max=CI(N_CELL_Change_Max)[1]- CI(N_CELL_Change_Max)[2],
                         group="ALL")
      sp_dis_all_extincted_se_all<-sp_dis_all_extincted_se_all%>%dplyr::select(names(sp_dis_all_extincted_se))
      sp_dis_all_extincted_se<-bind_rows(sp_dis_all_extincted_se, sp_dis_all_extincted_se_all)
      if (threshold==1){
        sp_dis_all_extincted_se$exposure<-" no exposure"
      }else{
        sp_dis_all_extincted_se$exposure<-"5-year exposure"
      }
      sp_dis_all_extincted_se$threshold<-ttt
      if (is.null(all_df)){
        all_df<-sp_dis_all_extincted_se
      }else{
        all_df<-bind_rows(all_df, sp_dis_all_extincted_se)
      }
    }
  }
  saveRDS(all_df, "../../Objects_Full_species/N_Cell_year/N_Cell_year.rda")
}
all_df<-readRDS("../../Objects_Full_species/N_Cell_year/N_Cell_year.rda")
all_df$dispersal<-ifelse(all_df$M==0, "no dispersal", "with dispersal")
ttt<-2
g<-"ALL"
for (g in c("ALL", "Amphibians", "Birds", "Mammals", "Reptiles")){
  for (ttt in c(2)){
    p<-ggplot(all_df%>%dplyr::filter((group==g)&(threshold==ttt)))+
      geom_ribbon(aes(x=year, ymin=mean_N_CELL_Change_Max-sd_N_CELL_Change_Max,
                      ymax=mean_N_CELL_Change_Max+sd_N_CELL_Change_Max, fill=SSP), alpha=0.2)+
      geom_line(aes(x=year, y=mean_N_CELL_Change_Max, color=SSP))+
      theme_bw()+
      #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
      facet_grid(exposure~dispersal, scale="free")+
      scale_color_manual(values=color_ssp)+
      scale_fill_manual(values=color_ssp)+
      labs(color = "SSP")+
      xlab("Year")+ylab("Number of cells/Max number of cells")
      #ggtitle(sprintf("Distribution>%d (%s)", ttt, g))
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_MAX_CELL_%s_%d.pdf", g, ttt), width=8, height=4)
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_MAX_CELL_%s_%d.png", g, ttt), width=8, height=4)
    p<-ggplot(all_df%>%dplyr::filter((group==g)&(threshold==ttt)))+
      geom_ribbon(aes(x=year, ymin=mean_N_CELL-sd_N_CELL,
                      ymax=mean_N_CELL+sd_N_CELL, fill=SSP), alpha=0.2)+
      geom_line(aes(x=year, y=mean_N_CELL, color=SSP))+
      theme_bw()+
      #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
      facet_grid(exposure~dispersal, scale="free")+
      scale_color_manual(values=color_ssp)+
      scale_fill_manual(values=color_ssp)+
      labs(color = "SSP")+
      xlab("Year")+ylab("Number of cells")
      #ggtitle(sprintf("Distribution>%d (%s)", ttt, g))
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_%s_%d.pdf", g, ttt), width=8, height=4)
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_%s_%d.png", g, ttt), width=8, height=4)
    p<-ggplot(all_df%>%dplyr::filter((group==g)&(threshold==ttt)))+
      #geom_ribbon(aes(x=year, ymin=mean_N_CELL_Change-N_CELL_Change,
      #                ymax=mean_N_CELL_Change+sd_N_CELL_Change, fill=SSP), alpha=0.2)+
      geom_line(aes(x=year, y=mean_N_CELL_Change, color=SSP))+
      theme_bw()+
      #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
      facet_grid(exposure~dispersal, scale="free")+
      scale_color_manual(values=color_ssp)+
      scale_fill_manual(values=color_ssp)+
      labs(color = "SSP")+
      xlab("Year")+ylab("Number of cells/Initial number of cells")
      #ggtitle(sprintf("Distribution>%d (%s)", ttt, g))
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_INIT_CELL_%s_%d.pdf", g, ttt))
    ggsave(p, filename=sprintf("../../Figures_Full_species/N_Cell_year/N_CELL_INIT_CELL_%s_%d.png", g, ttt))
  }
}
