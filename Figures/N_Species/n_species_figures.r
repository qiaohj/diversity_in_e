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

threshold<-5
if (F){
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  VARs<-c("pr", "tasmax")
  
  predict_range<-c(2021:2100)
  layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
  layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
  
  folders<-c(sprintf("Diversity_%d", threshold))
  #df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
  i=1
  j=4
  k=1
  #dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
  dispersals<-c(0:2)
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
          target_folder<-sprintf("../../Objects/%s/%s/%s_%d", folder, group, layer$LABEL, layer$M)
          target<-sprintf("%s/indices_df.rda", target_folder)
          
          print(paste("READING DATA", target_folder))
          
          indices_df<-readRDS(target)
          sp_dis<-readRDS(sprintf("%s/sp_dis.rda", target_folder))
          colnames(sp_dis)[which(colnames(sp_dis)=="N")]<-"N_CELL"
          
          keys<-seq(2040, 2100, by=20)
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
          #ggplot(sp_dis_key, aes(x=year, fill=factor(N_type)))+geom_bar()
        }
      }
    }
  }
  N_SP<-sp_dis_all%>%dplyr::group_by(group)%>%dplyr::summarise(N_SP=n_distinct(sp))
  sp_dis_all<-inner_join(sp_dis_all, N_SP, by=c("group"))
  sp_dis_all$Label1<-paste(sp_dis_all$GCM, sp_dis_all$SSP)
  sp_dis_all$Label2<-paste(sp_dis_all$M, sp_dis_all$N)
  saveRDS(sp_dis_all, "../../Figures/N_SPECIES_5/sp_dis_all.rda")
}


colnames(sp_dis_all)
sp_dis_all<-readRDS("../../Figures/N_SPECIES_5/sp_dis_all.rda")
sp_dis_all_sub<-sp_dis_all%>%dplyr::filter(year==2100)
sp_dis_all_sub<-sp_dis_all_sub%>%dplyr::filter(N_type=="EXTINCT")
sp_dis_all_sub_N<-sp_dis_all_sub%>%dplyr::group_by(group, Label1, Label2, GCM, SSP, M, N_type, N_SP, TYPE)%>%
  dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
sp_dis_all_sub_N$persentile<-sp_dis_all_sub_N$N_SP_EXTINCT/sp_dis_all_sub_N$N_SP
sp_dis_all_sub_N_mean<-sp_dis_all_sub_N%>%dplyr::group_by(group, Label2, SSP, M, N_type, N_SP, TYPE)%>%
  dplyr::summarise(persentile_MEAN=mean(persentile),
                   persentile_SD=sd(persentile))

head(sp_dis_all_sub_N)

p<-ggplot(sp_dis_all_sub_N_mean, aes(y=persentile_MEAN, x=SSP))+
  geom_bar(stat="identity", position=position_dodge(), aes(fill=factor(Label2)))+
  geom_errorbar(position=position_dodge(.9), width=0.2,
                aes(ymin=persentile_MEAN-persentile_SD, 
                    ymax=persentile_MEAN+persentile_SD,
                    group=factor(Label2))) +
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  facet_wrap( ~ group+TYPE, ncol=2)
ggsave(p, filename="../../Figures/N_SPECIES_5/Extinction.pdf", width=6, height=6)
ggsave(p, filename="../../Figures/N_SPECIES_5/Extinction.png", width=6, height=6)
