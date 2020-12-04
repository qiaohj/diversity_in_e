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

threshold<-5

if (F){
  for (threshold in c(1, 5)){
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
      sp_continent<-readRDS(sprintf("../../Objects/SP_Continent/%s.rda", group))
      sp_continent<-sp_continent%>%dplyr::filter(!is.na(continent))
      for (j in c(1:nrow(layer_df))){
        layer<-layer_df[j,]
        for (k in c(1:length(dispersals))){
          for (folder in folders){
            layer$M<-dispersals[k]
            layer$TYPE<-folder
            target_folder<-sprintf("../../Objects/%s/%s/%s_%d", folder, group, layer$LABEL, layer$M)
            target<-sprintf("%s/indices_df.rda", target_folder)
            
            print(paste("READING DATA", target_folder))
            
            #indices_df<-readRDS(target)
            sp_dis<-readRDS(sprintf("%s/sp_dis.rda", target_folder))
            colnames(sp_dis)[which(colnames(sp_dis)=="N")]<-"N_CELL"
            sp_dis<-full_join(sp_continent, sp_dis, by="sp")
            #View(sp_dis[which(sp_dis$sp=="Sachatamia_ilex"),])
            
            
            
            keys<-seq(2040, 2100, by=20)
            sp_dis_key<-NULL
            start_dis<-sp_dis%>%filter(year==2020)
            colnames(start_dis)[c(4, 5:ncol(start_dis))]<-
              paste("st", colnames(start_dis)[c(4, 5:ncol(start_dis))], sep="_")
            
            for (key in keys){
              sub<-sp_dis%>%filter(year==key)
              sub<-full_join(sub, start_dis, by=c("sp", "group", "continent"))
              sub$year<-key
              sp_dis_key<-bind_dplyr(sp_dis_key, sub)
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
            sp_dis_all<-bind_dplyr(sp_dis_all, sp_dis_key)
            #ggplot(sp_dis_key, aes(x=year, fill=factor(N_type)))+geom_bar()
          }
        }
      }
    }
    N_SP<-sp_dis_all%>%dplyr::group_by(group, continent)%>%dplyr::summarise(N_SP=n_distinct(sp))
    sp_dis_all<-inner_join(sp_dis_all, N_SP, by=c("group", "continent"))
    sp_dis_all$Label1<-paste(sp_dis_all$GCM, sp_dis_all$SSP)
    sp_dis_all<-sp_dis_all%>%dplyr::filter(!is.na(continent))
    saveRDS(sp_dis_all, sprintf("../../Figures/N_Extinction/sp_dis_continent_%d.rda", threshold))
  }
}

sp_dis_all_sub_N_all<-NULL
sp_dis_extinct<-NULL
for (threshold in c(1, 5)){
  rda<-sprintf("../../Figures/N_Extinction/sp_dis_continent_%d.rda", threshold)
  print(paste("Reading", rda))
  sp_dis_all<-readRDS(rda)
  sp_dis_all_sub_1<-sp_dis_all%>%dplyr::filter(year==2100)
  sp_dis_all_sub<-sp_dis_all_sub_1%>%dplyr::filter(N_type=="EXTINCT")
  sp_dis_all_sub_N<-sp_dis_all_sub%>%
    dplyr::group_by(group, Label1, GCM, SSP, M, N_type, N_SP, TYPE, continent)%>%
    dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
  sp_dis_all_sub_N$persentile<-sp_dis_all_sub_N$N_SP_EXTINCT/sp_dis_all_sub_N$N_SP
  if (threshold==1){
    sp_dis_all_sub_N$Label<-paste(sp_dis_all_sub_N$group, "  no exposure", sep="")
    sp_dis_all_sub_N$exposure<-" no exposure"
    sp_dis_all_sub_1$Label<-paste(sp_dis_all_sub_1$group, "  (no exposure)", sep="")
    sp_dis_all_sub_1$exposure<-" no exposure"
  }else{
    sp_dis_all_sub_N$Label<-paste(sp_dis_all_sub_N$group, " (5-year exposure)", sep="")
    sp_dis_all_sub_N$exposure<-"5-year exposure"
    sp_dis_all_sub_1$Label<-paste(sp_dis_all_sub_1$group, "  (5-year exposure)", sep="")
    sp_dis_all_sub_1$exposure<-"5-year exposure"
  }
  sp_dis_extinct<-bind_dplyr(sp_dis_extinct, sp_dis_all_sub_1)
  sp_dis_all_sub_N_all<-bind_dplyr(sp_dis_all_sub_N_all, sp_dis_all_sub_N)
  
}
sp_dis_extinct<-data.frame(sp_dis_extinct)
sp_dis_extinct<-sp_dis_extinct[!is.na(sp_dis_extinct$continent),]
sp_dis_extinct$continent_label<-NA
#euroasia, Africa, North America, South America, Antarctica, Europe, and Australia

sp_dis_extinct[which(sp_dis_extinct$continent==1), "continent_label"]<-"Africa"
sp_dis_extinct[which(sp_dis_extinct$continent==2), "continent_label"]<-"Euroasia"
sp_dis_extinct[which(sp_dis_extinct$continent==3), "continent_label"]<-"Australia"
sp_dis_extinct[which(sp_dis_extinct$continent==5), "continent_label"]<-"North America"
sp_dis_extinct[which(sp_dis_extinct$continent==6), "continent_label"]<-"South America"

sp_dis_all_sub_N_all<-data.frame(sp_dis_all_sub_N_all)
sp_dis_all_sub_N_all<-sp_dis_all_sub_N_all[!is.na(sp_dis_all_sub_N_all$continent),]
sp_dis_all_sub_N_all$continent_label<-NA
#euroasia, Africa, North America, South America, Antarctica, Europe, and Australia

sp_dis_all_sub_N_all[which(sp_dis_all_sub_N_all$continent==1), "continent_label"]<-"Africa"
sp_dis_all_sub_N_all[which(sp_dis_all_sub_N_all$continent==2), "continent_label"]<-"Euroasia"
sp_dis_all_sub_N_all[which(sp_dis_all_sub_N_all$continent==3), "continent_label"]<-"Australia"
sp_dis_all_sub_N_all[which(sp_dis_all_sub_N_all$continent==5), "continent_label"]<-"North America"
sp_dis_all_sub_N_all[which(sp_dis_all_sub_N_all$continent==6), "continent_label"]<-"South America"



sp_mean<-sp_dis_all_sub_N_all%>%dplyr::filter(M!=2)%>%
  dplyr::group_by(group, SSP, M, N_type, N_SP, TYPE, Label, exposure, continent, continent_label)%>%
  dplyr::summarise(persentile_MEAN=mean(persentile),
                   persentile_SD=sd(persentile))

write.csv(sp_mean, "../../Figures/N_Extinction/Extinction_continent.csv")

sp_mean$exposure<-gsub("\\(", "", sp_mean$exposure)
sp_mean$exposure<-gsub("\\)", "", sp_mean$exposure)
g<-"Amphibians"
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  sp_mean_item<-sp_mean%>%dplyr::filter(group==g)
  p<-ggplot(sp_mean_item, aes(y=persentile_MEAN, x=SSP))+
    geom_bar(stat="identity", position=position_dodge(), aes(fill=factor(M)))+
    geom_errorbar(position=position_dodge(.9), width=0.2,
                  aes(ymin=persentile_MEAN-persentile_SD, 
                      ymax=persentile_MEAN+persentile_SD,
                      group=factor(M))) +
    xlab("SSP scenario")+
    ggtitle(g)+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
    facet_grid(exposure~continent_label)+
    scale_fill_manual(values=color_dispersal, breaks=c(0:1), 
                      labels = c("no dispersal", "with dispersal"))+
    labs(fill = "Dispersal")+
    ylab("Mean extinction proportion")
  p
  
  
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_continent_%s.pdf", g),
         width=14, height=6)
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_continent_%s.png", g),
         width=14, height=6)
}

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
  xlab("Range size")+
  ylab("Number of species")+
  facet_grid(exposure~Label)
ggsave(p, filename="../../Figures/N_Extinction/Extinction_hist.pdf", width=12, height=6)
ggsave(p, filename="../../Figures/N_Extinction/Extinction_hist.png", width=12, height=6)


for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  sp_dis_extinct_item<-sp_dis_extinct%>%dplyr::filter(group==g)
  p<-ggplot(sp_dis_extinct_item)+
    geom_histogram(aes(x=st_N_CELL), fill=colors_black[4], bins=20)+
    geom_histogram(data=sp_dis_extinct_item%>%dplyr::filter(N_type=="EXTINCT"), 
                   aes(x=st_N_CELL), fill=colors_red[9], bins=20)+
    ggtitle(g)+
    scale_x_log10()+
    theme_bw()+
    xlab("Range size")+
    ylab("Number of species")+
    facet_grid(exposure~Label)
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%s.pdf", g), 
         width=12, height=6)
  ggsave(p, filename=sprintf("../../Figures/N_Extinction/Extinction_hist_%s.png", g),
         width=12, height=6)
  
}
