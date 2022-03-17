library(raster)
#library(rgdal)
#library(rgeos)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pROC)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

sp_dis_10km_all_sub_N_all<-NULL
sp_dis_10km_extinct<-NULL
sp_dis_10km_all_se_all<-NULL

for (exposure in c(0, 5)){
  rda<-sprintf("../../Figures/N_Extinction/sp_dis_all_%d_10km.rda", exposure)
  print(paste("Reading", rda))
  sp_dis_10km_all<-readRDS(rda)
  
  if (F){
    test<-sp_dis_10km_all[, .(N=.N), by=list(GCM, SSP, N_type, N_SP, M, group, year)]
    test<-test[N_type=="EXTINCT"]
    test<-test[year==2100]
    test<-test[group=="Mammals"]
  }
  
  colnames(sp_dis_10km_all)[which(colnames(sp_dis_10km_all)=="N_SP.x")]<-"N_SP"
  
  sp_dis_10km_all_sub_1<-sp_dis_10km_all%>%dplyr::filter(year==2100)
  sp_dis_10km_all_sub<-sp_dis_10km_all_sub_1%>%dplyr::filter((Extinct_Type!="")|(N_type=="EXTINCT"))
  sp_dis_10km_all_sub[which(sp_dis_10km_all_sub$N_type!="EXTINCT"), "N_type"]<-"ENDANGERED"
  if (F){
    item_0<-sp_dis_10km_all_sub%>%dplyr::filter(SSP=="SSP119"&M==0&group=="Birds"&N_type=="EXTINCT")
    item_1<-sp_dis_10km_all_sub%>%dplyr::filter(SSP=="SSP119"&M==1&group=="Birds"&N_type=="EXTINCT")
    readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Dispersal/Birds/Brachypteracias_leptosomus/EC-Earth3-Veg_SSP119_5_dispersal_1.rda")
    readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Dispersal/Birds/Brachypteracias_leptosomus/EC-Earth3-Veg_SSP119_5_dispersal_1_10km.rda")
    xxx<-data.table(sp_dis_10km_all)
    
  }
  sp_dis_10km_all_sub_N<-sp_dis_10km_all_sub%>%dplyr::group_by(group, Label1, GCM, SSP, M, N_type, N_SP, TYPE)%>%
    dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
  print(table(sp_dis_10km_all_sub_N$N_SP))
  sp_dis_10km_all_sub_N$persentile<-sp_dis_10km_all_sub_N$N_SP_EXTINCT/sp_dis_10km_all_sub_N$N_SP
  
  sp_dis_10km_all_se<-sp_dis_10km_all_sub%>%dplyr::group_by(group, SSP, M, N_SP, TYPE)%>%
    dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
  sp_dis_10km_all_se$persentile<-sp_dis_10km_all_se$N_SP_EXTINCT/sp_dis_10km_all_se$N_SP
  
  if (exposure==0){
    #sp_dis_10km_all_sub_N$Label<-paste(sp_dis_10km_all_sub_N$group, " (no climate resilience)", sep="")
    sp_dis_10km_all_sub_N$exposure<-" no climate resilience"
    #sp_dis_10km_all_sub_1$Label<-paste(sp_dis_10km_all_sub_1$group, " (no climate resilience)", sep="")
    sp_dis_10km_all_sub_1$exposure<-" no climate resilience"
    #sp_dis_10km_all_se$Label<-paste(sp_dis_10km_all_se$group, " (no climate resilience)", sep="")
    sp_dis_10km_all_se$exposure<-" no climate resilience"
    
  }else{
    #sp_dis_10km_all_sub_N$Label<-paste(sp_dis_10km_all_sub_N$group, " (climate resilience)", sep="")
    sp_dis_10km_all_sub_N$exposure<-"climate resilience"
    #sp_dis_10km_all_sub_1$Label<-paste(sp_dis_10km_all_sub_1$group, " (climate resilience)", sep="")
    sp_dis_10km_all_sub_1$exposure<-"climate resilience"
    #sp_dis_10km_all_se$Label<-paste(sp_dis_10km_all_se$group, " (climate resilience)", sep="")
    sp_dis_10km_all_se$exposure<-"climate resilience"
  }
  sp_dis_10km_extinct<-bind_dplyr(sp_dis_10km_extinct, sp_dis_10km_all_sub_1)
  sp_dis_10km_all_sub_N_all<-bind_dplyr(sp_dis_10km_all_sub_N_all, sp_dis_10km_all_sub_N)
  sp_dis_10km_all_se_all<-bind_dplyr(sp_dis_10km_all_se_all, sp_dis_10km_all_se)
  if (F){
    fff<-data.table(sp_dis_10km_all_sub_N)
    fff<-fff[SSP=="SSP119"]
    fff[, .(N=.N), by=list(group, GCM, SSP, M, N_type)]
  }
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
full_combinations<-expand.grid(group=c("Birds", "Mammals"),
                               GCM=GCMs, SSP=SSPs,
                               M=c(0, 1), 
                               N_type=c("EXTINCT", "ENDANGERED"),
                               exposure=c(" no climate resilience", "climate resilience"),
                               stringsAsFactors = F)
#full_combinations$Label<-sprintf("%s (%s)", full_combinations$group, full_combinations$exposure)
full_combinations$Label1<-sprintf("%s %s", full_combinations$GCM, full_combinations$SSP)
full_combinations$exposure_number<-0
full_combinations[which(full_combinations$exposure=="climate resilience"),]$exposure_number<-5

data.table(sp_dis_10km_all_sub_N_all)[, .(N=.N), by=list(group, GCM, SSP, M, N_type, exposure)]
sp_dis_10km_all_sub_N_all<-merge(sp_dis_10km_all_sub_N_all, full_combinations, 
                            by=c("group", "GCM", "SSP", "M", "N_type", "exposure", "Label1"), all=T)
sp_dis_10km_all_sub_N_all$TYPE<-sprintf("Diversity_exposure_%d_dispersal_%d_10km_2_100km",
                                   sp_dis_10km_all_sub_N_all$exposure_number, sp_dis_10km_all_sub_N_all$M)
sp_dis_10km_all_sub_N_all[which(sp_dis_10km_all_sub_N_all$group=="Birds"), "N_SP"]<-5191
sp_dis_10km_all_sub_N_all[which(sp_dis_10km_all_sub_N_all$group=="Mammals"), "N_SP"]<-2991
sp_dis_10km_all_sub_N_all[is.na(sp_dis_10km_all_sub_N_all)]<-0



sp_mean_10km<-sp_dis_10km_all_sub_N_all%>%
  dplyr::group_by(group, SSP, M, N_type, N_SP, TYPE, exposure)%>%
  dplyr::summarise(persentile_MEAN=mean(persentile),
                   persentile_SD=sd(persentile))

sp_mean_10km$exposure<-gsub("\\(", "", sp_mean$exposure)
sp_mean_10km$exposure<-gsub("\\)", "", sp_mean$exposure)
sp_mean_10km$x_label<-as.numeric(as.factor(sp_mean_10km$SSP))

sp_mean_2_10km<-sp_mean_10km%>%dplyr::filter(N_type=="EXTINCT")


p_10km<-ggplot()+
  geom_bar(data=sp_mean_2_10km%>%dplyr::filter(M==0), 
           stat="identity", position="stack", 
           aes(y=persentile_MEAN, x=x_label-0.23, fill=factor(M), group=N_type), width=0.45, color="grey")+
  geom_bar(data=sp_mean_2_10km%>%dplyr::filter(M==1), 
           stat="identity", position="stack", 
           aes(y=persentile_MEAN, x=x_label+0.23, fill=factor(M), group=N_type), width=0.45, color="grey")+
  geom_errorbar(data=sp_mean_2_10km%>%dplyr::filter(M==0), 
                position=position_dodge(0.1), width=0.1,
                aes(ymin=persentile_MEAN-persentile_SD, 
                    ymax=persentile_MEAN+persentile_SD,
                    y=persentile_MEAN, x=x_label-0.24,
                    group=N_type)) +
  geom_errorbar(data=sp_mean_2_10km%>%dplyr::filter(M==1), 
                position=position_dodge(0.1), width=0.1,
                aes(ymin=persentile_MEAN-persentile_SD, 
                    ymax=persentile_MEAN+persentile_SD,
                    y=persentile_MEAN, x=x_label+0.24,
                    group=N_type)) +
  xlab("SSP scenario")+
  scale_x_continuous(breaks=c(1:3), labels=SSPs)+
  #ggtitle(sprintf("Distribution>%d", ttt))+
  #ylim(c(0, 1))+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  facet_grid(exposure~group)+
  scale_fill_manual(values=color_dispersal, breaks=c(0:1), 
                    labels = c("no dispersal", "with dispersal"))+
  labs(fill = "Dispersal")+
  ylab("Mean extinction proportion")
p_10km



sp_dis_100km_all_sub_N_all<-NULL
sp_dis_100km_extinct<-NULL
sp_dis_100km_all_se_all<-NULL

for (exposure in c(0, 5)){
  rda<-sprintf("../../Figures/N_Extinction/sp_dis_all_%d_100km.rda", exposure)
  print(paste("Reading", rda))
  sp_dis_100km_all<-readRDS(rda)
  
  if (F){
    test<-sp_dis_100km_all[, .(N=.N), by=list(GCM, SSP, N_type, N_SP, M, group, year)]
    test<-test[N_type=="EXTINCT"]
    test<-test[year==2100]
    test<-test[group=="Mammals"]
  }
  
  colnames(sp_dis_100km_all)[which(colnames(sp_dis_100km_all)=="N_SP.x")]<-"N_SP"
  
  sp_dis_100km_all_sub_1<-sp_dis_100km_all%>%dplyr::filter(year==2100)
  sp_dis_100km_all_sub<-sp_dis_100km_all_sub_1%>%dplyr::filter((Extinct_Type!="")|(N_type=="EXTINCT"))
  sp_dis_100km_all_sub[which(sp_dis_100km_all_sub$N_type!="EXTINCT"), "N_type"]<-"ENDANGERED"
  if (F){
    item_0<-sp_dis_100km_all_sub%>%dplyr::filter(SSP=="SSP119"&M==0&group=="Birds"&N_type=="EXTINCT")
    item_1<-sp_dis_100km_all_sub%>%dplyr::filter(SSP=="SSP119"&M==1&group=="Birds"&N_type=="EXTINCT")
    readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Dispersal/Birds/Brachypteracias_leptosomus/EC-Earth3-Veg_SSP119_5_dispersal_1.rda")
    readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Dispersal/Birds/Brachypteracias_leptosomus/EC-Earth3-Veg_SSP119_5_dispersal_1_100km.rda")
    xxx<-data.table(sp_dis_100km_all)
    
  }
  sp_dis_100km_all_sub_N<-sp_dis_100km_all_sub%>%dplyr::group_by(group, Label1, GCM, SSP, M, N_type, N_SP, TYPE)%>%
    dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
  print(table(sp_dis_100km_all_sub_N$N_SP))
  sp_dis_100km_all_sub_N$persentile<-sp_dis_100km_all_sub_N$N_SP_EXTINCT/sp_dis_100km_all_sub_N$N_SP
  
  sp_dis_100km_all_se<-sp_dis_100km_all_sub%>%dplyr::group_by(group, SSP, M, N_SP, TYPE)%>%
    dplyr::summarise(N_SP_EXTINCT=n_distinct(sp))
  sp_dis_100km_all_se$persentile<-sp_dis_100km_all_se$N_SP_EXTINCT/sp_dis_100km_all_se$N_SP
  
  if (exposure==0){
    #sp_dis_100km_all_sub_N$Label<-paste(sp_dis_100km_all_sub_N$group, " (no climate resilience)", sep="")
    sp_dis_100km_all_sub_N$exposure<-" no climate resilience"
    #sp_dis_100km_all_sub_1$Label<-paste(sp_dis_100km_all_sub_1$group, " (no climate resilience)", sep="")
    sp_dis_100km_all_sub_1$exposure<-" no climate resilience"
    #sp_dis_100km_all_se$Label<-paste(sp_dis_100km_all_se$group, " (no climate resilience)", sep="")
    sp_dis_100km_all_se$exposure<-" no climate resilience"
    
  }else{
    #sp_dis_100km_all_sub_N$Label<-paste(sp_dis_100km_all_sub_N$group, " (climate resilience)", sep="")
    sp_dis_100km_all_sub_N$exposure<-"climate resilience"
    #sp_dis_100km_all_sub_1$Label<-paste(sp_dis_100km_all_sub_1$group, " (climate resilience)", sep="")
    sp_dis_100km_all_sub_1$exposure<-"climate resilience"
    #sp_dis_100km_all_se$Label<-paste(sp_dis_100km_all_se$group, " (climate resilience)", sep="")
    sp_dis_100km_all_se$exposure<-"climate resilience"
  }
  sp_dis_100km_extinct<-bind_dplyr(sp_dis_100km_extinct, sp_dis_100km_all_sub_1)
  sp_dis_100km_all_sub_N_all<-bind_dplyr(sp_dis_100km_all_sub_N_all, sp_dis_100km_all_sub_N)
  sp_dis_100km_all_se_all<-bind_dplyr(sp_dis_100km_all_se_all, sp_dis_100km_all_se)
  if (F){
    fff<-data.table(sp_dis_100km_all_sub_N)
    fff<-fff[SSP=="SSP119"]
    fff[, .(N=.N), by=list(group, GCM, SSP, M, N_type)]
  }
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
full_combinations<-expand.grid(group=c("Birds", "Mammals"),
                               GCM=GCMs, SSP=SSPs,
                               M=c(0, 1), 
                               N_type=c("EXTINCT", "ENDANGERED"),
                               exposure=c(" no climate resilience", "climate resilience"),
                               stringsAsFactors = F)
#full_combinations$Label<-sprintf("%s (%s)", full_combinations$group, full_combinations$exposure)
full_combinations$Label1<-sprintf("%s %s", full_combinations$GCM, full_combinations$SSP)
full_combinations$exposure_number<-0
full_combinations[which(full_combinations$exposure=="climate resilience"),]$exposure_number<-5

data.table(sp_dis_100km_all_sub_N_all)[, .(N=.N), by=list(group, GCM, SSP, M, N_type, exposure)]
sp_dis_100km_all_sub_N_all<-merge(sp_dis_100km_all_sub_N_all, full_combinations, 
                                  by=c("group", "GCM", "SSP", "M", "N_type", "exposure", "Label1"), all=T)
sp_dis_100km_all_sub_N_all$TYPE<-sprintf("Diversity_exposure_%d_dispersal_%d_100km_2_100km",
                                         sp_dis_100km_all_sub_N_all$exposure_number, sp_dis_100km_all_sub_N_all$M)
sp_dis_100km_all_sub_N_all[which(sp_dis_100km_all_sub_N_all$group=="Birds"), "N_SP"]<-5191
sp_dis_100km_all_sub_N_all[which(sp_dis_100km_all_sub_N_all$group=="Mammals"), "N_SP"]<-2991
sp_dis_100km_all_sub_N_all[is.na(sp_dis_100km_all_sub_N_all)]<-0



sp_mean_100km<-sp_dis_100km_all_sub_N_all%>%
  dplyr::group_by(group, SSP, M, N_type, N_SP, TYPE, exposure)%>%
  dplyr::summarise(persentile_MEAN=mean(persentile),
                   persentile_SD=sd(persentile))

sp_mean_100km$exposure<-gsub("\\(", "", sp_mean$exposure)
sp_mean_100km$exposure<-gsub("\\)", "", sp_mean$exposure)
sp_mean_100km$x_label<-as.numeric(as.factor(sp_mean_100km$SSP))

sp_mean_2_100km<-sp_mean_100km%>%dplyr::filter(N_type=="EXTINCT")


p_100km<-ggplot()+
  geom_bar(data=sp_mean_2_100km%>%dplyr::filter(M==0), 
           stat="identity", position="stack", 
           aes(y=persentile_MEAN, x=x_label-0.23, fill=factor(M), group=N_type), width=0.45, color="grey")+
  geom_bar(data=sp_mean_2_100km%>%dplyr::filter(M==1), 
           stat="identity", position="stack", 
           aes(y=persentile_MEAN, x=x_label+0.23, fill=factor(M), group=N_type), width=0.45, color="grey")+
  geom_errorbar(data=sp_mean_2_100km%>%dplyr::filter(M==0), 
                position=position_dodge(0.1), width=0.1,
                aes(ymin=persentile_MEAN-persentile_SD, 
                    ymax=persentile_MEAN+persentile_SD,
                    y=persentile_MEAN, x=x_label-0.24,
                    group=N_type)) +
  geom_errorbar(data=sp_mean_2_100km%>%dplyr::filter(M==1), 
                position=position_dodge(0.1), width=0.1,
                aes(ymin=persentile_MEAN-persentile_SD, 
                    ymax=persentile_MEAN+persentile_SD,
                    y=persentile_MEAN, x=x_label+0.24,
                    group=N_type)) +
  xlab("SSP scenario")+
  scale_x_continuous(breaks=c(1:3), labels=SSPs)+
  #ggtitle(sprintf("Distribution>%d", ttt))+
  #ylim(c(0, 1))+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  facet_grid(exposure~group)+
  scale_fill_manual(values=color_dispersal, breaks=c(0:1), 
                    labels = c("no dispersal", "with dispersal"))+
  labs(fill = "Dispersal")+
  ylab("Mean extinction proportion")
p_100km

p<-ggpubr::ggarrange(plotlist=list(p_10km, p_100km), nrow=2, ncol=1, labels=c("10km", "100km"))

ggsave(p, filename="../../Figures/N_Extinction/compare_10km_100km.png", width=7, height=10)
