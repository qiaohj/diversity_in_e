library(dplyr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  sp_dis_all<-readRDS("../../Figures/N_SPECIES_5/sp_dis_all.rda")
  extinct_sp<-sp_dis_all%>%dplyr::filter(year==2100)
  extinct_sp<-extinct_sp%>%dplyr::filter(N_type=="EXTINCT")
  extinct_sp<-extinct_sp%>%dplyr::filter(M==0)
  extinct_sp<-extinct_sp%>%dplyr::filter(TYPE=="Diversity_5")
  saveRDS(extinct_sp, "../../Objects/when_where_extinction_5/extinct_sp.rda")
}
extinct_sp<-readRDS("../../Objects/when_where_extinction_5/extinct_sp.rda")
i=1
args = commandArgs(trailingOnly=TRUE)
g<-args[1]
if (is.na(g)){
  g<-"Amphibians"
}
extinct_sp<-extinct_sp%>%filter(group==g)
source("commonFuns/functions.r")
df<-NULL
for (i in c(1:nrow(extinct_sp))){
  print(paste(i, nrow(extinct_sp), g, sep=" - "))
  item<-extinct_sp[i,]
  #item$sp<-"Bunomys_fratrorum"
  #item$group<-"Mammals"
  st_dis<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", 
                          item$group, item$sp))
  future_dis<-readRDS(sprintf("../../Objects/Niche_Models/%s/%s/dispersal_5/%s_%s_0.rda", 
                          item$group, item$sp, item$GCM, item$SSP))
  st_dis$group<-item$group
  st_dis$sp<-item$sp
  st_dis$GCM<-item$GCM
  st_dis$SSP<-item$SSP
  st_dis$extinct_year<-max(future_dis$YEAR)+1
  df<-bind_dplyr(df, st_dis)
}
saveRDS(df, sprintf("../../Objects/when_where_extinction_5/%s.rda", g))
