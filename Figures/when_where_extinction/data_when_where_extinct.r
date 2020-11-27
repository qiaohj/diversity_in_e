library(dplyr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
threshold<-1
if (F){
  sp_dis_all<-readRDS(sprintf("../../Figures/N_SPECIES_%d/sp_dis_all.rda", threshold))
  extinct_sp<-sp_dis_all%>%dplyr::filter(year==2100)
  extinct_sp<-extinct_sp%>%dplyr::filter(N_type=="EXTINCT")
  extinct_sp<-extinct_sp%>%dplyr::filter(M==0)
  extinct_sp<-extinct_sp%>%dplyr::filter(TYPE==sprintf("Diversity_%d", threshold))
  saveRDS(extinct_sp, sprintf("../../Objects/when_where_extinction_%d/extinct_sp.rda", threshold))
}
extinct_sp<-readRDS(sprintf("../../Objects/when_where_extinction_%d/extinct_sp.rda", threshold))
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
  if (threshold==5){
    future_dis<-readRDS(sprintf("../../Objects/Niche_Models/%s/%s/dispersal_%d/%s_%s_0.rda", 
                          item$group, item$sp, threshold, item$GCM, item$SSP))
  }else{
    future_dis<-readRDS(sprintf("../../Objects/Niche_Models/%s/%s/dispersal/%s_%s_0.rda", 
                                item$group, item$sp, item$GCM, item$SSP))
  }
  st_dis$group<-item$group
  st_dis$sp<-item$sp
  st_dis$GCM<-item$GCM
  st_dis$SSP<-item$SSP
  st_dis$extinct_year<-max(future_dis$YEAR)+1
  df<-bind_dplyr(df, st_dis)
}
saveRDS(df, sprintf("../../Objects/when_where_extinction_%d/%s.rda", threshold, g))
