library(dplyr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
threshold<-5
dispersal<-2
if (F){
  for (dispersal in c(0, 1)){
    for (threshold in c(1, 5)){
      print(paste(dispersal, threshold))
      sp_dis_all<-readRDS(sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", threshold))
      extinct_sp<-sp_dis_all%>%dplyr::filter(year==2100)
      extinct_sp<-extinct_sp%>%dplyr::filter(N_type=="EXTINCT")
      extinct_sp<-extinct_sp%>%dplyr::filter(M==dispersal)
      extinct_sp<-extinct_sp%>%dplyr::filter(TYPE==sprintf("Diversity_%d", threshold))
      saveRDS(extinct_sp, sprintf("../../Objects_Full_species/when_where_extinction_%d/extinct_sp_%d.rda", threshold, dispersal))
    }
  }
  
}
i=1
args = commandArgs(trailingOnly=TRUE)
g<-args[1]
if (is.na(g)){
  g<-"Amphibians"
}
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-1
}



source("commonFuns/functions.r")
df<-NULL
for (dispersal in (c(0:1))){
  extinct_sp<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/extinct_sp_%d.rda", threshold, dispersal))
  extinct_sp<-extinct_sp%>%filter(group==g)
  for (i in c(1:nrow(extinct_sp))){
    print(paste(i, nrow(extinct_sp), g, "threshold=", threshold, "dispersal=", dispersal, sep=" - "))
    item<-extinct_sp[i,]
    #item$sp<-"Bunomys_fratrorum"
    #item$group<-"Mammals"
    st_dis<-readRDS(sprintf("../../Objects_Full_species/IUCN_Distribution/%s/%s.rda", 
                            item$group, item$sp))
    future_dis<-readRDS(sprintf("../../Objects_Full_species/Niche_Models/%s/%s/dispersal_%d/%s_%s_%d.rda", 
                                item$group, item$sp, threshold, item$GCM, item$SSP, dispersal))
    st_dis$group<-item$group
    st_dis$sp<-item$sp
    st_dis$GCM<-item$GCM
    st_dis$SSP<-item$SSP
    st_dis$extinct_year<-max(future_dis$YEAR)+1
    st_dis$dispersal<-dispersal
    df<-bind_dplyr(df, st_dis)
  }
}
saveRDS(df, sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
