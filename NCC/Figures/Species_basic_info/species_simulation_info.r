library(dplyr)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
ttt<-2
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
full_df_all<-NULL
for (threshold in c(1,5)){
  for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
    sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
    sp_list$invoved<-ifelse(sp_list$area>ttt, "YES", "NO")
    df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
    sp2<-gsub(" ", "_", sp_list[which(sp_list$invoved=="YES"), ]$sp)
    full_df<-expand.grid(group=g, sp=sp2, GCM=GCMs, SSP=SSPs, 
                         dispersal=c("with dispersal", "no dispersal"),
                         exposure=c("no exposure", "5-year exposure"),
                         stringsAsFactors = F)
    df<-df%>%dplyr::filter(sp %in% sp2)%>%
      distinct(group, sp, GCM, SSP, extinct_year, dispersal)
    df$dispersal<-ifelse(df$dispersal==1, "with dispersal", "no dispersal")
    df$exposure<-ifelse(threshold==1, "no exposure", "5-year exposure")
    
    if (F){
      df<-df%>%group_by(group,sp,GCM,SSP,dispersal,exposure)%>%
        dplyr::summarise(max_extinct_year=max(extinct_year))
      df%>%filter((sp=="Ablepharus_rueppellii")&(GCM=="MRI-ESM2-0")&(SSP=="SSP585"))
      full_df_se<-full_df%>%dplyr::group_by(group,sp,GCM,SSP,dispersal,exposure)%>%
        dplyr::summarise(N=n())
      full_df_se[which(full_df_se$N>1),]
    }
    full_df<-left_join(full_df, df, by=c("group", "sp", "GCM", "SSP", "dispersal", "exposure"))
    
    full_df_all<-bind_dplyr(full_df_all, full_df)
    
  }
}

full_df_all[is.na(full_df_all$extinct_year), "extinct_year"]<-9999
full_df_all$sp<-gsub("_", " ", full_df_all$sp)
write.table(full_df_all, "../../Figures_Full_species/full_sp_list.csv", row.names = F, sep=",")


full_df_all$is_extinct<-ifelse(full_df_all$extinct_year==9999, "no", "yes")

full_df_all%>%group_by(group, is_extinct, SSP, dispersal, exposure)%>%
  dplyr::summarise(N=n())
