library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
df_sp_list<-list()
for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list)


table(df_sp_list$group)
table(df_sp_list$area)
table(df_sp_list[which(df_sp_list$area>2), "group"])
table(df_sp_list[which(df_sp_list$area>0), "group"])
