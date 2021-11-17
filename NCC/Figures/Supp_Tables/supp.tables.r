library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
bird_df<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")

bird_nb<-readRDS("../../Objects/Species_property/Birds_property.rda")

bird_nb<-bird_nb[!is.na(bio1_max)]
dim(bird_nb)

colnames(bird_nb)
cols<-c("range_bio1_sd_min", "range_bio1_sd_max",
        "range_bio5_sd_min", "range_bio5_sd_max",
        "range_bio6_sd_min", "range_bio6_sd_max",
        "range_bio12_sd_min", "range_bio12_sd_max",
        "range_bio13_sd_min", "range_bio13_sd_max",
        "range_bio14_sd_min", "range_bio14_sd_max",
        "N_CELL", "sp")
bird_nb<-bird_nb[, ..cols]

cols<-c("Order", "family", "sp", "HWI", "log_body_mass", "Diet", "Migration_3", "estimated_disp")

colnames(bird_df)

bird_df<-bird_df[, ..cols]

bird_property<-merge(bird_df, bird_nb, by="sp")
bird_property
write.csv(bird_property, "../../Objects/Supp.Tables/bird_property.csv", row.names=F)

when_extinct_bird_0<-readRDS("../../Objects/when_where_extinction_exposure_0/Birds.rda")
cols<-c("sp", "GCM", "SSP", "extinct_year", "dispersal")
when_extinct_bird_0<-unique(when_extinct_bird_0[, ..cols])
when_extinct_bird_0$exposure<-0
when_extinct_bird_5<-readRDS("../../Objects/when_where_extinction_exposure_5/Birds.rda")
when_extinct_bird_5<-unique(when_extinct_bird_5[, ..cols])
when_extinct_bird_5$exposure<-5

when_extinct_bird<-rbindlist(list(when_extinct_bird_0, when_extinct_bird_5))
cols<-c("sp", "Order", "family")
bird_property_item<-bird_property[, ..cols]

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
SSP_str<-SSPs[1]
GCM_str<-GCMs[1]
disp<-0
expo<-0
final_df_bird<-list()
for (SSP_str in SSPs){
  for (GCM_str in GCMs){
    for (disp in c(0, 1)){
      for (expo in c(0, 5)){
        item<-when_extinct_bird[(SSP==SSP_str)&(GCM==GCM_str)&(dispersal==disp)&(exposure==expo)]
        item<-merge(bird_property_item, item, by="sp", all=T)
        item$GCM<-GCM_str
        item$SSP<-SSP_str
        item$exposure<-expo
        item$dispersal<-disp
        final_df_bird[[paste(SSP_str, GCM_str, disp, expo)]]<-item
      }
    }
  }  
}
final_df_bird<-rbindlist(final_df_bird)
write.csv(final_df_bird, "../../Objects/Supp.Tables/bird_extinction.csv", row.names=F)

##For mammals

mammal_df<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")

mammal_nb<-readRDS("../../Objects/Species_property/Mammals_property.rda")

mammal_nb<-mammal_nb[!is.na(bio1_max)]
dim(mammal_nb)

colnames(mammal_nb)
cols<-c("range_bio1_sd_min", "range_bio1_sd_max",
        "range_bio5_sd_min", "range_bio5_sd_max",
        "range_bio6_sd_min", "range_bio6_sd_max",
        "range_bio12_sd_min", "range_bio12_sd_max",
        "range_bio13_sd_min", "range_bio13_sd_max",
        "range_bio14_sd_min", "range_bio14_sd_max",
        "N_CELL", "sp")
mammal_nb<-mammal_nb[, ..cols]

cols<-c("Order", "family", "sp", "ForStrat", "log_body_mass", "Diet", "estimated_disp")

colnames(mammal_df)

mammal_df<-mammal_df[, ..cols]

mammal_property<-merge(mammal_df, mammal_nb, by="sp")
mammal_property
write.csv(mammal_property, "../../Objects/Supp.Tables/mammal_property.csv", row.names=F)

when_extinct_mammal_0<-readRDS("../../Objects/when_where_extinction_exposure_0/Mammals.rda")
cols<-c("sp", "GCM", "SSP", "extinct_year", "dispersal")
when_extinct_mammal_0<-unique(when_extinct_mammal_0[, ..cols])
when_extinct_mammal_0$exposure<-0
when_extinct_mammal_5<-readRDS("../../Objects/when_where_extinction_exposure_5/Mammals.rda")
when_extinct_mammal_5<-unique(when_extinct_mammal_5[, ..cols])
when_extinct_mammal_5$exposure<-5

when_extinct_mammal<-rbindlist(list(when_extinct_mammal_0, when_extinct_mammal_5))
cols<-c("sp", "Order", "family")
mammal_property_item<-mammal_property[, ..cols]

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
SSP_str<-SSPs[1]
GCM_str<-GCMs[1]
disp<-0
expo<-0
final_df_mammal<-list()
for (SSP_str in SSPs){
  for (GCM_str in GCMs){
    for (disp in c(0, 1)){
      for (expo in c(0, 5)){
        item<-when_extinct_mammal[(SSP==SSP_str)&(GCM==GCM_str)&(dispersal==disp)&(exposure==expo)]
        item<-merge(mammal_property_item, item, by="sp", all=T)
        item$GCM<-GCM_str
        item$SSP<-SSP_str
        item$exposure<-expo
        item$dispersal<-disp
        final_df_mammal[[paste(SSP_str, GCM_str, disp, expo)]]<-item
      }
    }
  }  
}
final_df_mammal<-rbindlist(final_df_mammal)
write.csv(final_df_mammal, "../../Objects/Supp.Tables/mammal_extinction.csv", row.names=F)


