library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
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
  bird_property$sp<-gsub("_", " ", bird_property$sp)
  group_full_filter<-data.table(read.csv("../../Objects/species_list.csv", head=T, stringsAsFactors = F))
  bird_property<-merge(bird_property, group_full_filter, by=c("sp"), all.x=T)
  item1<-bird_property[sp=="Cuculus_saturatus"][1]
  item2<-bird_property[sp=="Corapipo leucorrhoa"][2]
  item3<-bird_property[sp=="Tangara cyanoptera"][2]
  group_full<-group_full[!(sp %in% unique_sp_t$sp)]
  group_full<-rbindlist(list(group_full, item1, item2, item3))
  
  table(group_full_filter$group)
}
group_full_filter<-data.table(read.csv("../../Objects/species_list.csv", head=T, stringsAsFactors = F))
bird<-group_full_filter[group=="Birds"]
bird$ForStrat<-NULL
i=1
bird$range_bio1_sd_min<--9999
bird$range_bio1_sd_max<--9999
bird$range_bio5_sd_min<--9999
bird$range_bio5_sd_max<--9999
bird$range_bio6_sd_min<--9999
bird$range_bio6_sd_max<--9999
bird$range_bio12_sd_min<--9999
bird$range_bio12_sd_max<--9999
bird$range_bio13_sd_min<--9999
bird$range_bio13_sd_max<--9999
bird$range_bio14_sd_min<--9999
bird$range_bio14_sd_max<--9999
bird$N_CELL<--9999
for (i in c(1:nrow(bird))){
  fit<-readRDS(sprintf("../../Objects/Dispersal/Birds/%s/fit.rda", gsub(" ", "_", bird[i]$sp)))
  bird[i]$range_bio1_sd_min<-fit$range_bio1_sd_min
  bird[i]$range_bio1_sd_max<-fit$range_bio1_sd_max
  bird[i]$range_bio5_sd_min<-fit$range_bio5_sd_min
  bird[i]$range_bio5_sd_max<-fit$range_bio5_sd_max
  bird[i]$range_bio6_sd_min<-fit$range_bio6_sd_min
  bird[i]$range_bio6_sd_max<-fit$range_bio6_sd_max
  bird[i]$range_bio12_sd_min<-fit$range_bio12_sd_min
  bird[i]$range_bio12_sd_max<-fit$range_bio12_sd_max
  bird[i]$range_bio13_sd_min<-fit$range_bio13_sd_min
  bird[i]$range_bio13_sd_max<-fit$range_bio13_sd_max
  bird[i]$range_bio14_sd_min<-fit$range_bio14_sd_min
  bird[i]$range_bio14_sd_max<-fit$range_bio14_sd_max
  bird[i]$N_CELL<-fit$N_CELL
}
bird_df<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")
col<-c("Order", "family", "sp")
bird_df<-bird_df[,..col]
bird_df$sp<-gsub("_", " ", bird_df$sp)
bird_with_family<-merge(bird, bird_df, by="sp", all.x=T)
bird_with_family[is.na(Order)]$Order<-"Unknown"
bird_with_family[is.na(family)]$family<-"Unknown"
write.csv(bird_with_family, "../../Objects/Supp.Tables/bird_property_10km.csv", row.names=F)

when_extinct_bird_0<-readRDS("../../Objects/when_where_extinction_exposure_0/Birds_10km.rda")
when_extinct_bird_0<-rbindlist(when_extinct_bird_0)
cols<-c("sp", "GCM", "SSP", "extinct_year", "dispersal")
when_extinct_bird_0<-unique(when_extinct_bird_0[, ..cols])
when_extinct_bird_0$exposure<-0
when_extinct_bird_5<-readRDS("../../Objects/when_where_extinction_exposure_5/Birds_10km.rda")
when_extinct_bird_5<-rbindlist(when_extinct_bird_5)
when_extinct_bird_5<-unique(when_extinct_bird_5[, ..cols])
when_extinct_bird_5$exposure<-5

when_extinct_bird<-rbindlist(list(when_extinct_bird_0, when_extinct_bird_5))
cols<-c("sp", "Order", "family")
bird_property_item<-bird_with_family[, ..cols]
when_extinct_bird$sp<-gsub("_", " ", when_extinct_bird$sp)
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
colnames(final_df_bird)[4]<-"ESM"
colnames(final_df_bird)[8]<-"climate_resilience"
final_df_bird[!(sp %in% bird_with_family$sp)]
final_df_bird[, .(N=.N), by=list(ESM, SSP, dispersal, climate_resilience)]

write.csv(final_df_bird, "../../Objects/Supp.Tables/bird_extinction_10km.csv", row.names=F)

##For mammals
if (F){
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
  colnames(final_df_mammal)[4]<-"ESM"
  colnames(final_df_mammal)[8]<-"climate_resilience"
  
  write.csv(final_df_mammal, "../../Objects/Supp.Tables/mammal_extinction.csv", row.names=F)
  
}


group_full_filter<-data.table(read.csv("../../Objects/species_list.csv", head=T, stringsAsFactors = F))
mammal<-group_full_filter[group=="Mammals"]

i=1
mammal$range_bio1_sd_min<--9999
mammal$range_bio1_sd_max<--9999
mammal$range_bio5_sd_min<--9999
mammal$range_bio5_sd_max<--9999
mammal$range_bio6_sd_min<--9999
mammal$range_bio6_sd_max<--9999
mammal$range_bio12_sd_min<--9999
mammal$range_bio12_sd_max<--9999
mammal$range_bio13_sd_min<--9999
mammal$range_bio13_sd_max<--9999
mammal$range_bio14_sd_min<--9999
mammal$range_bio14_sd_max<--9999
mammal$N_CELL<--9999

for (i in c(1:nrow(mammal))){
  fit<-readRDS(sprintf("../../Objects/Dispersal/Mammals/%s/fit.rda", gsub(" ", "_", mammal[i]$sp)))
  mammal[i]$range_bio1_sd_min<-fit$range_bio1_sd_min
  mammal[i]$range_bio1_sd_max<-fit$range_bio1_sd_max
  mammal[i]$range_bio5_sd_min<-fit$range_bio5_sd_min
  mammal[i]$range_bio5_sd_max<-fit$range_bio5_sd_max
  mammal[i]$range_bio6_sd_min<-fit$range_bio6_sd_min
  mammal[i]$range_bio6_sd_max<-fit$range_bio6_sd_max
  mammal[i]$range_bio12_sd_min<-fit$range_bio12_sd_min
  mammal[i]$range_bio12_sd_max<-fit$range_bio12_sd_max
  mammal[i]$range_bio13_sd_min<-fit$range_bio13_sd_min
  mammal[i]$range_bio13_sd_max<-fit$range_bio13_sd_max
  mammal[i]$range_bio14_sd_min<-fit$range_bio14_sd_min
  mammal[i]$range_bio14_sd_max<-fit$range_bio14_sd_max
  mammal[i]$N_CELL<-fit$N_CELL
}
mammal_df<-readRDS("../../Objects/IUCN_List/Mammals_df.rda")

mammal_df<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")
col<-c("Order", "family", "sp")
mammal_df<-mammal_df[,..col]
mammal_df$sp<-gsub("_", " ", mammal_df$sp)
mammal_with_family<-merge(mammal, mammal_df, by="sp", all.x=T)
mammal_with_family[is.na(Order)]$Order<-"Unknown"
mammal_with_family[is.na(family)]$family<-"Unknown"
mammal_with_family$Migration<-NULL
mammal_with_family$HWI<-NULL
mammal_with_family$in_ebird<-NULL
mammal_with_family$is_migratory_bird<-NULL

write.csv(mammal_with_family, "../../Objects/Supp.Tables/mammal_property_10km.csv", row.names=F)

when_extinct_mammal_0<-readRDS("../../Objects/when_where_extinction_exposure_0/Mammals_10km.rda")
when_extinct_mammal_0<-rbindlist(when_extinct_mammal_0)
cols<-c("sp", "GCM", "SSP", "extinct_year", "dispersal")
when_extinct_mammal_0<-unique(when_extinct_mammal_0[, ..cols])
when_extinct_mammal_0$exposure<-0
when_extinct_mammal_5<-readRDS("../../Objects/when_where_extinction_exposure_5/Mammals_10km.rda")
when_extinct_mammal_5<-rbindlist(when_extinct_mammal_5)
when_extinct_mammal_5<-unique(when_extinct_mammal_5[, ..cols])
when_extinct_mammal_5$exposure<-5

when_extinct_mammal<-rbindlist(list(when_extinct_mammal_0, when_extinct_mammal_5))
cols<-c("sp", "Order", "family")
mammal_property_item<-mammal_with_family[, ..cols]
when_extinct_mammal$sp<-gsub("_", " ", when_extinct_mammal$sp)
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
colnames(final_df_mammal)[4]<-"ESM"
colnames(final_df_mammal)[8]<-"climate_resilience"
final_df_mammal[!(sp %in% mammal_with_family$sp)]
final_df_mammal[, .(N=.N), by=list(ESM, SSP, dispersal, climate_resilience)]

write.csv(final_df_mammal, "../../Objects/Supp.Tables/mammal_extinction_10km.csv", row.names=F)
