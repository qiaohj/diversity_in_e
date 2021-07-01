library(dplyr)
library(data.table)
rm(list=ls())
threshold<-1
ttt<-2
g<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")
source("commonFuns/functions.r")
all_result<-NULL
for (threshold in c(1,5)){
  when_extinct<-NULL
  
  for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
    sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
    sp_list<-sp_list[which(sp_list$area>ttt),]
    
    df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
    sp_list$sp2<-gsub(" ", "_", sp_list$sp)
    df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
    when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
    when_extinct<-when_extinct%>%dplyr::filter(!is.infinite(extinct_year))
    for (i in c(1:nrow(when_extinct))){
      print(paste(threshold, g, i, nrow(when_extinct)))
      item<-when_extinct[i,]
      
      target_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", g, item$sp)
      fit<-readRDS(sprintf("%s/fit.rda", target_folder))
      dis<-readRDS(sprintf("%s/dispersal_%d/%s_%s_%d.rda",
                           target_folder, threshold, item$GCM, item$SSP, item$dispersal))
      dis<-dis%>%dplyr::filter(YEAR==(item$extinct_year-1))
      env_layer<-env_layers[[sprintf("%s_%s_%d", item$GCM, item$SSP, item$extinct_year)]]
      env_layer<-env_layer%>%dplyr::filter(mask_index %in% dis$mask_index)
      env_layer$is_temp_in<-between(env_layer$TEMP_MAX, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)&
        between(env_layer$TEMP_MIN, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)
      env_layer$is_prec_in<-between(env_layer$PR, fit$range_PR_sd_min, fit$range_PR_sd_max)
      for (nn in names(item)){
        env_layer[, nn]<-item[, nn]
      }
      env_layer$threshold<-threshold
      all_result<-bind_dplyr(all_result, env_layer)
    }
  }
}
saveRDS(all_result, "../../Objects_Full_species/why_extinct/why_extinct.rda")
