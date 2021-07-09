library(dplyr)
library(data.table)
rm(list=ls())
exposure<-0
g<-"Mammals"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")
source("commonFuns/functions.r")
all_result<-NULL
for (exposure in c(0,5)){
  when_extinct<-NULL
  
  for (g in c("Birds", "Mammals")){
    sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
    
    df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
    sp_list$sp<-sp_list$SP
    sp_list$sp2<-gsub(" ", "_", sp_list$sp)
    df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
    when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
    when_extinct<-when_extinct%>%dplyr::filter(!is.infinite(extinct_year))
    i=1
    for (i in c(1:nrow(when_extinct))){
      print(paste(exposure, g, i, nrow(when_extinct)))
      item<-when_extinct[i,]
      
      target_folder<-sprintf("../../Objects/Dispersal/%s/%s", g, item$sp)
      fit<-readRDS(sprintf("%s/fit.rda", target_folder))
      dis<-readRDS(sprintf("%s/%s_%s_%d_dispersal_%d.rda", 
                           target_folder, item$GCM, item$SSP, exposure, item$dispersal))
      dis<-rbindlist(dis)
      if (is.null(dis)){
        next()
      }
      if (nrow(dis)==0){
        next()
      }
      
      dis<-dis%>%dplyr::filter(YEAR==(item$extinct_year-1))
      env_layer<-env_layers[[sprintf("%s_%s", item$GCM, item$SSP)]]
      env_layer<-env_layer%>%dplyr::filter(year==item$extinct_year)
      env_layer<-env_layer%>%dplyr::filter(mask_100km %in% dis$mask_100km)
      env_layer$is_bio1<-between(env_layer$bio1, fit$range_bio1_sd_min, fit$range_bio1_sd_max)
      env_layer$is_bio5<-between(env_layer$bio5, fit$range_bio5_sd_min, fit$range_bio5_sd_max)
      env_layer$is_bio6<-between(env_layer$bio6, fit$range_bio6_sd_min, fit$range_bio6_sd_max)
      env_layer$is_bio12<-between(env_layer$bio12, fit$range_bio12_sd_min, fit$range_bio12_sd_max)
      env_layer$is_bio13<-between(env_layer$bio13, fit$range_bio13_sd_min, fit$range_bio13_sd_max)
      env_layer$is_bio14<-between(env_layer$bio14, fit$range_bio14_sd_min, fit$range_bio14_sd_max)
      env_layer$upper_bio1<-env_layer$bio1>fit$range_bio1_sd_max
      env_layer$lower_bio1<-env_layer$bio1<fit$range_bio1_sd_min
      env_layer$upper_bio5<-env_layer$bio5>fit$range_bio5_sd_max
      env_layer$lower_bio5<-env_layer$bio5<fit$range_bio5_sd_min
      env_layer$upper_bio6<-env_layer$bio6>fit$range_bio6_sd_max
      env_layer$lower_bio6<-env_layer$bio6<fit$range_bio6_sd_min
      env_layer$upper_bio12<-env_layer$bio12>fit$range_bio12_sd_max
      env_layer$lower_bio12<-env_layer$bio12<fit$range_bio12_sd_min
      env_layer$upper_bio13<-env_layer$bio13>fit$range_bio13_sd_max
      env_layer$lower_bio13<-env_layer$bio13<fit$range_bio13_sd_min
      env_layer$upper_bio14<-env_layer$bio14>fit$range_bio14_sd_max
      env_layer$lower_bio14<-env_layer$bio14<fit$range_bio14_sd_min
      
      env_layer<-as_tibble(env_layer)
      item<-data.frame(item)
      for (nn in names(item)){
        env_layer[, nn]<-item[, nn]
      }
      env_layer$exposure<-exposure
      all_result<-bind_dplyr(all_result, env_layer)
    }
  }
}
saveRDS(all_result, "../../Objects/why_extinct/why_extinct.rda")
