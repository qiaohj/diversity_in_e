library(dplyr)
library(data.table)
library(raster)
rm(list=ls())
exposure<-0
g<-"Mammals"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")
source("commonFuns/functions.r")
points_10km<-readRDS("../../Raster/points_10km.rda")
colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-data.table(rasterToPoints(mask_100km))
colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
points_10km<-merge(points_10km, points_100km, by="mask_100km")

all_result<-NULL
for (exposure in c(0,5)){
  when_extinct<-NULL
  
  for (g in c("Birds", "Mammals")){
    sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
    
    df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s_10km.rda", exposure, g))
    df<-rbindlist(df)
    sp_list$sp<-sp_list$SP
    sp_list$sp2<-gsub(" ", "_", sp_list$sp)
    df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
    when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
    when_extinct<-when_extinct%>%dplyr::filter(!is.infinite(extinct_year))
    i=1
    for (i in c(35:nrow(when_extinct))){
      print(paste(exposure, g, i, nrow(when_extinct)))
      item<-when_extinct[i,]
      
      target_folder<-sprintf("../../Objects/Dispersal/%s/%s", g, item$sp)
      if (!file.exists(sprintf("%s/fit.rda", target_folder))){
        next()
      }
      fit<-readRDS(sprintf("%s/fit.rda", target_folder))
      f<-sprintf("%s/%s_%s_%d_dispersal_%d_10km.rda", 
                 target_folder, item$GCM, item$SSP, exposure, item$dispersal)
      res<-"10km"
      if (!file.exists(f)){
        f<-sprintf("%s/%s_%s_%d_dispersal_%d.rda", 
                   target_folder, item$GCM, item$SSP, exposure, item$dispersal)
        res<-"100km"
      }
      if (!file.exists(f)){
        next()
      }
      dis<-readRDS(f)
      if (length(dis)==0){
        next()
      }
      dis<-rbindlist(dis)
      if (res=="10km"){
        dis<-merge(dis, points_10km, by.x="mask_10km", by.y="mask_10km")
        dis<-dis[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR", "exposure")]
        colnames(dis)[1:3]<-c("x", "y", "mask_100km")
      }
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
saveRDS(all_result, "../../Objects/why_extinct/why_extinct_10km.rda")
