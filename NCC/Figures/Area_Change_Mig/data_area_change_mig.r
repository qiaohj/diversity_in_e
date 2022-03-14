library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Birds"
}

esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")

env<-NULL
print(group)

if (group=="Birds"){
  group_df<-readRDS("../../Data/Birds/bird_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  #group_disp$estimated_disp
  #group_disp2<-readRDS("../../Objects_PNAS/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  #dd<-merge(data.frame(group_disp[, c("iucn_name", "estimated_disp")]), 
  #          data.frame(group_disp2[, c("iucn_name", "estimated_disp")]), 
  #          by.x="iucn_name", by.y="iucn_name", all.x=T, all.y=F)
  group_full<-merge(group_df, group_disp, by.x="SCINAME", by.y="iucn_name", all=F)
  
}else{
  group_df<-readRDS("../../Data/Mammals/mammal_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  group_full<-merge(group_df, group_disp, by.x="binomial", by.y="Scientific", all=F)
  
  colnames(group_full)[1]<-"SCINAME"
  colnames(group_full)[27]<-"Shape_Area"
  colnames(group_disp)[1]<-"iucn_name"
  
}

PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(2)

predict_range<-c(2021:2100)

bi<-"Phylloscopus borealoides"

group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]
coms<-expand.grid(exposure_threshold=c(0, 5), dispersal=c(0, 1))
i=1
j=1
item_str<-esm_ssp[1]
all<-list()
for (i in 1:length(group_full_sum_area$SCINAME)) {
  print(paste(i, nrow(group_full_sum_area)))
  for (j in 1:nrow(coms)){
    dispersal<-coms[j, "dispersal"]
    exposure_threshold<-coms[j, "exposure_threshold"]
    bi<-group_full_sum_area$SCINAME[i]
    
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
    fit_str<-sprintf("%s/fit_seasonal_2.rda", target_folder)
    if (!file.exists(fit_str)){
      
      next()
    }
    
    for (item_str in esm_ssp){
      target<-sprintf("%s/%s_%d_dispersal_%d_10km_seasonal_2.rda", target_folder, item_str,
                      exposure_threshold, dispersal)
      
      if (file.exists(target)){
        res<-"10km"
        df_mig<-readRDS(target)
        if (!file.exists(sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, item_str,
                                 exposure_threshold, dispersal))){
          print("skip 2")
          next()
        }
        df_all<-readRDS(sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, item_str,
                                exposure_threshold, dispersal))
      }else{
        res<-"100km"
        target<-sprintf("%s/%s_%d_dispersal_%d_100km_seasonal_2.rda", target_folder, item_str,
                        exposure_threshold, dispersal)
        if (!file.exists(target)){
          print("skip 3")
          next()
        }
        df_mig<-readRDS(target)
        if (!file.exists(sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                                 exposure_threshold, dispersal))){
          print("skip 4")
          next()
        }
        df_all<-readRDS(sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                                exposure_threshold, dispersal))
        
      }
      if (length(df_mig)>0){
        df_mig<-rbindlist(df_mig)
        df_mig_se<-df_mig[, .(N_mig=.N), by=list(YEAR)]
        df_all<-rbindlist(df_all)
        df_allg_se<-df_all[, .(N_all=.N), by=list(YEAR)]
        df<-merge(df_mig_se, df_allg_se, by="YEAR", all=T)
        df$dispersal<-dispersal
        df$exposure<-exposure_threshold
        eee<-strsplit(item_str, "_")[[1]]
        df$esm<-eee[1]
        df$ssp<-eee[2]
        df$sp<-bi
        df$res<-res
        all[[length(all)+1]]<-df
      }
    }
  }
}

all<-rbindlist(all)
all[is.na(N_all)]$N_all<-0
all[is.na(N_mig)]$N_mig<-0
saveRDS(all, "../../Objects/mig_bird_log.rda")
