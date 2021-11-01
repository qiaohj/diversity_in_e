library(raster)
library(data.table)

rm(list=ls())
setDTthreads(threads=1)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Birds"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
i=1
j=1
k=2

dispersals<-c(0:1)
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df_100km.rda")
env_layers<-env_layers[year==2020]
mask<-readRDS("../../Objects/mask.rda")
source("commonFuns/functions.r")
nb<-NULL
target<-"Dispersal"
i=1
for (target in c("1850", "1970")){
  
  for (i in c(1:nrow(df_list))){
    print(paste(i, nrow(df_list), group))
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$SP)
    
    #print(paste(Sys.time(), 1))
    #enm_folder<-sprintf("../../Objects/Niche_Models/%s/%s/dispersal_%d", group, item$sp, threshold)
    
    
    source_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, item$sp)
    start_dis<-readRDS(sprintf("%s/initial_disp_exposure_0_dispersal_0.rda", source_folder))
    start_dis<-merge(start_dis, env_layers, by=c("x", "y", "mask_100km"))
    if (target=="1850"){
      fit<-readRDS(sprintf("%s/fit.rda", source_folder))
    }else{
      fit<-readRDS(sprintf("%s/fit_1970.rda", source_folder))
    }
    
    fit$bio1_max<-max(start_dis$bio1, na.rm = T)
    fit$bio1_min<-min(start_dis$bio1, na.rm = T)
    fit$bio5_max<-max(start_dis$bio5, na.rm = T)
    fit$bio5_min<-min(start_dis$bio5, na.rm = T)
    fit$bio6_max<-max(start_dis$bio6, na.rm = T)
    fit$bio6_min<-min(start_dis$bio6, na.rm = T)
    fit$bio12_max<-max(start_dis$bio12, na.rm = T)
    fit$bio12_min<-min(start_dis$bio12, na.rm = T)
    fit$bio13_max<-max(start_dis$bio13, na.rm = T)
    fit$bio13_min<-min(start_dis$bio13, na.rm = T)
    fit$bio14_max<-max(start_dis$bio14, na.rm = T)
    fit$bio14_min<-min(start_dis$bio14, na.rm = T)
    
    
    
    fit$sp<-item$sp
    
    fit$nb_bio1_sd<-fit$range_bio1_sd_max-fit$range_bio1_sd_min
    fit$nb_bio5_sd<-fit$range_bio5_sd_max-fit$range_bio5_sd_min
    fit$nb_bio6_sd<-fit$range_bio6_sd_max-fit$range_bio6_sd_min
    fit$nb_bio12_sd<-fit$range_bio12_sd_max-fit$range_bio12_sd_min
    fit$nb_bio13_sd<-fit$range_bio13_sd_max-fit$range_bio13_sd_min
    fit$nb_bio14_sd<-fit$range_bio14_sd_max-fit$range_bio14_sd_min
    fit$target<-target
    nb<-bind(nb, fit)
  }
}
saveRDS(nb, sprintf("../../Objects/Species_property/%s_property_compared.rda", group))
