library(raster)
library(data.table)

rm(list=ls())
setDTthreads(threads=1)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
i=1
j=1
k=2

dispersals<-c(0:1)
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

mask<-readRDS("../../Objects_Full_species/mask.rda")
source("commonFuns/functions.r")
nb<-NULL
for (i in c(1:nrow(df_list))){
  print(paste(i, nrow(df_list), group))
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  
  if (item$area<=0){
    next()
  }
  #print(paste(Sys.time(), 1))
  #enm_folder<-sprintf("../../Objects/Niche_Models/%s/%s/dispersal_%d", group, item$sp, threshold)
  
  
  source_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", group, item$sp)
  start_dis<-readRDS(sprintf("%s/occ_with_env.rda", source_folder))
  fit<-readRDS(sprintf("%s/fit.rda", source_folder))
  
  fit$t_max_max<-max(start_dis$TEMP_MAX, na.rm = T)
  fit$t_max_min<-min(start_dis$TEMP_MAX, na.rm = T)
  fit$t_min_max<-max(start_dis$TEMP_MIN, na.rm = T)
  fit$t_min_min<-min(start_dis$TEMP_MIN, na.rm = T)
  fit$pr_max<-max(start_dis$PR, na.rm = T)
  fit$pr_min<-min(start_dis$PR, na.rm = T)
  
  
  fit$sp<-item$sp
  
  fit$nb_TEMP_sd<-fit$range_TEMP_sd_max-fit$range_TEMP_sd_min
  fit$nb_PR_sd<-fit$range_PR_sd_max-fit$range_PR_sd_min
  nb<-bind(nb, fit)
}
saveRDS(nb, sprintf("../../Objects_Full_species/Species_property/%s_property.rda", group))
  