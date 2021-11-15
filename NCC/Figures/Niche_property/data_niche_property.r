library(raster)
library(data.table)

rm(list=ls())
setDTthreads(threads=1)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Mammals"
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

mask<-readRDS("../../Objects/mask.rda")
source("commonFuns/functions.r")

nb_old<-readRDS(sprintf("../../Objects/Species_property/%s_property_with_range.rda", group))

nb<-NULL
for (i in c(1:nrow(df_list))){
  print(paste(i, nrow(df_list), group))
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$SP)
  item_old<-nb_old[sp==item$sp]
  if (nrow(item_old)==0){
    next()
  }
  if (item$N_CELL<=0){
    next()
  }
  #print(paste(Sys.time(), 1))
  #enm_folder<-sprintf("../../Objects/Niche_Models/%s/%s/dispersal_%d", group, item$sp, threshold)
  
  
  source_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, item$sp)
  #start_dis<-readRDS(sprintf("../../Objects_PNAS/Dispersal/%s/%s/occ_with_env.rda", group, item$sp))
  fit<-readRDS(sprintf("%s/fit.rda", source_folder))
  
  fit$t_max_max<-item_old$t_max_max
  fit$t_max_min<-item_old$t_max_min
  fit$t_min_max<-item_old$t_min_max
  fit$t_min_min<-item_old$t_min_min
  fit$pr_max<-item_old$pr_max
  fit$pr_min<-item_old$pr_min
  
  
  fit$sp<-item$sp
  
  fit$nb_bio1_sd<-fit$range_bio1_sd_max-fit$range_bio1_sd_min
  fit$nb_bio12_sd<-fit$range_bio12_sd_max-fit$range_bio12_sd_min
  nb<-bind(nb, fit)
}
saveRDS(nb, sprintf("../../Objects/Species_property/%s_property.rda", group))
  