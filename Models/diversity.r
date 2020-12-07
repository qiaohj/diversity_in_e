library(raster)
library(data.table)
library(ggplot2)
rm(list=ls())
setDTthreads(threads=1)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Amphibians"
}
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-5
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
j=1
k=2

dispersals<-c(0:1)
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

mask<-readRDS("../../Objects/mask.rda")

for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:length(dispersals))){
    layer$M<-dispersals[k]
    target_folder<-sprintf("../../Objects/Diversity_%d/%s/%s_%d", threshold, group, layer$LABEL, layer$M)
    if (dir.exists(target_folder)){
      print(paste("Skip ", target_folder))
      next()
    }
    dir.create(target_folder, showWarnings = F, recursive = T)
    
    i=247
    diversity_df<-list()
    for (i in c(1:nrow(df_list))){
      item<-df_list[i,]
      item$sp<-gsub(" ", "_", item$sp)
      
      if (item$area<=0){
        next()
      }
      #print(paste(Sys.time(), 1))
      enm_folder<-sprintf("../../Objects/Niche_Models/%s/%s/dispersal_%d", group, item$sp, threshold)
      
      #print(sprintf("%s/%s_%s_%d.rda", enm_folder, layer$GCM, layer$SSP, layer$M))
      env_item_all<-readRDS(sprintf("%s/%s_%s_%d.rda", enm_folder, layer$GCM, layer$SSP, layer$M))
      env_item_all<-env_item_all[, c("x", "y", "mask_index", "YEAR")]
      
      source_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", source_folder))
      selected_cols<-c("x", "y", "mask_index")
      start_dis<-unique(start_dis[, ..selected_cols])
      start_dis$YEAR<-2020
      
      
      if (is.null(env_item_all)){
        env_item_all<-start_dis
      }else{
        env_item_all<-rbindlist(list(env_item_all, start_dis))
      }
      colnames(env_item_all)[4]<-"year"
      print(paste(group, layer$LABEL, layer$M, item$sp, i, nrow(df_list), "threshold", threshold))
      splited<-split(env_item_all, by="year")
      YYYY=2020
      for (YYYY in c(2020:2100)){
        diversity<-diversity_df[[as.character(YYYY)]]
        env_item<-splited[[as.character(YYYY)]]
        if (is.null(env_item)){
          next()
        }
        if (nrow(env_item)==0){
          next()
        }
        env_item$sp<-item$sp
        if (is.null(diversity)){
          diversity<-list(env_item)
        }else{
          diversity[[length(diversity)+1]]<-env_item
        }
        diversity_df[[as.character(YYYY)]]<-diversity
      }
    }
    saveRDS(diversity_df, sprintf("%s/diversity_df.rda", target_folder))
  }
}

if (F){
  tt<-"../../Objects/Diversity_1/Amphibians/EC-Earth3-Veg_SSP119_1"
  temp<-readRDS(sprintf("%s/diversity_df.rda", tt))
  
  temp_sub<-temp[["2020"]]
  temp_sub<-rbindlist(temp_sub)
  
  
  length(unique(temp_sub$sp))
}