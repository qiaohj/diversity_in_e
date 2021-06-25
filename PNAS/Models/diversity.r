library(raster)
library(data.table)
library(ggplot2)
rm(list=ls())
setDTthreads(threads=1)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Mammals"
}
exposure<-as.numeric(args[2])
if (is.na(exposure)){
  exposure<-0
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

#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

mask<-readRDS("../../Objects/mask.rda")

for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  
    target_folder<-sprintf("../../Objects/Diversity_exposure_%d/%s/%s", exposure, group, layer$LABEL)
    if (dir.exists(target_folder)){
      print(paste("Skip ", target_folder))
      next()
    }
    dir.create(target_folder, showWarnings = F, recursive = T)
    
    i=1
    diversity_df<-list()
    for (i in c(1:nrow(df_list))){
      print(paste(j, nrow(layer_df), i, nrow(df_list)))
      item<-df_list[i,]
      item$SP<-gsub(" ", "_", item$SP)
      
      if (item$N_CELL<=0){
        next()
      }
      #print(paste(Sys.time(), 1))
      source_folder<-sprintf("../../Objects/%s/%s", group, item$SP)
      start_dis<-readRDS(sprintf("%s/initial_disp_exposure_%d.rda", source_folder, exposure))
      #print(sprintf("%s/%s_%s_%d.rda", enm_folder, layer$GCM, layer$SSP, layer$M))
      all_dis<-readRDS(sprintf("%s/%s_%s_%d.rda", source_folder, layer$GCM, layer$SSP, exposure))
      all_dis[["2020"]]<-start_dis
      
      YYYY=2024
      for (YYYY in c(2020:2100)){
        #print(YYYY)
        diversity<-diversity_df[[as.character(YYYY)]]
        env_item<-all_dis[[as.character(YYYY)]]
        
        if (is.null(env_item)){
          next()
        }
        if (nrow(env_item)==0){
          next()
        }
        selected_cols<-c("x", "y", "mask_100km")
        env_item<-env_item[, ..selected_cols]
        env_item$YEAR<-YYYY
        env_item$sp<-item$SP
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

if (F){
  tt<-"../../Objects/Diversity_1/Amphibians/EC-Earth3-Veg_SSP119_1"
  temp<-readRDS(sprintf("%s/diversity_df.rda", tt))
  
  temp_sub<-temp[["2020"]]
  temp_sub<-rbindlist(temp_sub)
  
  
  length(unique(temp_sub$sp))
}