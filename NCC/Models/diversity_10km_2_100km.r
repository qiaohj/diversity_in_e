library(raster)
library(data.table)
library(ggplot2)
#rm(list=ls())
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

dispersal<-as.numeric(args[3])
if (is.na(dispersal)){
  dispersal<-0
}

if (group=="Birds"){
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  
}else{
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  colnames(group_disp)[1]<-"iucn_name"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
points_10km<-readRDS("../../Raster/points_10km.rda")
colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-data.table(rasterToPoints(mask_100km))
colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
points_10km<-merge(points_10km, points_100km, by="mask_100km")
i=1
j=1
k=2

#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]

#mask<-readRDS("../../Objects/mask.rda")
j=1
for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  target_folder<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d_10km_2_100km/%s/%s", exposure, dispersal, group, layer$LABEL)
  
  
  if (dir.exists(target_folder)){
    print(paste("Skip ", target_folder))
    next()
  }
  dir.create(target_folder, showWarnings = F, recursive = T)
  
  i=1
  diversity_df<-list()
  item<-df_list[SP=="Zosterops_olivaceus"]
  for (i in c(1:nrow(df_list))){
    print(paste("exposure:", exposure, "dispersal:", dispersal, j, nrow(layer_df), i, nrow(df_list)))
    item<-df_list[i,]
    item$SP<-gsub(" ", "_", item$SP)
    
    if (item$N_CELL<=0){
      next()
    }
    #print(paste(Sys.time(), 1))
    source_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, item$SP)
    start_dis_str<-sprintf("%s/initial_disp_10km_exposure_%d_dispersal_%d.rda", source_folder, exposure, dispersal)
    item_str<-sprintf("%s/%s_%s_%d_dispersal_%d_10km.rda", source_folder, layer$GCM, layer$SSP, exposure, dispersal)
    res<-"10km"
    if ((!file.exists(start_dis_str))|(!file.exists(item_str))){
      start_dis_str<-sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda", source_folder, exposure, dispersal)
      item_str<-sprintf("%s/%s_%s_%d_dispersal_%d.rda", source_folder, layer$GCM, layer$SSP, exposure, dispersal)
      res<-"100km"
      if (!file.exists(start_dis_str)){
        next()
      }
    }
    
    start_dis<-readRDS(start_dis_str)
    colnames(start_dis)[3]<-"mask"
    #print(sprintf("%s/%s_%s_%d.rda", enm_folder, layer$GCM, layer$SSP, layer$M))
    
    
    all_dis<-readRDS(item_str)
    
    all_dis[["2020"]]<-start_dis
    
    YYYY=2021
    for (YYYY in c(2020:2100)){
      #print(YYYY)
      diversity<-diversity_df[[as.character(YYYY)]]
      env_item<-all_dis[[as.character(YYYY)]]
      if (is.null(env_item)){
        next()
      }
      colnames(env_item)[3]<-"mask"
      
      #if (YYYY!=2020){
      #env_item<-env_item[suitable==1]
      #}
      if (nrow(env_item)==0){
        next()
      }
      selected_cols<-c("x", "y", "mask")
      env_item<-env_item[, ..selected_cols]
      env_item$YEAR<-YYYY
      env_item$sp<-item$SP
      env_item$res<-res
      if (res=="10km"){
        env_item<-merge(env_item, points_10km, by.x="mask", by.y="mask_10km")
        
        env_item<-env_item[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR", "sp", "res")]
        
      }else{
        colnames(env_item)[c(1:3)]<-c("x_100km", "y_100km", "mask_100km")
        env_item$N<-1
      }
      if (is.null(diversity)){
        diversity<-list(env_item)
        names(diversity)<-item$SP
      }else{
        diversity[[item$SP]]<-env_item
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