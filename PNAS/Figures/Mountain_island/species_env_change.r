library(sf)
library(raster)
library(data.table)
library(ggplot2)
library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

envs<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")

birds<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")

sp<-birds$SP[1]
i=1
all_info<-list()
sps<-unique(birds$SP)
for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/Dispersal/Birds/%s/initial_disp_exposure_0_dispersal_0.rda", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  if (nrow(dis)==0){
    next()
  }
  com<-names(envs)[1]
  item_set<-NULL
  for (com in names(envs)){
    env_items<-envs[[com]][mask_100km %in% dis$mask_100km]
    v<-"bio1"
    for (v in c("bio1", "bio5", "bio6",  "bio12", "bio13", "bio14")){
      range<-max(env_items[, ..v], na.rm = T) - min(env_items[, ..v], na.rm = T)
      sd<-sd(pull(env_items[, ..v]), na.rm = T)
      item<-data.frame(sp=sp, range=range, sd=sd, var=v, 
                       GCM=env_items[1, ]$GCM, SSP=env_items[1, ]$SSP, group="Bird")
      if (is.null(item_set)){
        item_set<-item
      }else{
        item_set<-rbind(item_set, item)
      }
    }
  }
  all_info[[sp]]<-item_set
}


mammals<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")
sp<-mammals$SP[1]
i=1
sps<-unique(mammals$SP)
for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/Dispersal/Mammals/%s/initial_disp_exposure_0_dispersal_0.rda", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  if (nrow(dis)==0){
    next()
  }
  com<-names(envs)[1]
  item_set<-NULL
  for (com in names(envs)){
    env_items<-envs[[com]][mask_100km %in% dis$mask_100km]
    v<-"bio1"
    for (v in c("bio1", "bio5", "bio6",  "bio12", "bio13", "bio14")){
      range<-max(env_items[, ..v], na.rm = T) - min(env_items[, ..v], na.rm = T)
      sd<-sd(pull(env_items[, ..v]), na.rm = T)
      item<-data.frame(sp=sp, range=range, sd=sd, var=v, 
                       GCM=env_items[1, ]$GCM, SSP=env_items[1, ]$SSP, group="Mammal")
      if (is.null(item_set)){
        item_set<-item
      }else{
        item_set<-rbind(item_set, item)
      }
    }
  }
  all_info[[sp]]<-item_set
}

all_info<-rbindlist(all_info)

saveRDS(all_info, "../../Objects/Island/species_env_change.rda")

