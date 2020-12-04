library(raster)
library(dplyr)
#library(alphahull)
library(concaveman)
library(sf)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
if (is.na(group)){
  group<-"Amphibians"
}


source("commonFuns/functions.r")


df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=3369

final_df<-NULL
#beginCluster()
continent<-raster("../../Raster/Continent_ect4.tif")
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  
  if (item$area<=0){
    next()
  }
  
  target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
  target<-sprintf("%s/occ_with_env_and_continent.rda", target_folder)
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  print(paste(i, nrow(df_list), item$sp, target_folder))
  start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
  start_dis$continent<-raster::extract(continent, data.frame(start_dis[, c("x", "y")]))
  saveRDS(start_dis, target)
}