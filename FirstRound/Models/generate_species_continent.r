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

for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  
  if (item$area<=0){
    next()
  }
  print(paste(i, nrow(df_list), item$sp))
  target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
  target<-sprintf("%s/occ_with_env_and_continent.rda", target_folder)
  df<-readRDS(target)
  sp<-data.frame(sp=item$sp, group=group, continent=unique(df$continent))
  final_df<-bind_dplyr(final_df, sp)
}
saveRDS(final_df, sprintf("../../Objects/SP_Continent/%s.rda", group))
