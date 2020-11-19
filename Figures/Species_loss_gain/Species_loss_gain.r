library(raster)
#library(rgdal)
#library(rgeos)
#library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/functions.r")
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Amphibians"
}


GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
j=1
k=1
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
dispersals<-c(0:2)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
add_location<-function(indices, location, type){
  location$metric<-indices
  location$type<-type
  location
}
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}
threshold<-5
for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:length(dispersals))){
    layer$M<-dispersals[k]
    target_folder<-sprintf("../../Objects/Diversity_%d/%s/%s_%d", threshold, group, layer$LABEL, layer$M)
    target<-sprintf("%s/loss_gain.rda", target_folder)
    if (file.exists(target)){
      next()
    }
    saveRDS(NULL, target)
    print(paste("READING DATA", target_folder))
    diversity_df<-readRDS(sprintf("%s/diversity_df.rda", target_folder))
    YYYY<-names(diversity_df)[1]
    loss_gain_df<-NULL
    
    for (YYYY in names(diversity_df)){
      
      diversity<-diversity_df[[YYYY]]
      print(paste("BINDING DATA", YYYY, target_folder))
      diversity<-rbindlist(diversity)
      diversity<-diversity%>%distinct()
      
      print(paste("EXTRACTING DATA", YYYY, target_folder))
      if (F){
        system.time({
          diversity$mask_index<-raster::extract(mask, diversity[, c("x", "y")])  
        })
      }
      #system.time({
      #  diversity<-inner_join(diversity, points, by=c("x", "y", "mask_index"))
      #})
      print(paste("SPLITING DATA", YYYY, target_folder))
      diversity<-diversity[complete.cases(diversity),]
      diversity$mask_index_str<-paste("i", diversity$mask_index, sep="_")
      diversity_split<-diversity %>% 
        named_group_split(mask_index_str)
      print(paste("CALCULATING LOSS/GAIN DATA", YYYY, target_folder))
      if (YYYY=="2020"){
        base_diversity<-diversity_split
      }else{
        id<-names(diversity_split)[1000]
        for (id in names(diversity_split)){
          if (id %in% names(base_diversity)){
            base_sp<-base_diversity[[id]]$sp
          }else{
            base_sp<-c()
          }
          
          new_sp<-diversity_split[[id]]$sp
          intersect_sp<-intersect(base_sp, new_sp)
          n_overlap<-length(intersect_sp)
          loss<-length(base_sp)-n_overlap
          gain<-length(new_sp)-n_overlap
          item<-data.frame(YEAR=YYYY, 
                           mask_index=diversity_split[[id]][1, "mask_index"],
                           x=diversity_split[[id]][1, "x"],
                           y=diversity_split[[id]][1, "y"],
                           n_overlap=n_overlap,
                           n_loss=loss,
                           n_gain=gain
                           )
          loss_gain_df<-bind(loss_gain_df, item)
        }
      }
      
    }
    saveRDS(loss_gain_df, target)

    
  }
}
