library(raster)
#library(rgdal)
library(MASS)
#library(cluster)
library(dplyr)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("functions.r")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
Labels<-expand.grid(GCM=GCMs, SSP=SSPs)


df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
exposure_year<-5
year_range<-c(2015:2100)
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folders<-c(sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/exposure", target_folder)
    source_rda<-sprintf("%s/exposure_end.rda", target)
    target_rda<-sprintf("%s/exposure_acc.rda", target)
    if (file.exists(target_rda)){
      if (file.size(target_rda)>50){
        next()
      }
    }    
    saveRDS(NULL, target_rda)
    df<-readRDS(source_rda)
    
    print(paste(group, i, nrow(df_list), item$sp, target, nrow(df)))
    min_years<-df%>%dplyr::filter(is_exposure==1)%>%dplyr::ungroup()%>%dplyr::group_by(mask_index, GCM, SSP)%>%
      dplyr::summarise(min_year=min(year))
    j=2
    df<-data.table(df)
    setindex(df, mask_index, GCM, SSP, year)
    for (j in c(1:nrow(min_years))){
      #print(paste(j, nrow(min_years)))
      index_item<-min_years[j,]
      if (F){
        test<-df%>%dplyr::filter((mask_index==index_item$mask_index)&
                           (GCM==index_item$GCM)&
                           (SSP==index_item$SSP)&
                           (year>=index_item$min_year))
        test_new<-new_df%>%dplyr::filter((mask_index==index_item$mask_index)&
                                   (GCM==index_item$GCM)&
                                   (SSP==index_item$SSP)&
                                   (year>=index_item$min_year))
        test_new2<-new_df%>%dplyr::filter((mask_index==index_item$mask_index)&
                                           (GCM==index_item$GCM)&
                                           (SSP==index_item$SSP)&
                                           (year<=index_item$min_year))
        
        system.time({
          df<-df %>%
            mutate(is_exposure=replace(is_exposure, 
                                       (mask_index==index_item$mask_index)&
                                         (GCM==index_item$GCM)&
                                         (SSP==index_item$SSP)&
                                         (year>=index_item$min_year), 
                                       1))
        })
        
      }
      
      #system.time({
        df[((mask_index==index_item$mask_index)&
             (GCM==index_item$GCM)&
             (SSP==index_item$SSP)&
             (year %in% index_item$min_year:2100)), "is_exposure"]<-1
      #})
    }
    saveRDS(df, target_rda)
  }
}


if (F){
  df<-readRDS(target_rda)
  test<-df[2,]
  t_mask_index<-test$mask_index
  t_GCM<-test$GCM
  t_SSP<-test$SSP
  t_year<-test$year
  t_is_out<-test$is_out
}

