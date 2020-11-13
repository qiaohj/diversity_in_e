library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
source("functions.r")
fun_exposure<-function(t_mask_index, t_year, t_GCM, t_SSP, t_is_out, df, exposure_year){
  if (t_is_out==0){
    return(0)
  }else{
    sub_item<-df%>%dplyr::filter((mask_index==t_mask_index)&
                                   (GCM==t_GCM)&
                                   (SSP==t_SSP)&
                                   (between(year, t_year, t_year+exposure_year-1)))
    sum_out<-sum(sub_item$is_out)
    if (sum_out==nrow(sub_item)){
      return(1)
    }else{
      return(0)
    }
  }
}

fun_exposure_sub<-function(t_mask_index, t_year, t_is_out, df_sub, exposure_year){
  if (t_is_out==0){
    return(0)
  }else{
    sub_item<-df_sub%>%dplyr::filter((mask_index==t_mask_index)&
                                   (between(year, t_year, t_year+exposure_year-1)))
    sum_out<-sum(sub_item$is_out)
    if (sum_out==nrow(sub_item)){
      return(1)
    }else{
      return(0)
    }
  }
}

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
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
    target_rda<-sprintf("%s/exposure_end.rda", target)
    if (file.exists(target_rda)){
      next()
    }
    saveRDS(NULL, target_rda)
    
    df<-readRDS(sprintf("%s/exposure.rda", target))
    print(paste(group, i, nrow(df_list), item$sp, target, nrow(df)))
    #ll<-df%>%dplyr::distinct(x, y)
    #plot(ll$x, ll$y)
    #0: in range 1: out of mve but in range box 2: out because of temp 3: out because of prec 4: out because of temp and prec
    #table(df$range_type)
    df$is_out<-1
    df[which(df$range_type==0), "is_out"]<-0
    print(system.time({
    final_df<-NULL
    for (label_i in c(1:nrow(Labels))){
      label<-Labels[label_i,]
      df_item<-df%>%dplyr::filter((GCM==label$GCM)&(SSP==label$SSP))
      df_item<-df_item %>%dplyr::rowwise()%>% 
        mutate(is_exposure=fun_exposure_sub(mask_index, year, is_out, df_item, exposure_year))
      final_df<-bind(final_df, df_item)
    }
    }))
    
    saveRDS(final_df, target_rda)
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

