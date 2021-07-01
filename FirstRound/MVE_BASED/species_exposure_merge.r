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
  group<-"Mammals"
}
source("functions.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
Labels<-expand.grid(GCM=GCMs, SSP=SSPs)


df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
exposure_year<-5
year_range<-c(2015:2100)
final_df<-NULL
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
    print(paste(group, i, nrow(df_list), item$sp, target, nrow(df)))
    df<-readRDS(target_rda)
    df_sub<-df[,c("mask_index", "year", "GCM", "SSP", "range_type", "is_exposure")]
    df_sub$N_SP<-1
    df_sub$N_exposure<-df_sub$is_exposure
    df_sub$range_type_0<-0
    df_sub$range_type_1<-0
    df_sub$range_type_2<-0
    df_sub$range_type_3<-0
    df_sub$range_type_4<-0
    df_sub[which(df_sub$range_type==0), "range_type_0"]<-1
    df_sub[which((df_sub$range_type==1)&(df_sub$is_exposure==1)), "range_type_1"]<-1
    df_sub[which((df_sub$range_type==2)&(df_sub$is_exposure==1)), "range_type_2"]<-1
    df_sub[which((df_sub$range_type==3)&(df_sub$is_exposure==1)), "range_type_3"]<-1
    df_sub[which((df_sub$range_type==4)&(df_sub$is_exposure==1)), "range_type_4"]<-1
    
    
    if (is.null(final_df)){
      final_df<-df_sub
    }else{
      final_df<-full_join(final_df, df_sub, by=c("mask_index", "year", "GCM", "SSP"))
      final_df[is.na(final_df)]<-0
      final_df$N_SP.x<-final_df$N_SP.x+final_df$N_SP.y
      final_df$N_exposure.x<-final_df$N_exposure.x+final_df$N_exposure.y
      final_df$range_type_0.x<-final_df$range_type_0.x+final_df$range_type_0.y
      final_df$range_type_1.x<-final_df$range_type_1.x+final_df$range_type_1.y
      final_df$range_type_2.x<-final_df$range_type_2.x+final_df$range_type_2.y
      final_df$range_type_3.x<-final_df$range_type_3.x+final_df$range_type_3.y
      final_df$range_type_4.x<-final_df$range_type_4.x+final_df$range_type_4.y
      final_df<-final_df[, c("mask_index", "year", "GCM", "SSP", "N_SP.x", "N_exposure.x",
                             "range_type_0.x", "range_type_1.x", "range_type_2.x",
                             "range_type_3.x", "range_type_4.x")]
      colnames(final_df)<-c("mask_index", "year", "GCM", "SSP", "N_SP", "N_exposure",
                            "range_type_0", "range_type_1", "range_type_2",
                            "range_type_3", "range_type_4")
    }
  }
}

saveRDS(final_df, sprintf("../../Objects/Species_exposure/%s.rda", group))
if (F){
  test<-df[2,]
  t_mask_index<-test$mask_index
  t_GCM<-test$GCM
  t_SSP<-test$SSP
  t_year<-test$year
  t_is_out<-test$is_out
}

