library(raster)
#library(rgdal)
#library(rgeos)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

if (F){
  sp_dis_all_sub_N_all<-NULL
  sp_dis_extinct<-NULL
  sp_dis_all_se_all<-NULL
  g<-"Birds"
  print("Loading ENV DATA ...")
  env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")
  
  result<-NULL
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    for (threshold in c(1, 5)){
      df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
      sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
      #sp_list<-sp_list[which(sp_list$area>ttt),]
      sp_list$sp2<-gsub(" ", "_", sp_list$sp)
      df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
      when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
      for (i in c(1:nrow(when_extinct))){
        print(paste(i, nrow(when_extinct), g, threshold))
        item<-when_extinct[i,]
        dis<-readRDS(sprintf("../../Objects_Full_species/Niche_Models/%s/%s/dispersal_%d/%s_%s_%d.rda",
                             item$group, item$sp, threshold, item$GCM, item$SSP, item$dispersal))
        if (is.null(dis)){
          min_distance<-NA
          TYPE<-"Unsuitable"
        }else{
          dis_prev<-dis[YEAR==(item$extinct_year-1)]
          model<-readRDS(sprintf("../../Objects_Full_species/Niche_Models/%s/%s/fit.rda", 
                                 item$group, item$sp, threshold, item$GCM, item$SSP))
          
          env_item<-data.table(env_layers[[sprintf("%s_%s_%d", item$GCM, item$SSP, item$extinct_year)]])
          
          env_item<-env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                                (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                                (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]
          if (nrow(env_item)>0){
            min_distance<-env_item[, min_dist(x, y, dis_prev), by = 1:nrow(env_item)]
            min_distance<-min(min_distance$V1)/1000
            if (F){
              ggplot()+
                geom_tile(data=env_item, aes(x=x, y=y))+
                #geom_tile(data=dis_prev, aes(x=x, y=y), fill="red")
                geom_point(data=dis_prev, aes(x=x, y=y), color="red")
              
            }
            TYPE<-"Unreachable"
          }else{
            min_distance<-NA
            TYPE<-"Unsuitable"
          }
        }
        item$N_Prev=nrow(dis_prev)
        item$N_Extinct=nrow(env_item)
        item$min_distance=min_distance
        item$type<-TYPE
        item$threshold<-threshold
        result<-bind_dplyr(result, item)
      }
    }
  }
  
  saveRDS(result, "../../Objects_Full_species/Extinction_type/Extinction_type.rda")
}
df_sp_list<-list()

for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list)

ttt<-2
df_sp_list<-df_sp_list[area>ttt]
df_sp_list$sp<-gsub(" ", "_", df_sp_list$sp)


result<-readRDS("../../Objects_Full_species/Extinction_type/Extinction_type.rda")
result<-result%>%dplyr::filter(sp %in% df_sp_list$sp)
result_count<-result%>%dplyr::group_by(type, dispersal, threshold, SSP, GCM)%>%
  dplyr::summarise(N=n())
result_count_se<-result_count%>%dplyr::group_by(type, dispersal, threshold, SSP)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N))
result_count_se$dispersal<-ifelse(result_count_se$dispersal==0, "no dispersal", "with dispersal")
result_count_se$exposure<-ifelse(result_count_se$threshold==1, " no exposure", "5-year exposure")
result_count_se<-result_count_se%>%ungroup()%>%dplyr::select(-c(threshold))
write.table(result_count_se, "../../Objects_Full_species/Extinction_type/Extinction_type.csv", row.names = F, sep=",")
result_unreachable<-result%>%filter((type=="Unreachable")&(dispersal==1))%>%
  dplyr::group_by(SSP, dispersal)%>%
  dplyr::summarise(mean_distance=mean(min_distance),
                   sd_distance=sd(min_distance))

