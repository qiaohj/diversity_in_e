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
  env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")
  exposure=0
  result<-NULL
  for (g in c("Birds", "Mammals")){
    for (exposure in c(0, 5)){
      df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
      sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
      sp_list$sp<-sp_list$SP
      sp_list$sp2<-gsub(" ", "_", sp_list$sp)
      df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
      when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
      i=1
      for (i in c(1:nrow(when_extinct))){
        print(paste(i, nrow(when_extinct), g, exposure))
        item<-when_extinct[i,]
        dis<-readRDS(sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d.rda",
                             item$group, item$sp, item$GCM, item$SSP, exposure, item$dispersal))
        if (is.null(dis)){
          min_distance<-NA
          TYPE<-"Unsuitable"
        }else{
          dis<-rbindlist(dis)
          if (nrow(dis)==0){
            min_distance<-NA
            TYPE<-"Unsuitable"
          }else{
            dis_prev<-dis[YEAR==(item$extinct_year-1)]
            fit<-readRDS(sprintf("../../Objects/Dispersal/%s/%s/fit.rda", 
                                   item$group, item$sp))
            
            env_item<-data.table(env_layers[[sprintf("%s_%s", item$GCM, item$SSP)]])
            env_item<-env_item[year==item$extinct_year]
            env_item<-env_item[between(bio1, fit$range_bio1_sd_min, fit$range_bio1_sd_max)&
                             between(bio5, fit$range_bio5_sd_min, fit$range_bio5_sd_max)&
                             between(bio6, fit$range_bio6_sd_min, fit$range_bio6_sd_max)&
                             between(bio12, fit$range_bio12_sd_min, fit$range_bio12_sd_max)&
                             between(bio13, fit$range_bio13_sd_min, fit$range_bio13_sd_max)&
                             between(bio14, fit$range_bio14_sd_min, fit$range_bio14_sd_max)]
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
        }
        item$N_Prev=nrow(dis_prev)
        item$N_Extinct=nrow(env_item)
        item$min_distance=min_distance
        item$type<-TYPE
        item$exposure<-exposure
        result<-bind_dplyr(result, item)
      }
    }
  }
  
  saveRDS(result, "../../Objects/Extinction_type/Extinction_type.rda")
}
df_sp_list<-list()

for (group in c("Birds", "Mammals")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list, fill=T)

df_sp_list$sp<-gsub(" ", "_", df_sp_list$SP)


result<-readRDS("../../Objects/Extinction_type/Extinction_type.rda")
#result<-result%>%dplyr::filter(sp %in% df_sp_list$sp)
result_count_y_year<-result%>%dplyr::group_by(type, dispersal, exposure, SSP, GCM, extinct_year)%>%
  dplyr::summarise(N=n())
result_count_y_year<-inner_join(result_count_y_year, result_count_extinct, by=c("dispersal", "exposure", "SSP", "GCM"))
result_count_y_year$proportion<-result_count_y_year$N/result_count_y_year$N_all_extinction
ggplot(result_count_y_year)+geom_point(aes(x=extinct_year, y=N, color=SSP, shape=GCM))+
  facet_grid(dispersal+type~exposure)


result_count<-result%>%dplyr::group_by(type, dispersal, exposure, SSP, GCM)%>%
  dplyr::summarise(N=n())

result_count_extinct<-result%>%dplyr::group_by(dispersal, exposure, SSP, GCM)%>%
  dplyr::summarise(N_all_extinction=n())

result_count<-inner_join(result_count, result_count_extinct, by=c("dispersal", "exposure", "SSP", "GCM"))
result_count$proportion<-result_count$N/result_count$N_all_extinction

ggplot(result_count)+geom_point(aes(x=SSP, y=N, color=GCM))+
  facet_grid(dispersal~exposure)
  
View(result_count%>%dplyr::filter((dispersal==1)&exposure==5))
result_count_se<-result_count%>%dplyr::group_by(type, dispersal, exposure, SSP)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   mean_proportion=mean(proportion),
                   sd_proportion=sd(proportion))
result_count_se[is.na(result_count_se)]<-0
result_count_se$dispersal<-ifelse(result_count_se$dispersal==0, "no dispersal", "with dispersal")
result_count_se$exposure<-ifelse(result_count_se$exposure==0, " no exposure", "5-year exposure")
result_count_se$label<-sprintf("%.1f%% ± %.1f%%", result_count_se$mean_proportion*100,
                               result_count_se$sd_proportion*100)
result_count_se$label2<-sprintf("%.1f ± %.1f", result_count_se$mean_N,
                               result_count_se$sd_N)

write.table(result_count_se, "../../Objects/Extinction_type/Extinction_type.csv", row.names = F, sep=",")
result_unreachable<-result%>%filter((type=="Unreachable")&(dispersal==1))%>%
  dplyr::group_by(SSP, dispersal)%>%
  dplyr::summarise(mean_distance=mean(min_distance),
                   sd_distance=sd(min_distance))

