library(dplyr)
library(data.table)
library(raster)
rm(list=ls())
exposure<-5
g<-"Mammals"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
alt<-raster("../../Raster/ALT/alt_eck4.tif")
source("commonFuns/functions.r")
all_result<-NULL
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
i=2

points_10km<-readRDS("../../Raster/points_10km.rda")
colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-data.table(rasterToPoints(mask_100km))
colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
points_10km<-merge(points_10km, points_100km, by="mask_100km")
i=16
SSP_i<-SSPs[1]
GCM_i<-GCMs[1]
da<-0
for (exposure in c(0)){
  when_extinct<-NULL
  
  for (g in c("Mammals")){
    sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
    
    sp_list$sp<-gsub(" ", "_", sp_list$SP)
    for (i in c(1:nrow(sp_list))){
      print(paste(exposure, g, i, nrow(sp_list)))
      for (SSP_i in SSPs){
        for (GCM_i in GCMs){
          
          for (da in c(0:1)){
            
            item<-sp_list[i,]
            item$GCM<-GCM_i
            item$SSP<-SSP_i
            item$dispersal<-da
            target_folder<-sprintf("../../Objects/Dispersal/%s/%s", g, item$sp)
            f<-sprintf("%s/%s_%s_%d_dispersal_%d_10km.rda",
                       target_folder, item$GCM, item$SSP, exposure, item$dispersal)
            res<-"10km"
            if (!file.exists(f)){
              f<-sprintf("%s/%s_%s_%d_dispersal_%d.rda",
                         target_folder, item$GCM, item$SSP, exposure, item$dispersal)
              res<-"100km"
              
            }
            if (!file.exists(f)){
              
              next()
            }
            dis<-readRDS(f)
            
            if (is.null(dis)){
              next()
            }
            if (length(dis)==0){
              next()
            }
            dis<-rbindlist(dis)
            if (nrow(dis)==0){
              next()
            }
            if (res=="10km"){
              dis<-merge(dis, points_10km, by.x="mask_10km", by.y="mask_10km")
              dis<-dis[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR")]
              colnames(dis)[1:3]<-c("x", "y", "mask_100km")
            }
            if (nrow(dis)==0){
              next()
            }
            item$extinct_year<-max(dis$YEAR)+1
            dis_start<-dis%>%dplyr::filter(YEAR==2021)
            end_year<-item$extinct_year
            if (is.infinite(end_year)){
              end_year<-2101
            }
            dis_end<-dis%>%dplyr::filter(YEAR==(end_year-1))
            
            dis_start$type<-"start"
            dis_end$type<-"end"
            dis_all<-bind_rows(dis_start, dis_end)
            dis_all$alt<-raster::extract(alt, dis_all[, c("x", "y")])
            dis_all<-dis_all%>%dplyr::filter(!is.na(alt))
            dis_all_se<-dis_all%>%dplyr::group_by(type)%>%
              dplyr::summarise(mean_alt=mean(alt),
                               median_alt=quantile(alt, 0.5),
                               mean_y=mean(y),
                               mean_abs_y=mean(abs(y)))
            if (nrow(dis_all_se)==1){
              next()
            }
            if (F){
              plot(alt)
              points(dis_all$x, dis_all$y, col=factor(dis_all$type))
            }
            item<-data.frame(item)
            for (nn in names(item)){
              dis_all_se[, nn]<-item[, nn]
            }
            dis_all_se$exposure<-exposure
            all_result<-bind_dplyr(all_result, dis_all_se)
          }
        }
      }
    }
  }
}
saveRDS(all_result, sprintf("../../Objects/dispersal_trait_10km/dispersal_trait_exposure_%d_%s.rda", exposure, g))
