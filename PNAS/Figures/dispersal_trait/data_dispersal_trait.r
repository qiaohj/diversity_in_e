library(dplyr)
library(data.table)
library(raster)
rm(list=ls())
threshold<-1
ttt<-2
g<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
alt<-raster("../../Raster/ALT/alt_eck4.tif")
source("commonFuns/functions.r")
all_result<-NULL
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
i=2
for (threshold in c(5)){
  when_extinct<-NULL
  
  for (g in c("Mammals")){
    sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
    sp_list<-sp_list[which(sp_list$area>ttt),]
    
    sp_list$sp<-gsub(" ", "_", sp_list$sp)
    for (i in c(1:nrow(sp_list))){
      print(paste(threshold, g, i, nrow(sp_list)))
      for (SSP_i in SSPs){
        for (GCM_i in GCMs){
          
          for (da in c(0:1)){
            
            item<-sp_list[i,]
            item$GCM<-GCM_i
            item$SSP<-SSP_i
            item$dispersal<-da
            
            target_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", g, item$sp)
            dis<-readRDS(sprintf("%s/dispersal_%d/%s_%s_%d.rda",
                                 target_folder, threshold, item$GCM, item$SSP, item$dispersal))
            
            if (is.null(dis)){
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
            for (nn in names(item)){
              dis_all_se[, nn]<-item[, nn]
            }
            dis_all_se$threshold<-threshold
            all_result<-bind_dplyr(all_result, dis_all_se)
          }
        }
      }
    }
  }
}
saveRDS(all_result, sprintf("../../Objects_Full_species/dispersal_trait/dispersal_trait_exposure_%d_%s.rda", threshold, g))
