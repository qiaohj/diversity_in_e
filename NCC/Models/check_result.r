library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
library(rmapshaper)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

bird_df<-readRDS("../../Data/Birds/bird_df.rda")
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
bird_full<-merge(bird_df, bird_disp, by.x="SCINAME", by.y="iucn_name", all=F)

unique <- unique(bird_full$SCINAME)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
print("reading env layers")

mask_100km<-raster("../../Raster/mask_100km.tif")

bi<-bird_full_sum_area[bird_full_sum_area$sum_are<=1.5*min(bird_full_sum_area$sum_are)]$SCINAME[1]
bird_full_sum_area<-bird_full_sum_area[sample(nrow(bird_full_sum_area), nrow(bird_full_sum_area))]
i=1
bi="Amazona xantholora"
for (i in 1:length(bird_full_sum_area$SCINAME)) {
  bi<-bird_full_sum_area$SCINAME[i]
  print(paste(i, length(bird_full_sum_area$SCINAME), bi))
  target_folder<-sprintf("../../Objects/Birds/%s", gsub(" ", "_", bi))
  item_str<-"EC-Earth3-Veg_SSP245"
  exposure_threshold<-5
  dispersal<-0
  if (!file.exists(sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                          exposure_threshold, dispersal))){
    next()
  }
  dispersal_log_0<-readRDS(sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                                 exposure_threshold, dispersal))
  dispersal_log_0<-rbindlist(dispersal_log_0)
  
  dispersal_log_0_s<-dispersal_log_0[,.(N=.N), by=list(YEAR)]
  dispersal_log_1<-readRDS(sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                                   exposure_threshold, 1))
  dispersal_log_1<-rbindlist(dispersal_log_1)
  
  dispersal_log_1_s<-dispersal_log_1[,.(N=.N), by=list(YEAR)]
  
  xx<-merge(dispersal_log_0_s, dispersal_log_1_s, by="YEAR", all=T)
  if (nrow(xx[N.x>N.y])!=0){
    asdf
  }
}
