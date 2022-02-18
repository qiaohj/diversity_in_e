library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Birds"
}

esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
print(paste("Reading", esm_ssp[esm_i]))
if (group=="Birds"){
  group_df<-readRDS("../../Data/Birds/bird_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  #group_disp$estimated_disp
  #group_disp2<-readRDS("../../Objects_PNAS/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  #dd<-merge(data.frame(group_disp[, c("iucn_name", "estimated_disp")]), 
  #          data.frame(group_disp2[, c("iucn_name", "estimated_disp")]), 
  #          by.x="iucn_name", by.y="iucn_name", all.x=T, all.y=F)
  group_full<-merge(group_df, group_disp, by.x="SCINAME", by.y="iucn_name", all=F)
  unique <- unique(group_full$SCINAME)
  unique<-as.character(unique)
  
}else{
  group_df<-readRDS("../../Data/Mammals/mammal_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  group_full<-merge(group_df, group_disp, by.x="binomial", by.y="Scientific", all=F)
  
  unique <- unique(group_full$binomial)
  unique<-as.character(unique)
  colnames(group_full)[1]<-"SCINAME"
  colnames(group_full)[27]<-"Shape_Area"
  colnames(group_disp)[1]<-"iucn_name"
  
}
group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]


bi<-group_full_sum_area[group_full_sum_area$sum_are<=1.5*min(group_full_sum_area$sum_are)]$SCINAME[1]
group_full_sum_area<-group_full_sum_area[sample(nrow(group_full_sum_area), nrow(group_full_sum_area))]
i=1
bi="Poicephalus rufiventris"
coms<-expand.grid(exposure_threshold=c(0, 5), dispersal=c(0, 1))
j=1
distribution_threshold<-100
dispersal_threshold<-50
#dispersal_threshold<-890
#distribution_threshold<-2000

distribution_all<-readRDS("../../Objects/N_cell_init_distribution.rda")
distribution_all<-distribution_all[N<=distribution_threshold]
group_full_sum_area<-group_full_sum_area[SCINAME  %in% group_full[estimated_disp<=dispersal_threshold]$SCINAME]
group_full_sum_area<-group_full_sum_area[SCINAME %in% distribution_all$sp]


for (i in 1:length(group_full_sum_area$SCINAME)) {
  start_time<-Sys.time()
  #print(paste("Start time:", start_time))
  for (j in 1:nrow(coms)){
    dispersal<-coms[j, "dispersal"]
    exposure_threshold<-coms[j, "exposure_threshold"]
    bi<-group_full_sum_area$SCINAME[i]
    bi_x<-group_full[SCINAME==bi]
    bi_y<-distribution_all[sp==bi]
    if (bi=="Aratinga maculata"){
      #asdf
    }
    
    
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
    fit_str<-sprintf("%s/fit.rda", target_folder)
    if (!file.exists(fit_str)){
      next()
    }
      
    
    
    for (item_str in esm_ssp){
      target<-sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, item_str,
                      exposure_threshold, dispersal)
      if (!file.exists(target)){
        next()
      }
      if (file.size(target)<1000){
        xx<-readRDS(target)
        if (length(xx)==0){
          print(paste(i, length(group_full_sum_area$SCINAME), bi, 
                      "exposure", exposure_threshold, "dispersal", dispersal,
                      "Distribution area", bi_x[1]$estimated_disp,
                      "Dispersal distance", bi_y[1]$N))
          print(paste("removing", target))
          unlink(target)
        }
      }
    }
  }
}

