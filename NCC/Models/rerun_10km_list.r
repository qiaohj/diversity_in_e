library(data.table)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
#rm(list=ls())
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

env<-NULL
print(group)

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
coms<-expand.grid(e_str=esm_ssp, exposure=c(0, 5), dispersal=c(0, 1))
group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]
rerun_list<-list()
for (i in 1:length(group_full_sum_area$SCINAME)) {
  print(paste(i, length(group_full_sum_area$SCINAME)))
  bi<-group_full_sum_area$SCINAME[i]
  bi<-gsub(" ", "_", bi)
  source_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, bi)
  e_item<-coms[1,]
  e_str<-1
  
  for (e_str in c(1:nrow(coms))){
    e_item<-coms[e_str,]
    item_str<-sprintf("%s/%s_%d_dispersal_%d_10km.rda", source_folder, e_item$e_str, e_item$exposure, e_item$dispersal)
    e_item$run_10km<-file.exists(item_str)
    item_str<-sprintf("%s/%s_%d_dispersal_%d.rda", source_folder, e_item$e_str, e_item$exposure, e_item$dispersal)
    e_item$run_100km<-file.exists(item_str)
    e_item$sp<-bi
    e_item$group<-group
    rerun_list[[length(rerun_list)+1]]<-e_item
  }
  
}

rerun_list<-rbindlist(rerun_list)
saveRDS(rerun_list, sprintf("../../Objects/rerun_10km_%s.rda", group))

View(rerun_list[, .(N=.N), by=list(e_str, exposure, dispersal, run_10km, run_100km)])

