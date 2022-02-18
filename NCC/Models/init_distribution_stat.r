library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
library(rmapshaper)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Mammals"
}
all<-list()
for (group in c("Birds", "Mammals")){
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
  PRESENCE<-c(1,2,3,4,5)
  ORIGIN<-c(1,2,3,5,6)
  SEASONAL<-c(1,2)
  print("reading env layers")
  i=1
  group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
  group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]
  
  
  bi<-group_full_sum_area[group_full_sum_area$sum_are<=1.5*min(group_full_sum_area$sum_are)]$SCINAME[1]
  group_full_sum_area<-group_full_sum_area[sample(nrow(group_full_sum_area), nrow(group_full_sum_area))]
  exposure_threshold=0
  dispersal=0
  
  for (i in 1:length(group_full_sum_area$SCINAME)) {
    bi<-group_full_sum_area$SCINAME[i]
    if (bi=="Aratinga maculata"){
      #asdf
    }
    print(paste(i, length(group_full_sum_area$SCINAME), bi, 
                "exposure", exposure_threshold, "dispersal", dispersal))
    
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
    check_point<- sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda", 
                          target_folder, exposure_threshold, dispersal)
    if (!file.exists(check_point)){
      next()
    }
    item<-readRDS(check_point)
    df_item<-data.frame(sp=bi, group=group, N=nrow(item))
    all[[paste(group, bi)]]<-df_item
  }
}

all<-rbindlist(all)
saveRDS(all, "../../Objects/N_cell_init_distribution.rda")
all<-readRDS("../../Objects/N_cell_init_distribution.rda")
quantile(all$N, c(seq(0, 1, 0.1), 0.99))
if (F){
  group_disp_bird<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  group_disp_mammal<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  estimated_disp<-c(group_disp_bird$estimated_disp, group_disp_mammal$estimated_disp)
  hist(estimated_disp)
  quantile(estimated_disp, c(seq(0, 1, 0.1), 0.99))
}
