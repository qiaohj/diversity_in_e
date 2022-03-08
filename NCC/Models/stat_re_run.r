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
sp_list<-readRDS("../../Objects/extincted_species_list.rda")
sp_list<-sp_list[N>200|estimated_disp>100]
cols<-c("Scientific", "group.x")
table(unique(sp_list[, ..cols])$group.x)




group<-"Birds"
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


unique<-unique[sample(length(unique), length(unique))]

group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]


bi<-group_full_sum_area[group_full_sum_area$sum_are<=1.5*min(group_full_sum_area$sum_are)]$SCINAME[1]
group_full_sum_area<-group_full_sum_area[sample(nrow(group_full_sum_area), nrow(group_full_sum_area))]
i=1
bi="Poicephalus rufiventris"
coms<-expand.grid(exposure_threshold=c(0, 5), dispersal=c(0, 1))
j=1
distribution_threshold<-200
dispersal_threshold<-100
#dispersal_threshold<-890
#distribution_threshold<-2000

distribution_all<-readRDS("../../Objects/N_cell_init_distribution.rda")
group_full_sum_area<-merge(group_full_sum_area, distribution_all, by.x="SCINAME", by.y="sp")
group_full_sum_area<-merge(group_full_sum_area, group_disp, by.x="SCINAME", by.y="iucn_name")


final_df<-list()
for (i in 1:length(group_full_sum_area$SCINAME)) {
  item<-group_full_sum_area[i]
  bi<-group_full_sum_area$SCINAME[i]
  target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
  fit_str<-sprintf("%s/fit_seasonal_2.rda", target_folder)
  if (file.exists(fit_str)){
    
    final_df[[as.character(i)]]<-item
  }
}
final_df<-rbindlist(final_df)

dim(final_df[N<=200&estimated_disp<=100])
dim(final_df[N>200|estimated_disp>100])

