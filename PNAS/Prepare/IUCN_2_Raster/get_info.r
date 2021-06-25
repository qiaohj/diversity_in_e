library(dplyr)
library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
groups<-c("Birds", "Mammals", "Amphibians", "Reptiles")
group<-"Birds"
for (group in groups){
  print(group)
  
  f<-list.files(sprintf("../../Objects_Full_species/IUCN_Distribution/%s", group))
  i<-f[2]
  result<-NULL
  for (i in f){
    
    sp<-gsub("\\.rda", "", i)
    df<-readRDS(sprintf("../../Objects_Full_species/IUCN_Distribution/%s/%s", group, i))
    area<-nrow(df)
    if (is.null(area)){
      area<--1
    }
    item<-data.frame(sp=sp, group=group, area=area)
    if (is.null(result)){
      result<-item
    }else{
      result<-bind_rows(result, item)
    }
  }
  if (F){
    result_bak<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
    table(result_bak[which(result_bak$area<=0), "area"])
    table(result[which(result_bak$area<=0), "area"])
    result[which(result$area<=0), ]
  }
  
  db<-readRDS(sprintf("../../Data/Raw/IUCN/%s/%s_df.rda", toupper(group), tolower(group)))
  area_db<-db%>%dplyr::select(binomial, SHAPE_Area)%>%dplyr::group_by(binomial)%>%dplyr::summarise(area_shape=sum(SHAPE_Area))
  
  result$sp<-gsub("_", " ", result$sp)
  result<-left_join(result, area_db, by=c("sp"="binomial"))
  
  saveRDS(result, sprintf("../../Objects/IUCN_List/%s.rda", group))
  
  #hist(result$area)
  #data.frame(table(result$area))
  #result[which(result$area==0),]
}
