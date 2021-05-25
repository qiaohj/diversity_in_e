library(rgdal)
library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
shape<-readOGR(dsn="../../Data/Raw/IUCN/MAMMALS", layer="MAMMALS_TERRESTRIAL_ONLY")


sp_list<-read.table("../../Primates.log", head=F, sep=",", stringsAsFactors = F)
sp<-sp_list$V1[1]
sp_list$area<-NA
for (i in c(1:nrow(sp_list))){
  sp<-sp_list[i, "V1"]
  shp<-shape[(shape$binomial==sp),]
  if (nrow(shp)==0){
    asdf
    next()
  }
  area<-sum(area(shp))/1e6
  sp_list[i, "area"]<-area
}
write.table(sp_list, "../../Primates.csv", row.names=F, sep=",")


sub_shp<-shape[(shape$order_=="PRIMATES"),]
PRIMATES_list<-unique(sub_shp$binomial)
write.table(PRIMATES_list, "../../PRIMATES_list.csv", row.names=F, sep=",")
