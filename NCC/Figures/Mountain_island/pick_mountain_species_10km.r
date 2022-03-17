library(rgdal)
library(raster)
library(data.table)
library(ggplot2)
library(sf)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_10km.tif")

if (F){
  #Detect the islands
  points<-readRDS("../../Objects/Island/islands.rda")
  mountain<-readOGR("../../Shape/GMBA mountain inventory V1.2/GMBA Mountain Inventory_v1.2-World.shp")
  mountain<-spTransform(mountain, crs(mask))
  
  mountain_blank<-mask
  values(mountain_blank)[!is.na(values(mountain_blank))]<-0
  df_index<-list()
  for (i in c(1:nrow(mountain))){
    print(i)
    item<-mountain[i,]
    mountain_item<-mask(mask, item)
    index<-!is.na(values(mountain_item))
    mountain_blank[index]<-i
    item_df<-data.frame(Name=item$Name, Country=item$Country, index=i)
    df_index[[length(df_index)+1]]<-item_df
  }
  #mountains<-crop(mask, extent(mountain))
  mountains<-mask(mask, mountain)
  writeRaster(mountains, "../../Objects/Island/mountain_10km.tif", overwrite=T)
  
  mountain<-raster("../../Objects/Island/mountain_10km.tif")
  plot(mountain)
  index<-!is.na(values(mountain))
  mountain<-mask
  values(mountain)[!is.na(values(mountain))]<-0
  plot(mountain)
  values(mountain)[index]<-1
  plot(mountain)
  writeRaster(mountain, "../../Objects/Island/mountain_10km.tif", overwrite=T)
}

points<-readRDS("../../Objects/Island/islands.rda")
birds<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")
mountain<-raster("../../Objects/Island/mountain_10km.tif")
sp<-birds$SP[1]
all_info<-list()
sps<-unique(birds$SP)
i=1548
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
extant_mask<-raster::extent(mountain)
sp<-"Acropternis_orthonyx"
for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", "Birds", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  dis<-dis[which((dis$PRESENCE %in% PRESENCE)&
                   (dis$ORIGIN %in% ORIGIN)&
                   (dis$SEASONAL %in% SEASONAL)),]
  if (nrow(dis)==0){
    next()
  }
  extend<-st_bbox(dis)
  extend<-c(extend[1], extend[3], extend[2], extend[4])
  
  if (between(extend[1], extant_mask[1], extant_mask[2])&
      between(extend[2], extant_mask[1], extant_mask[2])&
      between(extend[3], extant_mask[3], extant_mask[4])&
      between(extend[4], extant_mask[3], extant_mask[4])){
    dis_c<-crop(mountain, extend)
    dis_c<-mask(dis_c, dis)
    v<-values(dis_c)
    v_plain<-length(v[which(v==0)])
    v_mountain<-length(v[which(v==1)])
    if (F){
      plot(st_geometry(dis))
      plot(dis_c, add=T)
    }
    item<-data.frame(sp=sp, plain=v_plain, mountain=v_mountain, group="Birds")
    all_info[[sp]]<-item
  }else{
    print("no overlap, skip")
  }
}


mammals<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")
sp<-mammals$SP[1]

sps<-unique(mammals$SP)

for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", "Mammals", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  dis<-dis[which((dis$presence %in% PRESENCE)&
                   (dis$origin %in% ORIGIN)&
                   (dis$seasonal %in% SEASONAL)),]
  if (nrow(dis)==0){
    next()
  }
  extend<-st_bbox(dis)
  extend<-c(extend[1], extend[3], extend[2], extend[4])
  
  if (between(extend[1], extant_mask[1], extant_mask[2])&
      between(extend[2], extant_mask[1], extant_mask[2])&
      between(extend[3], extant_mask[3], extant_mask[4])&
      between(extend[4], extant_mask[3], extant_mask[4])){
    dis_c<-crop(mountain, extend)
    dis_c<-mask(dis_c, dis)
    v<-values(dis_c)
    v_plain<-length(v[which(v==0)])
    v_mountain<-length(v[which(v==1)])
    if (F){
      plot(st_geometry(dis))
      plot(dis_c, add=T)
    }
    item<-data.frame(sp=sp, plain=v_plain, mountain=v_mountain, group="Mammals")
    all_info[[sp]]<-item
  }else{
    print("no overlap, skip")
  }
}

alal_info<-rbindlist(all_info)
alal_info$type<-"mixed"
alal_info[sp=="Acropternis_orthonyx"]
all_ibnfo[type=="mountain"]
alal_info[plain==0]$type<-"mountain"
alal_info[mountain==0]$type<-"plain"
table(alal_info$type)
nrow(alal_info[mountain_p>0.9])
saveRDS(alal_info, "../../Objects/Island/species_mountain_10km.rda")

all_ibnfo<-readRDS("../../Objects/Island/species_mountain.rda")
