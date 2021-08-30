library(rgdal)
library(raster)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_100km.tif")

if (F){
  #Detect the islands
  points<-readRDS("../../Objects/Island/islands.rda")
  mountain<-readOGR("../../Shape/GMBA mountain inventory V1.2/GMBA Mountain Inventory_v1.2-World.shp")
  mountain<-spTransform(mountain, crs(mask))
  
  mountains<-crop(mask, extent(mountain))
  mountains<-mask(mountains, mountain)
  plot(mountains)
  mountain_index<-extract(mountains, points[, c("x", "y")])
  points$is_mountain<-ifelse(is.na(mountain_index), F, T)
  saveRDS(points, "../../Objects/Island/islands.rda")
  values(mask)[!is.na(values(mask))]<-points$is_mountain
  plot(mask)
  writeRaster(mask, "../../Objects/Island/mountain.tif")
  
  mountain_buff<-buffer(mountain, width=50000)
  
  mountains<-crop(mask, extent(mountain_buff))
  mountains<-mask(mountains, mountain_buff)
  plot(mountains)
  mountain_index<-extract(mountains, points[, c("x", "y")])
  points$is_mountain_buffer<-ifelse(is.na(mountain_index), F, T)
  saveRDS(points, "../../Objects/Island/islands.rda")
  values(mask)[!is.na(values(mask))]<-points$is_mountain_buffer
  plot(mask)
  writeRaster(mask, "../../Objects/Island/mountain_buffer.tif")
  
}

points<-readRDS("../../Objects/Island/islands.rda")
birds<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")
moubtain<-raster("../../Objects/Island/mountain.tif")
sp<-birds$SP[1]
all_info<-list()
sps<-unique(birds$SP)

for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/Dispersal/Birds/%s/initial_disp_exposure_0_dispersal_0.rda", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  if (nrow(dis)==0){
    next()
  }
  lands<-extract(moubtain, dis[, c("x", "y")])
  xx<-data.table(table(lands))
  plain<-ifelse(nrow(xx[lands==0])>0, xx[lands==0]$N, 0)
  mountain<-ifelse(nrow(xx[lands==1])>0, xx[lands==1]$N, 0)
  item<-data.frame(sp=sp, plain=plain, mountain=mountain, group="Bird")
  all_info[[sp]]<-item
}


mammals<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")
sp<-mammals$SP[1]

sps<-unique(mammals$SP)

for (i in c(1:length(sps))){
  print(paste(i, length(sps)))
  sp<-sps[i]
  sp<-gsub(" ", "_", sp)
  target<-sprintf("../../Objects/Dispersal/Mammals/%s/initial_disp_exposure_0_dispersal_0.rda", sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  if (nrow(dis)==0){
    next()
  }
  lands<-extract(moubtain, dis[, c("x", "y")])
  xx<-data.table(table(lands))
  plain<-ifelse(nrow(xx[lands==0])>0, xx[lands==0]$N, 0)
  mountain<-ifelse(nrow(xx[lands==1])>0, xx[lands==1]$N, 0)
  item<-data.frame(sp=sp, plain=plain, mountain=mountain, group="Mammals")
  all_info[[sp]]<-item
}

all_info<-rbindlist(all_info)
all_info$type<-"mixed"
all_info[plain==0]$type<-"mountain"
all_info[mountain==0]$type<-"plain"
table(all_info$type)

saveRDS(all_info, "../../Objects/Island/species_mountain.rda")

all_info<-readRDS("../../Objects/Island/species_mountain.rda")
