library(sf)
library(raster)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_100km.tif")

if (F){
  #Detect the islands
  points<-data.table(rasterToPoints(mask))
  xsize<-dim(mask)[2]
  points$x_int<-points$mask_100km %% xsize
  points$y_int<-floor(points$mask_100km/xsize)
  
  plot(points$x_int, points$y_int)
  
  mark_neighbers<-function(p_x, p_y, p_group){
    neighbers<-points[between(x_int, p_x-1, p_x+1)&
                        between(y_int, p_y-1, p_y+1)&
                        (group==-1)]
    #print(neighbers)
    points[mask_100km %in% neighbers$mask_100km]$group<<-p_group
    if (nrow(neighbers)>0){
      for (i in c(1:nrow(neighbers))){
        it<-points[mask_100km==neighbers[i,]$mask_100km]
        #points[mask_100km==it$mask_100km, "group"]<-p_group
        mark_neighbers(it$x_int, it$y_int, p_group)
      }
    }
    
  }
  group_index<-1
  points$group<--1
  while (T){
    item<-points[group==-1]
    if (nrow(item)==0){
      break()
    }
    item<-item[1,]
    points[mask_100km==item$mask_100km, "group"]<-group_index
    p_x<-item$x_int
    p_y<-item$y_int
    p_group<-group_index
    mark_neighbers(item$x_int, item$y_int, group_index)
    group_index<-group_index+1
    print(paste(group_index, nrow(points[group==-1])))
  }
  ggplot(points)+geom_tile(aes(x=x, y=y, fill=factor(group)))+
    theme(legend.position = NULL)
  saveRDS(points, "../../Objects/Island/islands.rda")
  values(mask)[!is.na(values(mask))]<-points$group
  plot(mask)
  writeRaster(mask, "../../Objects/Island/islands_index.tif")
  threshold_group<-54  
  point_group<-data.table(table(points$group))
  threshold<-point_group[V1==threshold_group]$N
  points$is_island<-ifelse(points$group %in% point_group[N<=threshold]$V1, T, F)
  ggplot(points)+geom_tile(aes(x=x, y=y, fill=factor(is_island)))+
    theme(legend.position = NULL)
  saveRDS(points, "../../Objects/Island/islands.rda")
  values(mask)[!is.na(values(mask))]<-points$is_island
  plot(mask)
  writeRaster(mask, "../../Objects/Island/islands.tif", overwrite=T)
}

points<-readRDS("../../Objects/Island/islands.rda")
birds<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")
islands<-raster("../../Objects/Island/islands.tif")
sp<-birds$SP[1]
all_info<-list()
sps<-unique(birds$SP)
i=1
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
  lands<-extract(islands, dis[, c("x", "y")])
  xx<-data.table(table(lands))
  continent<-ifelse(nrow(xx[lands==0])>0, xx[lands==0]$N, 0)
  island<-ifelse(nrow(xx[lands==1])>0, xx[lands==1]$N, 0)
  item<-data.frame(sp=sp, continent=continent, island=island, group="Bird")
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
  lands<-extract(islands, dis[, c("x", "y")])
  xx<-data.table(table(lands))
  continent<-ifelse(nrow(xx[lands==0])>0, xx[lands==0]$N, 0)
  island<-ifelse(nrow(xx[lands==1])>0, xx[lands==1]$N, 0)
  item<-data.frame(sp=sp, continent=continent, island=island, group="Mammals")
  all_info[[sp]]<-item
}

all_info<-rbindlist(all_info)
all_info$type<-"mixed"
all_info[continent==0]$type<-"island"
all_info[island==0]$type<-"continent"
table(all_info$type)
saveRDS(all_info, "../../Objects/Island/species_island.rda")

all_info<-readRDS("../../Objects/Island/species_island.rda")



