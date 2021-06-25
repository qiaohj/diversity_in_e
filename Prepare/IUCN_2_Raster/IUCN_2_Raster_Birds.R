library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(sf)
library(fasterize)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#for BIRDS
mask_bak<-raster("../../Raster/mask.tif")
if (F){
  mask_bak<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/EC-Earth3-Veg/SSP119/1850/bio12_eck4_1km.tif")
  mask_high<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/mask_high_res.tif")
  vessel <- sf::st_read(dsn = "/media/huijieqiao/QNAS/Sp_Richness_GCM/Data/Raw/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
  sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_bak))
  #st_write(sp_df_eck4, "../../Data/Raw/IUCN/BIRDS/BIRDS_ECK4.gdb")
  bird_df<-data.frame(sp_df_eck4)
  bird_df$Shape<-NULL
  
  colnames(bird_df)
  saveRDS(bird_df, "../../Data/Raw/IUCN/BIRDS/bird_df.rda")
}

vessel <- sf::st_read(dsn = "../../Data/Raw/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_bak))
bird_df<-readRDS("../../Data/Raw/IUCN/BIRDS/bird_df.rda")

unique <- unique(bird_df$SCINAME)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
i=511

if (F){
  stat<-NULL
  for (i in 1:length(unique)) {
    bi<-unique[i]
    print(paste(i, length(unique), bi))
    tmp_sf <- sp_df_eck4[which((bird_df$SCINAME == bi)&(bird_df$PRESENCE %in% PRESENCE)&
                                 (bird_df$ORIGIN %in% ORIGIN)&(bird_df$SEASONAL %in% SEASONAL)), ]
    if (class(tmp_sf$Shape)[1]=="sfc_MULTISURFACE"){
      asdfasdf
    }
    plot(st_geometry(tmp_sf))
    
    tmp_sf %>% st_collection_extract("MULTIPOLYGON")
    st_write(tmp_sf, sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)), append=FALSE)
    poly<-st_read(sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)))
    head(st_as_text(tmp_sf$Shape[1]))
    st_cast(st_sfc(tmp_sf), "MULTIPOLYGON")
    st_cast(st_sf(a = 1, st_sfc(tmp_sf$Shape)), "MULTIPOLYGON")
    
    sum_area<-sum(tmp_sf$Shape_Area)
    n<-nrow(tmp_sf)
    type<-class(tmp_sf$Shape)[1]
    item<-data.frame(sp=bi, sum_area=sum_area, n=n, type=type)
    if (is.null(stat)){
      stat<-item
    }else{
      stat<-bind_rows(stat, item)
    }
    unique(stat$type)
    stat%>%filter(type=="sfc_GEOMETRY")
    stat%>%filter(type=="sfc_MULTISURFACE")
    stat%>%filter(type=="sfc_MULTIPOLYGON")
    mask<-raster(extent(tmp_sf), res=res(mask_bak), crs=crs(mask_bak))
    st_crs(tmp_sf)<-crs(mask_bak)
    extract(mask_bak, tmp_sf)
  }
}
for (i in 1:length(unique)) {
  
  bi<-unique[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target<-sprintf("../../Objects/IUCN_Distribution/Birds/%s.rda", gsub(" ", "_", bi))
  if (file.exists(target)){
    n<-nrow(readRDS(target))
    if (is.null(n)){
      n<-0
    }
    if (n>0){
      next()
    }
  }
  saveRDS(NA, file=target)
  
  print("extracting the matched polygons")
  tmp_sf <- sp_df_eck4[which((bird_df$SCINAME == bi)&(bird_df$PRESENCE %in% PRESENCE)&
                               (bird_df$ORIGIN %in% ORIGIN)&(bird_df$SEASONAL %in% SEASONAL)), ]
  if ((class(tmp_sf$Shape)[1]=="sfc_GEOMETRY")|(class(tmp_sf$Shape)[1]=="sfc_MULTISURFACE")){
    st_write(tmp_sf, sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)), append=FALSE)
    tmp_sf<-st_read(sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)))
  }
  if (nrow(tmp_sf)==0){
    next()
  }
  mask<-mask_bak
  print("rasterizing to raster")
  mask<-raster(extent(tmp_sf), res=res(mask_bak), crs=crs(mask_bak))
  rp <- fasterize(tmp_sf, mask)
  v<-values(rp)
  if (length(v[!is.na(v)])==0){
    values(rp)<-1
  }
  #
  #plot(tmp_sf)
  #plot(tmp, add=T, col="red")
  #plot(rp), add=T, col="blue")
  
  print("saving result")
  
  #n_pixel<-length(values(rp))
  
  ppp<-data.frame(rasterToPoints(rp))
  ppp$layer<-extract(mask_bak, ppp[, c("x", "y")])
  if (nrow(ppp[which(!is.na(ppp$layer)),])==0){
    x<-ppp$x+c(100000, 0, -100000)
    y<-ppp$y+c(100000, 0, -100000)
    ppp<-expand.grid(x=x, y=y)
    ppp$layer<-extract(mask_bak, ppp[, c("x", "y")])
    ppp<-ppp[complete.cases(ppp),]
    if (nrow(ppp)>1){
      ppp<-ppp[1,]
    }
  }
  ppp<-ppp[which(!is.na(ppp$layer)),]
  #plot(ppp$x, ppp$y)
  print(target)
  saveRDS(ppp, target)
}

