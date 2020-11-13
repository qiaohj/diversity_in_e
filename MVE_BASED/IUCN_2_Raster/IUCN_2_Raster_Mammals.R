library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(sf)
library(fasterize)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#for MAMMALS
mask_bak<-raster("../../Raster/mask.tif")
if (F){
  sp_df<-readOGR("../../Data/Raw/IUCN/MAMMALS", "MAMMALS_TERRESTRIAL_ONLY") 
  
  sp_df_eck4<-spTransform(sp_df, crs(mask_bak))
  writeOGR(sp_df_eck4, "../../Data/Raw/IUCN/MAMMALS", "MAMMALS_TERRESTRIAL_ONLY_ECK4", driver="ESRI Shapefile")
  
  mammals_df<-data.frame(sp_df_eck4)
  mammals_df$geometry<-NULL
  colnames(mammals_df)
  saveRDS(mammals_df, "../../Data/Raw/IUCN/MAMMALS/mammals_df.rda")
  
}
sp_df_eck4<-sf::st_read("../../Data/Raw/IUCN/MAMMALS", layer="MAMMALS_TERRESTRIAL_ONLY_ECK4") 
mammals_df<-readRDS("../../Data/Raw/IUCN/MAMMALS/mammals_df.rda")

unique <- unique(mammals_df$binomial)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
i=100
bi<-"Acerodon humilis"
for (i in 1:length(unique)) {
  
  bi<-unique[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target<-sprintf("../../Objects/IUCN_Distribution/Mammals/%s.rda", gsub(" ", "_", bi))
  if (file.exists(target)){
    tt<-readRDS(target)
    nr<-nrow(tt)
    if (is.null(nr)){
      nr<-0
    }
    if (nr>0){
      next()  
    }
  }
  saveRDS(NA, file=target)
  
  print("extracting the matched polygons")
  tmp_sf <- sp_df_eck4[which((mammals_df$binomial == bi)&(mammals_df$presence %in% PRESENCE)&
                               (mammals_df$origin %in% ORIGIN)&(mammals_df$seasonal %in% SEASONAL)), ]
  if (nrow(tmp_sf)==0){
    next()
  }
  mask<-mask_bak
  print("rasterizing to raster")
  mask<-raster(extent(tmp_sf), res=res(mask_bak), crs=crs(mask_bak))
  rp <- fasterize(tmp_sf, mask, field="presence")
  v<-values(rp)
  if (length(v[!is.na(v)])==0){
    values(rp)<-1
  }
  #
  #plot(tmp_sf)
  #plot(tmp, add=T, col="red")
  #plot(rp)
  
  #plot(mask_bak)
  #plot(mask, add=T, col="blue")
  #plot(rp, add=T, col="red")
  
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

