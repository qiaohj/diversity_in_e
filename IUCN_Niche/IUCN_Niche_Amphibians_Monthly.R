library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(ntbox)
setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
#for AMPHIBIANS
sp_df<-readOGR("../Shape/iucn_species_Ranges/AMPHIBIANS", "AMPHIBIANS_eck4") 
mask_bak<-raster("../Raster/mask.tif")

unique <- unique(sp_df@data$binomial)
unique<-as.character(unique)


NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}

tmax<-list()
tmin<-list()
prec<-list()
for (month in c(1:12)){
  tmax[[month]]<-raster(sprintf("../Raster/Bioclim 2.1/Present/variable/tmax_%d.tif", month))
}
i=1
for (i in 1:length(unique)) {
  
  bi<-unique[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target<-sprintf("../Object/IUCN_Distribution/Monthly/Amphibians/%s.rda", gsub(" ", "_", bi))
  if (file.exists(target)){
    next()
  }
  saveRDS(NA, file=target)
  
  print("extracting the matched polygons")
  tmp <- sp_df[sp_df$binomial == bi, ]
  
  mask<-mask_bak
  print("rasterizing to raster")
  mask<-raster(extent(tmp), res=res(mask_bak), crs=crs(mask_bak))
  rp <- rasterize(tmp, mask)
  #
  #plot(sp_df)
  #plot(tmp, add=T)
  #plot(rp)
  
  print("saving result")
  no_na<-!is.na(values(rp))
  n_pixel<-length(values(rp)[no_na])
  if (n_pixel==0){
    mask<-raster(extent(tmp@bbox), res=res(mask_bak), crs=crs(mask_bak))
    rp <- rasterize(tmp@bbox, mask)
    no_na<-!is.na(values(rp))
    n_pixel<-length(values(rp)[no_na])
  }
  if (n_pixel>0){
    values(rp)[no_na]<-1
    ppp<-data.frame(rasterToPoints(rp))
    print(target)
    vartypes<-c("bioc", "prec", "tmax", "tmin")
    month=1
    for (month in (1:12)){
      ppp1<-ppp
      ppp1$temp<-extract(), ppp[, c("x", "y")])
      ppp1$prec<-extract(bio13, ppp[, c("x", "y")])
    }
    #ppp$bio1<-extract(bio1, ppp[, c("x", "y")])
    
    
    
    #ppp$bio12<-extract(bio12, ppp[, c("x", "y")])
    ppp2<-ppp
    ppp2$temp<-extract(bio6, ppp[, c("x", "y")])
    ppp2$prec<-extract(bio14, ppp[, c("x", "y")])
    
    ppp<-ppp[complete.cases(ppp),]
    saveRDS(ppp, target)
    #writeRaster(rp, sprintf("../Raster/IUCN_Distribution/Amphibians/%s.tif", gsub(" ", "_", bi)), overwrite=T)
  }
  n_pixel<-nrow(ppp)
  if (n_pixel>7){
    fit <- cov.rob(ppp[, c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")], 
                   quantile.used=NDquntil(nrow(ppp), 0.95),  method = "mve")
    saveRDS(fit, gsub("\\.rda", "\\.fit\\.rda", target))
    
  }
  rp<-NA
  tmp<-NA
  gc()
}


r<-stack("/Users/huijieqiao/temp/share 2/spatial03/worldclim/cmip6/7_fut/10m/IPSL-CM6A-LR/ssp370/wc2.1_10m_tmax_IPSL-CM6A-LR_ssp370_2081-2100.tif")
xxx<-extract(r, ppp[, c("x", "y")])
