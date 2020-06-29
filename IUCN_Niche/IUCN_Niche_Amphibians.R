library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(ntbox)

#for AMPHIBIANS
sp_df<-readOGR("../Shape/iucn_species_Ranges/AMPHIBIANS", "AMPHIBIANS_eck4") 
mask_bak<-raster("../Raster/mask.tif")

unique <- unique(sp_df@data$binomial)
unique<-as.character(unique)


pc1<-raster("../Raster/Bioclim/PCs/Present/pc1.tif")
pc2<-raster("../Raster/Bioclim/PCs/Present/pc2.tif")
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}

i=1
for (i in 1:length(unique)) {
  
  bi<-unique[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target<-sprintf("../Object/IUCN_Distribution/Amphibians/%s.rda", gsub(" ", "_", bi))
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
    saveRDS(ppp, target)
    writeRaster(rp, sprintf("../Raster/IUCN_Distribution/Amphibians/%s.tif", gsub(" ", "_", bi)), overwrite=T)
  }
  if (n_pixel>3){
    ppp$pc1<-extract(pc1, ppp[, c("x", "y")])
    ppp$pc2<-extract(pc2, ppp[, c("x", "y")])
    ppp<-ppp[complete.cases(ppp),]
    saveRDS(ppp, target)
    fit <- cov.rob(ppp[, c("pc1", "pc2")], quantile.used=NDquntil(nrow(ppp), 0.95),  method = "mve")
    saveRDS(fit, gsub("\\.rda", "\\.fit\\.rda", target))
    mve<-cov_center(data = ppp, mve = T, level = 0.95, vars = 4:5)
    saveRDS(mve, gsub("\\.rda", "\\.mve\\.rda", target))
    best_ellipse <- ellipsoidhull(as.matrix(ppp[fit$best, c("pc1", "pc2")] ))
    saveRDS(best_ellipse, gsub("\\.rda", "\\.best_ellipse\\.rda", target))
    #plot(ppp$pc1, ppp$pc2, pch=".")
    #lines(predict(best_ellipse), col="blue")
  }
  rp<-NA
  tmp<-NA
  gc()
}
