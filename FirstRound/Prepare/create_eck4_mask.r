library(raster)
eck4_prj<-"+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
mask<-raster("../../Raster/mask_high_res.tif")
crs(mask)<-eck4_prj
new_mask<-projectRaster(mask, res=c(100000, 100000), crs=eck4_prj, method="ngb")
writeRaster(new_mask, "../../Raster/mask.tif", overwrite=T)
plot(new_mask)

bio1<-raster("../../Raster/bioclim/wc2.0_bio_30s_01.tif")
new_bio1<-projectRaster(bio1, res=c(100000, 100000), crs=eck4_prj)
writeRaster(new_bio1, "../../Raster/bioclim/bio1_eck4.tif", overwrite=T)

new_bio1<-projectRaster(bio1, res=c(10000, 10000), crs=eck4_prj)
writeRaster(new_bio1, "../../Raster/bioclim/bio1_eck4_10km.tif", overwrite=T)

bio12<-raster("../../Raster/bioclim/wc2.0_bio_30s_12.tif")
new_bio12<-projectRaster(bio12, res=c(10000, 10000), crs=eck4_prj)
writeRaster(new_bio12, "../../Raster/bioclim/bio12_eck4_10km.tif", overwrite=T)
