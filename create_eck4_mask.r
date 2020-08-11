library(raster)
eck4_prj<-"+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
mask<-raster("../../Raster/mask_high_res.tif")
crs(mask)<-eck4_prj
new_mask<-projectRaster(mask, res=c(100000, 100000), crs=eck4_prj, method="ngb")
writeRaster(new_mask, "../../Raster/mask.tif", overwrite=T)
plot(new_mask)
