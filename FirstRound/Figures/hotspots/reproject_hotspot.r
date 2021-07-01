library(raster)
library(rgdal)
p <- shapefile('../../Shape/hotspots_2016_1/hotspots_2016_1.shp')

p_sub<-p[which(p$Type=="hotspot area"),]
pgeo <- spTransform(p_sub, CRS('+proj=longlat +datum=WGS84'))
ext <- floor(extent(pgeo))
rr <- raster(ext, res=0.1)
crs(rr)<-'+proj=longlat +datum=WGS84'
rr <- rasterize(pgeo, rr, field=1)
plot(rr)
writeRaster(rr, "../../Raster/hotspots_2016_1/hotspots_2016_1_ll_0.1.tif", overwrite=T)
mask<-raster("../../Raster/mask_index.tif")

rr_2<-projectRaster(rr, crs=crs(mask), res=res(mask))
writeRaster(rr_2, "../../Raster/hotspots_2016_1/hotspots_2016_eck4.hotspot_area.tif", overwrite=T)


pgeo <- spTransform(p, CRS('+proj=longlat +datum=WGS84'))
ext <- floor(extent(pgeo))
rr <- raster(ext, res=0.1)
crs(rr)<-'+proj=longlat +datum=WGS84'
rr <- rasterize(pgeo, rr, field=1)
plot(rr)
writeRaster(rr, "../../Raster/hotspots_2016_1/hotspots_2016_1_ll_all_0.1.tif", overwrite=T)
mask<-raster("../../Raster/mask_index.tif")

rr_2<-projectRaster(rr, crs=crs(mask), res=res(mask))
writeRaster(rr_2, "../../Raster/hotspots_2016_1/hotspots_2016_eck4.all.tif", overwrite=T)
