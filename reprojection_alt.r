library(raster)
library(p_r$gain_loss)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
#proj4string(r)
proj.re<-proj4string(mask)

r<-raster("../../Raster/ALT/alt.tif")
r2<-projectRaster(r, crs=proj.re, res=res(mask))

writeRaster(r2, 
            "../../Raster/ALT/alt_eck4.tif", 
            format="GTiff", overwrote=T)


slope<-terrain(r, opt="slope")
plot(slope)

r2<-projectRaster(slope, crs=proj.re, res=res(mask))

writeRaster(r2, 
            "../../Raster/ALT/slope_eck4.tif", 
            format="GTiff", overwrote=T)
