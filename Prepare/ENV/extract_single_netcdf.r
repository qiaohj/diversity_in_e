library(raster)
library(sp)
library(ncdf4)
library(maptools)
setwd("~/Downloads")
woa13<-nc_open("pr_day_UKESM1-0-LL_historical_ssp119_r1i1p1f2_gn_18500101-21001230.nc")
lat<-woa13$dim$lat$vals
lon<-woa13$dim$lon$vals
time<-woa13$dim$time$vals


lon_lat<-expand.grid(lon, lat)
v_woa13<-ncvar_get(woa13, varid="pr")
t=1
d_all<-data.frame(lon=as.vector(lon_lat$Var1),
                  lat=as.vector(lon_lat$Var2))

points<-SpatialPointsDataFrame(coords=data.frame(lon=d_all$lon, lat=d_all$lat),
                               data=d_all)

for (t in c(1:length(time))){
  
  tt<-time[t]
  print(paste(t, length(time), tt))
  
  r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
  r <- rasterize(points, r, as.vector(v_woa13[,,t]), fun=mean)
  r  <- rotate(r)
  #plot(r)
  projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
  writeRaster(r, filename=sprintf("~/Downloads/prec/%d.tif", tt*10), format="GTiff", overwrite=TRUE)
}
t=1
for (t in c(1:length(time))){
  tt<-time[t]
  print(paste(t, length(time), tt))
  pngfile<-sprintf("~/Downloads/prec_png/%d.png", tt*10)
  if (file.exists(pngfile)){
    next()
  }
  r<-raster(sprintf("~/Downloads/prec/%d.tif", tt*10), format="GTiff", overwrite=TRUE)
  png(filename=pngfile, width = 800, height = 480)
  plot(r)
  dev.off()
}
