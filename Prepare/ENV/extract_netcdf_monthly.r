library(raster)
library(sp)
library(ncdf4)
library(maptools)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
base<-"../../Data/Raw/ENV"
target<-"../../Raster"
files<-list.files(base, pattern="\\.nc", recursive = T, full.names = T)
files<-read.table(text = files, sep="/", stringsAsFactors=FALSE)
colnames(files)[6:8]<-c("GCM", "SSP", "VAR")
files
args = commandArgs(trailingOnly=TRUE)

index<-as.numeric(args[1])

#index<-1
netcdf_f<-sprintf("%s/%s/%s/%s", base, files[index, "GCM"], files[index, "SSP"], files[index, "VAR"])
woa13<-nc_open(netcdf_f)
lat<-woa13$dim$lat$vals
lon<-woa13$dim$lon$vals
time<-woa13$dim$time$vals
varname<-names(woa13$var)[length(names(woa13$var))]
date1<-as.Date("1850/1/1", format="%Y/%m/%d")
date2<-as.Date("2100/12/31", format="%Y/%m/%d")
date2-date1
length(time)
year_length<-2100-1850+1
if (files[index, "GCM"]=="UKESM1"){
  days_in_month<-(max(time)-min(time)+1)/year_length/12  
}else{
  days_in_month<-NA
}

lon_lat<-expand.grid(lon, lat)
v_woa13<-ncvar_get(woa13, varid=varname)
t=1
d_all<-data.frame(lon=as.vector(lon_lat$Var1),
                  lat=as.vector(lon_lat$Var2))

points<-SpatialPointsDataFrame(coords=data.frame(lon=d_all$lon, lat=d_all$lat),
                               data=d_all)
y=0
m=11
for (y in c(0:(year_length-1))){
  target_folder<-sprintf("%s/%s/%s/%s/%d", target, files[index, "GCM"], files[index, "SSP"], 
                         varname, y+1850)
  dir.create(target_folder, showWarnings = F, recursive = T)
  for (m in c(0:11)){
    if (is.na(days_in_month)){
      from_day<-as.numeric(as.Date(sprintf("%d/%d/1",y+1850,m+1), format="%Y/%m/%d") - date1 + 1)
      if (m==11){
        to_day<-as.numeric(as.Date(sprintf("%d/%d/1",y+1851, 1), format="%Y/%m/%d") - date1)
      }else{
        to_day<-as.numeric(as.Date(sprintf("%d/%d/1",y+1850,m+2), format="%Y/%m/%d") - date1)
      }
    }else{
      from_day<-(y*12+m)*days_in_month+1
      to_day<-(y*12+m+1)*days_in_month
    }
    print(paste("GCM:", files[index, "GCM"], 
                "SSP:", files[index, "SSP"],
                "VAR:", varname,
                "Y:", y, "M:", m, "FROM:", from_day, "TO:", to_day))
    v_matrix<-v_woa13[,,c(from_day:to_day)]
    if (varname=="pr"){
      sum_v<-apply(v_matrix, 1:2, sum, na.rm = TRUE)
      
      r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
      r <- rasterize(points, r, as.vector(sum_v), fun=mean)
      r  <- rotate(r)
      #plot(r)
      projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
      writeRaster(r, filename=sprintf("%s/sum_%d.tif", target_folder, m+1), format="GTiff", overwrite=TRUE)
    }
    if (varname=="tasmax"){
      max_v<-apply(v_matrix, 1:2, max, na.rm = TRUE)
      
      r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
      r <- rasterize(points, r, as.vector(max_v), fun=mean)
      r  <- rotate(r)
      #plot(r)
      projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
      writeRaster(r, filename=sprintf("%s/max_%d.tif", target_folder, m+1), format="GTiff", overwrite=TRUE)
      
      mean_v<-apply(v_matrix, 1:2, mean, na.rm = TRUE)
      
      r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
      r <- rasterize(points, r, as.vector(mean_v), fun=mean)
      r  <- rotate(r)
      #plot(r)
      projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
      writeRaster(r, filename=sprintf("%s/mean_%d.tif", target_folder, m+1), format="GTiff", overwrite=TRUE)
    }
    if (varname=="tasmin"){
      min_v<-apply(v_matrix, 1:2, min, na.rm = TRUE)
      
      r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
      r <- rasterize(points, r, as.vector(min_v), fun=mean)
      r  <- rotate(r)
      #plot(r)
      projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
      writeRaster(r, filename=sprintf("%s/min_%d.tif", target_folder, m+1), format="GTiff", overwrite=TRUE)
      
      mean_v<-apply(v_matrix, 1:2, mean, na.rm = TRUE)
      
      r <- raster(ncols=length(lon), nrows= length(lat), xmn=0, xmx=360, ymn=-90, ymx=90)
      r <- rasterize(points, r, as.vector(mean_v), fun=mean)
      r  <- rotate(r)
      #plot(r)
      projection(r) <- "+proj=longlat +datum=WGS84 +no_defs+towgs84=0,0,0"
      writeRaster(r, filename=sprintf("%s/mean_%d.tif", target_folder, m+1), format="GTiff", overwrite=TRUE)
    }
  }
}
