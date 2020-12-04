library(raster)
library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}
eck4_prj<-"+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")
Ys<-c(1850:2100)
Ms<-c(1:12)
i=10
COMs<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=Ys, stringsAsFactors = F)
template_year<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s.tif"
template_year_eck4<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"
COMs<-COMs[sample(nrow(COMs)),]
for (i in c(1:nrow(COMs))){
  print(paste(i, nrow(COMs)))
  com<-COMs[i,]
  label<-""
  if (com$VAR=="pr"){
    label<-"sum"
  }
  if (com$VAR=="tasmax"){
    label<-"max"
  }
  if (com$VAR=="tasmin"){
    label<-"min"
  }
  target<-sprintf(template_year_eck4, com$GCM, com$SSP, com$VAR, com$Y, label)
  if (file.exists(target)){
    next()
  }
  saveRDS(NA, target)
  r_file<-sprintf(template_year, com$GCM, com$SSP, com$VAR, com$Y, label)
  r<-raster(r_file)
  r_eck4<-projectRaster(r, res=c(1e+05, 1e+05), crs=eck4_prj)
  writeRaster(r_eck4, target, overwrite=T)
}