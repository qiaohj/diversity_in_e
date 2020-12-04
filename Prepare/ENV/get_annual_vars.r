#To convert to mm is very easy: simply multiply the field by the amount of seconds for the desired time window wanted, 
#e.g. for mm/day times by 86400.
#Once in mm/day, if you wanted to get mm/yr you could sum the total of precip through the year to get annual precip. 
#Likewise, if you wanted max or min monthly temperatures you would sum the daily temps for the month 
#(usually always 30 days in the protocol, but make sure of this!) and then divide it by 30 to get the mean of the month.
library(raster)
library(dplyr)
numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}
eck4_prj<-"+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")
Ys<-c(1850:2100)
Ms<-c(1:12)
i=10
COMs<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=Ys, stringsAsFactors = F)
template<-"../../Raster/ENV/Monthly/%s/%s/%s/%d/%s_%d.tif"
template_year<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s.tif"
df<-NULL
COMs<-COMs[sample(nrow(COMs)),]
for (i in c(1:nrow(COMs))){
  com<-COMs[i,]
  labels<-""
  if (com$VAR=="pr"){
    labels<-c("sum")
  }
  if (com$VAR=="tasmax"){
    labels<-c("max")
  }
  if (com$VAR=="tasmin"){
    labels<-c("min")
  }
  label<-labels[1]
  if (file.exists(sprintf(template_year, com$GCM, com$SSP, com$VAR, com$Y, label))){
    next()
  }
  for (label in labels){
    stacked<-c()
    mask<-NULL
    for (m in Ms){
      r_file<-sprintf(template, com$GCM, com$SSP, com$VAR, com$Y, label, m)
      print(paste(i, nrow(COMs), r_file))
      stacked<-c(stacked, r_file)
      if (is.null(mask)){
        mask<-raster(r_file)
      }
    }
    r<-stack(stacked)
    r_year<-NULL
    if (com$VAR=="pr"){
      r_year<-calc(r, sum)
      values(r_year)<-86400*values(r_year)
    }
    if (com$VAR=="tasmax"){
      r_year<-calc(r, max)
      values(r_year)<-values(r_year)/10
    }
    if (com$VAR=="tasmin"){
      r_year<-calc(r, min)
      values(r_year)<-values(r_year)/10
    }
    writeRaster(r_year, sprintf(template_year, com$GCM, com$SSP, com$VAR, com$Y, label), overwrite=T)
  }
}