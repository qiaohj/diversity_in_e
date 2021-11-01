library(dismo)
library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")

start_range<-c(1850:2100)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=start_range)
i=1
mask<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/mask_index.tif")
2616/6777
for (i in c(2501:2616)){
#for (i in c(1:nrow(start_layer_df))){
  print(paste("Init layer list:", i, nrow(start_layer_df)))
  item<-start_layer_df[i,]
  var_tamplate<-"/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/ENV/Monthly/%s/%s/%s/%d/%s_%d.tif"
  
  pr<-stack(sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum", c(1:12)))
  values(pr)<-86400*values(pr)
  tmax<-stack(sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max", c(1:12)))
  values(tmax)<-values(tmax)-273.16
  tmin<-stack(sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min", c(1:12)))
  values(tmin)<-values(tmin)-273.16
  bioc<-biovars(pr, tmin, tmax)
  target_folder<-sprintf("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d", 
                         item$GCM, item$SSP, item$Y)
  dir.create(target_folder, recursive = T, showWarnings = F)
  j=1
  for (j in c(1:19)){
    rrr<-bioc[[j]]
    writeRaster(rrr, sprintf("%s/bio%d.tif", target_folder, j), overwrite=T)
    rr<-projectRaster(rrr, mask)
    values(rr)<-round(values(rr) * 100)
    writeRaster(rr, sprintf("%s/bio%d_eck4.tif", target_folder, j), overwrite=T)
  }
}
