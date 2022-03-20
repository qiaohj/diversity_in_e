library(VoCC)
library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_100km.tif")
esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
base<-"../../Raster/Bioclim/%s/%s/%d/bio%d_eck4.tif"
y_range<-c(2021:2100)
bio<-c(1:19)
coms<-expand.grid(esm_ssp=esm_ssp, bio=bio, stringsAsFactors = F)
i=1

for (i in c(1:nrow(coms))){
  item<-coms[i,]
  print(paste(i, nrow(coms), item$esm_ssp, "bio:", item$bio))
  sss_str<-strsplit(item$esm_ssp, "_")[[1]]
  target<-sprintf("../../Raster/VoCC_mask/%s/%s/bio%d_voccMag.tif", sss_str[1], sss_str[2], item$bio)
  if (file.exists(target)){
    next()
  }
  folder<-sprintf("../../Raster/VoCC_mask/%s/%s", sss_str[1], sss_str[2])
  if (!dir.exists(folder)){
    dir.create(folder, recursive = T)
  }
  saveRDS(NULL, target)
  ff<-c()
  for (y in y_range){
    f<-sprintf(base, sss_str[1], sss_str[2], y, item$bio)
    ff<-c(ff, f)
  }
  
  yr<-stack(ff)
  yr<-mask(yr, mask)
  yr<-yr/100
  tr <- tempTrend(yr, th = 10)
  sg <- spatGrad(yr, th = 0.0001, projected = FALSE)
  
  # Magnitude and angle of the climate velocity (km/yr) 1960-2009
  
  v <- gVoCC(tr,sg)
  
  if (F){
    rrr<-raster(gsub("VoCC_mask", "VoCC", target))
    rrr<-mask(rrr, mask)
    plot(rrr)
    plot(v[[1]])
  }
  writeRaster(v[[1]], target, overwrite=T)
  writeRaster(v[[2]], sprintf("../../Raster/VoCC_mask/%s/%s/bio%d_voccAng.tif", sss_str[1], sss_str[2], item$bio), overwrite=T)
}
