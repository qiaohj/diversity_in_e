library(raster)
library(data.table)
library(gdalUtilities)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
if (F){
  r<-raster("../../Raster/2020_LC_Type5.tif")
  r_sinu<-projectRaster(r, res=c(500, 500), crs=crs(mask_10km), method="ngb")
  writeRaster(r_sinu, filename="../../Raster/2020_LC_Type5_sinu.tif")
}
r_sinu<-raster("../../Raster/2020_LC_Type5_sinu_1km.tif")
print(sprintf("Current core number is %d", getDTthreads()))
esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
group_df_birds<-readRDS("../../Data/Birds/bird_df.rda")
group_df_mammals<-readRDS("../../Data/Mammals/mammal_df.rda")
group_df_birds<-data.table(group="Birds", sp=unique(group_df_birds$SCINAME))
group_df_mammals<-data.table(group="Mammals", sp=unique(group_df_mammals$binomial))
group_df<-rbindlist(list(group_df_birds, group_df_mammals))
print("reading mask layers")
mask_10km<-raster("../../Raster/mask_10km.tif")
no_na_mask_10km<-!is.na(values(mask_10km))
mask_points_10km<-data.table(rasterToPoints(mask_10km))
#mask_points_10km<-st_as_sf(mask_points_10km, coords = c("x", "y"), crs = st_crs(mask_10km))

mask_100km<-raster("../../Raster/mask_100km.tif")
mask_points_100km<-data.table(rasterToPoints(mask_100km))
no_na_mask_100km<-!is.na(values(mask_100km))
#mask_points_100km<-st_as_sf(mask_points_100km, coords = c("x", "y"), crs = st_crs(mask_100km))

bi="Pachyptila crassirostris"
coms<-expand.grid(exposure_threshold=c(0, 5), dispersal=c(0, 1))

i=1
j=1
group_df<-group_df[sample(nrow(group_df), nrow(group_df))]
for (i in 1:length(group_df$sp)) {
  start_time<-Sys.time()
  bi<-group_df$sp[i]
  group<-group_df$group[i]
  print(paste(bi, group, i, length(group_df$sp), "Start time:", start_time))
  for (j in 1:nrow(coms)){
    dispersal<-coms[j, "dispersal"]
    exposure_threshold<-coms[j, "exposure_threshold"]
    
    
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
    fit_str<-sprintf("%s/fit.rda", target_folder)
    if (!file.exists(fit_str)){
      next()
    }
    item_str<-esm_ssp[1]
    for (item_str in esm_ssp){
      print(paste(i, nrow(group_df), "disp", dispersal, "expo", exposure_threshold, item_str))
      target<-sprintf("%s/%s_%d_dispersal_%d_lc.rda", target_folder, item_str,
                      exposure_threshold, dispersal)
      if (file.exists(target)){
        next()
      }
      saveRDS(NULL, target)
      source<-sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, item_str,
                      exposure_threshold, dispersal)
      disp_item<-NULL
      res<-NULL
      if (file.exists(source)){
        disp_item<-readRDS(source)
        res<-"10km"
        mask<-mask_10km
        mask_p<-mask_points_10km
        disp_2020<- readRDS(sprintf("%s/initial_disp_10km_exposure_%d_dispersal_%d.rda", 
                                    target_folder, exposure_threshold, dispersal))
        mask_p$v<-NA
        mask_p[mask_10km %in% disp_2020$mask_10km]$v<-1
        values(mask)[no_na_mask_10km]<-mask_p$v
        mask<-crop(mask, c(min(disp_2020$x - 5000), max(disp_2020$x + 5000), 
                           min(disp_2020$y - 5000), max(disp_2020$y + 5000)))
      }
      if (is.null(disp_item)){
        source<-sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                        exposure_threshold, dispersal)
        if (file.exists(source)){
          disp_item<-readRDS(source)
          if (length(disp_item)==0){
            next()
          }
          res<-"100km"
          mask<-mask_100km
          mask_p<-mask_points_100km
          disp_2020<- readRDS(sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda", 
                                target_folder, exposure_threshold, dispersal))
          mask_p$v<-NA
          mask_p[mask_100km %in% disp_2020$mask_100km]$v<-1
          values(mask)[no_na_mask_100km]<-mask_p$v
          mask<-crop(mask, c(min(disp_2020$x - 50000), max(disp_2020$x + 50000), 
                             min(disp_2020$y - 50000), max(disp_2020$y + 50000)))
        }
      }
      if (is.null(disp_item)){
        next()
      }
      mask_2020_file<-sprintf("%s/exposure_%d_dispersal_%d_disp_2020.tif", 
                              target_folder, exposure_threshold, dispersal)
      mask_2020_file_1km<-sprintf("%s/exposure_%d_dispersal_%d_disp_2020_1km.tif", 
                              target_folder, exposure_threshold, dispersal)
      writeRaster(mask, mask_2020_file, datatype="INT1U",
                  overwrite=T)
      gdalwarp(mask_2020_file, mask_2020_file_1km, tr=c(1000, 1000), r="near", ot="Byte", 
               srcnodata=255, dstnodata=255, overwrite=T)
      
      mask_1km<-raster(mask_2020_file_1km)
      lc_2020<-crop(r_sinu, extent(mask_1km))
      extent(lc_2020)<-extent(mask_1km)
      lc_2020<-mask(lc_2020, mask_1km)
      dispersal_log<-list()
      dispersal_log[["lc_2020"]]<-lc_2020
      writeRaster(lc_2020, sprintf("%s/exposure_%d_dispersal_%d_lc_2020.tif",
                                   target_folder, exposure_threshold, dispersal), datatype="INT1U",
                  overwrite=T)
      disp_2100<-disp_item[["2100"]]
      if (is.null(disp_2100)){
        next()
      }
      if (nrow(disp_2100)==0){
        next()
      }
      disp_2100<-disp_2100[suitable==1]
      if (nrow(disp_2100)>0){
        if (res=="100km"){
          mask<-mask_100km
          mask_p<-mask_points_100km
          mask_p$v<-NA
          mask_p[mask_100km %in% disp_2100$mask_100km]$v<-1
          values(mask)[no_na_mask_100km]<-mask_p$v
          mask<-crop(mask, c(min(disp_2100$x)-50000, max(disp_2100$x)+50000, 
                             min(disp_2100$y)-50000, max(disp_2100$y)+50000))
        }else{
          mask<-mask_10km
          mask_p<-mask_points_10km
          mask_p$v<-NA
          mask_p[mask_10km %in% disp_2100$mask_10km]$v<-1
          values(mask)[no_na_mask_10km]<-mask_p$v
          mask<-crop(mask, c(min(disp_2100$x)-5000, max(disp_2100$x)+5000, 
                             min(disp_2100$y)-5000, max(disp_2100$y)+5000))
        }
        
        mask_2100_file<-sprintf("%s/%s_exposure_%d_dispersal_%d_disp_2100.tif", 
                                target_folder, item_str, exposure_threshold, dispersal)
        mask_2100_file_1km<-sprintf("%s/%s_exposure_%d_dispersal_%d_disp_2100_1km.tif", 
                                     target_folder, item_str, exposure_threshold, dispersal)
        writeRaster(mask, mask_2100_file,datatype="INT1U",
                    overwrite=T)
        gdalwarp(mask_2100_file, mask_2100_file_1km, tr=c(1000, 1000), r= "near", ot="Byte", 
                 srcnodata=255, dstnodata=255, overwrite=T)
        
        mask_1km<-raster(mask_2100_file_1km)
        lc_2100<-crop(r_sinu, extent(mask_1km))
        extent(lc_2100)<-extent(mask_1km)
        lc_2100<-mask(lc_2100, mask_1km)
        writeRaster(lc_2100, sprintf("%s/%s_exposure_%d_dispersal_%d_lc_2100.tif", 
                                     target_folder, item_str, exposure_threshold, dispersal),datatype="INT1U",
                    overwrite=T)
        dispersal_log[["lc_2100"]]<-lc_2100
      }
      print("Writing result")
      saveRDS(dispersal_log, target)
      print("Done! Writing result")
    }
  }
  end_time<-Sys.time()
  print(paste("end time", end_time))
}

