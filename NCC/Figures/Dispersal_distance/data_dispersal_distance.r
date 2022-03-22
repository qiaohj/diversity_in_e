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

esm_ssp<-c("MRI-ESM2-0_SSP245")
group_df_birds<-readRDS("../../Data/Birds/bird_df.rda")
group_df_mammals<-readRDS("../../Data/Mammals/mammal_df.rda")
group_df_birds<-data.table(group="Birds", sp=unique(group_df_birds$SCINAME))
group_df_mammals<-data.table(group="Mammals", sp=unique(group_df_mammals$binomial))
group_df<-rbindlist(list(group_df_birds, group_df_mammals))
print("reading mask layers")

bi="Glyphonycteris sylvestris"
coms<-expand.grid(exposure_threshold=c(0), dispersal=c(1))

i=1
j=1
group_df<-group_df[sample(nrow(group_df), nrow(group_df))]
#group_df<-group_df[group=="Mammals"]
all_disp<-list()
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
      source<-sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, item_str,
                      exposure_threshold, dispersal)
      disp_item<-NULL
      res<-NULL
      if (file.exists(source)){
        disp_item<-readRDS(source)
        res<-"10km"
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
        }
      }
      if (is.null(disp_item)){
        next()
      }
      if (length(disp_item)==0){
        next()
      }
      disp_item<-rbindlist(disp_item)
      disp_dist<-disp_item[, .(mean_disp=mean(disp_item$disp), sd_disp=sd(disp_item$disp)),
                           by=list(exposure)]
      disp_dist$res<-res
      disp_dist$sp<-bi
      disp_dist$group<-group
      eeee<-strsplit(item_str, "_")[[1]]
      disp_dist$esm<-eeee[1]
      disp_dist$ssp<-eeee[2]
      all_disp[[length(all_disp)+1]]<-disp_dist
    }
  }
  end_time<-Sys.time()
  print(paste("end time", end_time))
}

all_disp<-rbindlist(all_disp)
saveRDS(all_disp, "../../Objects/Dispersal_distances/all_mean_disp_dist_new.rda")