library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
library(rmapshaper)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))
mask_10km<-raster("../../Raster/mask_10km.tif")
#mask_1km<-raster("../../Raster/mask_1km.tif")

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2))
}

bird_df<-readRDS("../../Data/Birds/bird_df.rda")
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
bird_full<-merge(bird_df, bird_disp, by.x="SCINAME", by.y="iucn_name", all=F)

unique <- unique(bird_full$SCINAME)
unique<-as.character(unique)


print("reading env layers")

start_env_layers_se<-readRDS("../../Objects/stacked_layers_1850_2020_df_10km.rda")
x_size<-dim(mask_10km)[2]
bi<-"Tangara rufigenis"
bird_full_sum_area<-bird_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
bird_full_sum_area<-bird_full_sum_area[order(sum_are),]
bi<-bird_full_sum_area[bird_full_sum_area$sum_are<=1.5*min(bird_full_sum_area$sum_are)]$SCINAME[1]
bird_full_sum_area<-bird_full_sum_area[sample(nrow(bird_full_sum_area), nrow(bird_full_sum_area))]
ebird_base<-"/media/huijieqiao/SSD_Fast/ES50_eBird/Tables/Speccies_Split_202112"
i=2
for (i in 5269:length(bird_full_sum_area$SCINAME)) {
  
  bi<-bird_full_sum_area$SCINAME[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  ebird_full<-sprintf("%s/%s", ebird_base, bi)
  if (!file.exists(ebird_full)){
    next()
  }
  
  target<-sprintf("../../Objects/Dispersal/Birds/%s/fit_ebird.rda", gsub(" ", "_", bi))
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  ff<-""
  cols_ll<-c("SCIENTIFIC_NAME", "LATITUDE", "LONGITUDE", "YEAR")
  all_items<-list()
  for (ff in list.files(ebird_full, pattern="\\.rda")){
    item<-readRDS(sprintf("%s/%s/%s", ebird_base, bi, ff))
    item<-item[APPROVED==1|REVIEWED==1]
    item<-item[, ..cols_ll]
    item<-unique(item)
    if (nrow(item)>0){
      all_items[[ff]]<-item
    }
  }
  all_items<-rbindlist(all_items)
  cols_ll<-c("LONGITUDE", "LATITUDE")
  points<-st_as_sf(x=all_items, coords = cols_ll, 
                   crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  points<-st_transform(points, crs=st_crs(mask_10km))
  all_items$mask_10km<-raster::extract(mask_10km, st_coordinates(points))
  all_items_n<-all_items[, .(N=.N), by=list(mask_10km, YEAR)]
  
  #nb_v<-start_env_layers_se[mask_10km %in% nb_p$mask_10km]
  nb_v<-merge(all_items_n, start_env_layers_se, by.x=c("mask_10km", "YEAR"),
              by.y=c("mask_10km", "year"))
  
  cols<-c("x", "y", "mask_10km")
  mean_bio1<-mean(nb_v$bio1, na.rm=T)
  sd_bio1<-sd(nb_v$bio1, na.rm=T)
  quantile_bio1<-quantile(nb_v$bio1, c(0.25, 0.75), na.rm=T)
  IQR_bio1<-quantile_bio1[2]-quantile_bio1[1]
  range_bio1_sd_min<-mean_bio1-3*sd_bio1
  range_bio1_sd_max<-mean_bio1+3*sd_bio1
  bio1<-nb_v[between(bio1, range_bio1_sd_min, range_bio1_sd_max)]$bio1
  range_bio1_sd_min<-min(bio1)
  range_bio1_sd_max<-max(bio1)
  
  range_bio1_IQR_min<-mean_bio1-1.5*IQR_bio1
  range_bio1_IQR_max<-mean_bio1+1.5*IQR_bio1
  
  mean_bio5<-mean(nb_v$bio5, na.rm=T)
  sd_bio5<-sd(nb_v$bio5, na.rm=T)
  quantile_bio5<-quantile(nb_v$bio5, c(0.25, 0.75), na.rm=T)
  IQR_bio5<-quantile_bio5[2]-quantile_bio5[1]
  range_bio5_sd_min<-mean_bio5-3*sd_bio5
  range_bio5_sd_max<-mean_bio5+3*sd_bio5
  bio5<-nb_v[between(bio5, range_bio5_sd_min, range_bio5_sd_max)]$bio5
  range_bio5_sd_min<-min(bio5)
  range_bio5_sd_max<-max(bio5)
  
  range_bio5_IQR_min<-mean_bio5-1.5*IQR_bio5
  range_bio5_IQR_max<-mean_bio5+1.5*IQR_bio5
  
  mean_bio6<-mean(nb_v$bio6, na.rm=T)
  sd_bio6<-sd(nb_v$bio6, na.rm=T)
  quantile_bio6<-quantile(nb_v$bio6, c(0.25, 0.75), na.rm=T)
  IQR_bio6<-quantile_bio6[2]-quantile_bio6[1]
  range_bio6_sd_min<-mean_bio6-3*sd_bio6
  range_bio6_sd_max<-mean_bio6+3*sd_bio6
  bio6<-nb_v[between(bio6, range_bio6_sd_min, range_bio6_sd_max)]$bio6
  range_bio6_sd_min<-min(bio6)
  range_bio6_sd_max<-max(bio6)
  
  range_bio6_IQR_min<-mean_bio6-1.5*IQR_bio6
  range_bio6_IQR_max<-mean_bio6+1.5*IQR_bio6
  
  mean_bio12<-mean(nb_v$bio12, na.rm=T)
  sd_bio12<-sd(nb_v$bio12, na.rm=T)
  quantile_bio12<-quantile(nb_v$bio12, c(0.25, 0.75), na.rm=T)
  IQR_bio12<-quantile_bio12[2]-quantile_bio12[1]
  range_bio12_sd_min<-mean_bio12-3*sd_bio12
  range_bio12_sd_max<-mean_bio12+3*sd_bio12
  bio12<-nb_v[between(bio12, range_bio12_sd_min, range_bio12_sd_max)]$bio12
  range_bio12_sd_min<-min(bio12)
  range_bio12_sd_max<-max(bio12)
  
  range_bio12_IQR_min<-mean_bio12-1.5*IQR_bio12
  range_bio12_IQR_max<-mean_bio12+1.5*IQR_bio12
  
  mean_bio13<-mean(nb_v$bio13, na.rm=T)
  sd_bio13<-sd(nb_v$bio13, na.rm=T)
  quantile_bio13<-quantile(nb_v$bio13, c(0.25, 0.75), na.rm=T)
  IQR_bio13<-quantile_bio13[2]-quantile_bio13[1]
  range_bio13_sd_min<-mean_bio13-3*sd_bio13
  range_bio13_sd_max<-mean_bio13+3*sd_bio13
  bio13<-nb_v[between(bio13, range_bio13_sd_min, range_bio13_sd_max)]$bio13
  range_bio13_sd_min<-min(bio13)
  range_bio13_sd_max<-max(bio13)
  
  range_bio13_IQR_min<-mean_bio13-1.5*IQR_bio13
  range_bio13_IQR_max<-mean_bio13+1.5*IQR_bio13
  
  mean_bio14<-mean(nb_v$bio14, na.rm=T)
  sd_bio14<-sd(nb_v$bio14, na.rm=T)
  quantile_bio14<-quantile(nb_v$bio14, c(0.25, 0.75), na.rm=T)
  IQR_bio14<-quantile_bio14[2]-quantile_bio14[1]
  range_bio14_sd_min<-mean_bio14-3*sd_bio14
  range_bio14_sd_max<-mean_bio14+3*sd_bio14
  bio14<-nb_v[between(bio14, range_bio14_sd_min, range_bio14_sd_max)]$bio14
  range_bio14_sd_min<-min(bio14)
  range_bio14_sd_max<-max(bio14)
  
  range_bio14_IQR_min<-mean_bio14-1.5*IQR_bio14
  range_bio14_IQR_max<-mean_bio14+1.5*IQR_bio14
  
  min_x<-min(nb_v$x, na.rm=T)
  max_x<-max(nb_v$x, na.rm=T)
  mean_x<-mean(nb_v$x, na.rm=T)
  
  min_y<-min(nb_v$y, na.rm=T)
  max_y<-max(nb_v$y, na.rm=T)
  max_abs_y<-max(abs(nb_v$y), na.rm=T)
  mean_y<-mean(nb_v$y, na.rm=T)
  
  
  N_CELL<-nrow(all_items_n)
  fit<-data.frame(
    mean_bio1=mean_bio1,
    sd_bio1=sd_bio1,
    quantile_bio1_low=quantile_bio1[1],
    quantile_bio1_high=quantile_bio1[2],
    IQR_bio1=IQR_bio1,
    range_bio1_sd_min=range_bio1_sd_min,
    range_bio1_sd_max=range_bio1_sd_max,
    range_bio1_IQR_min=range_bio1_IQR_min,
    range_bio1_IQR_max=range_bio1_IQR_max,
    
    mean_bio5=mean_bio5,
    sd_bio5=sd_bio5,
    quantile_bio5_low=quantile_bio5[1],
    quantile_bio5_high=quantile_bio5[2],
    IQR_bio5=IQR_bio5,
    range_bio5_sd_min=range_bio5_sd_min,
    range_bio5_sd_max=range_bio5_sd_max,
    range_bio5_IQR_min=range_bio5_IQR_min,
    range_bio5_IQR_max=range_bio5_IQR_max,
    
    mean_bio6=mean_bio6,
    sd_bio6=sd_bio6,
    quantile_bio6_low=quantile_bio6[1],
    quantile_bio6_high=quantile_bio6[2],
    IQR_bio6=IQR_bio6,
    range_bio6_sd_min=range_bio6_sd_min,
    range_bio6_sd_max=range_bio6_sd_max,
    range_bio6_IQR_min=range_bio6_IQR_min,
    range_bio6_IQR_max=range_bio6_IQR_max,
    
    mean_bio12=mean_bio12,
    sd_bio12=sd_bio12,
    quantile_bio12_low=quantile_bio12[1],
    quantile_bio12_high=quantile_bio12[2],
    IQR_bio12=IQR_bio12,
    range_bio12_sd_min=range_bio12_sd_min,
    range_bio12_sd_max=range_bio12_sd_max,
    range_bio12_IQR_min=range_bio12_IQR_min,
    range_bio12_IQR_max=range_bio12_IQR_max,
    
    mean_bio13=mean_bio13,
    sd_bio13=sd_bio13,
    quantile_bio13_low=quantile_bio13[1],
    quantile_bio13_high=quantile_bio13[2],
    IQR_bio13=IQR_bio13,
    range_bio13_sd_min=range_bio13_sd_min,
    range_bio13_sd_max=range_bio13_sd_max,
    range_bio13_IQR_min=range_bio13_IQR_min,
    range_bio13_IQR_max=range_bio13_IQR_max,
    
    mean_bio14=mean_bio14,
    sd_bio14=sd_bio14,
    quantile_bio14_low=quantile_bio14[1],
    quantile_bio14_high=quantile_bio14[2],
    IQR_bio14=IQR_bio14,
    range_bio14_sd_min=range_bio14_sd_min,
    range_bio14_sd_max=range_bio14_sd_max,
    range_bio14_IQR_min=range_bio14_IQR_min,
    range_bio14_IQR_max=range_bio14_IQR_max,
    min_x=min_x,
    max_x=max_x,
    mean_x=mean_x,
    min_y=min_y,
    max_y=max_y,
    mean_y=mean_y,
    max_abs_y=max_abs_y,
    N_CELL=N_CELL
  )
  
  
  print("saving result")
  saveRDS(fit, target)
  
  
  
  #saveRDS(p_buffer, sprintf("%s/p_buffer.rda", target_folder))
}
