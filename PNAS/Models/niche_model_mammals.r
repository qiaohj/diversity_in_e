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
setDTthreads(3)
print(sprintf("Current core number is %d", getDTthreads()))
#for MAMMALS

mask_10km<-raster("../../Raster/mask_10km.tif")
#mask_1km<-raster("../../Raster/mask_1km.tif")

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2))
}

PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)


mammal_df<-readRDS("../../Data/Mammals/mammal_df.rda")
mammal_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
mammal_full<-merge(mammal_df, mammal_disp, by.x="binomial", by.y="Scientific", all=F)

unique <- unique(mammal_full$binomial)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
i=511

print("reading env layers")

start_env_layers_se<-readRDS("../../Objects/stacked_layers_1850_2020_df_10km.rda")
#future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list.rda")

is_edge<-function(index, all_index, xsize){
  #print(index)
  neighbors<-c(index-1, index+1, index-xsize, index+xsize)
  if (sum(neighbors %in% all_index * 1)==4){
    F
  }else{
    T
  }
}


predict_range<-c(2021:2100)
exposure_threshold<-0
i=2
unique<-unique[sample(length(unique), length(unique))]
x_size<-dim(mask_10km)[2]

mammal_full_sum_area<-mammal_full[, .(sum_are=sum(SHAPE_Area)), by="binomial"]
mammal_full_sum_area<-mammal_full_sum_area[order(sum_are),]
bi<-mammal_full_sum_area[mammal_full_sum_area$sum_are<=1.5*min(mammal_full_sum_area$sum_are)]$binomial[1]
mammal_full_sum_area<-mammal_full_sum_area[sample(nrow(mammal_full_sum_area), nrow(mammal_full_sum_area))]
for (i in 1:length(mammal_full_sum_area$binomial)) {
  
  bi<-mammal_full_sum_area$binomial[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target_folder<-sprintf("../../Objects/Mammals/%s", gsub(" ", "_", bi))
  if (dir.exists(target_folder)){
    next()
  }
  dir.create(target_folder)
  tmp_sf<-NULL
  print("extracting the matched polygons")
  if (file.exists(sprintf("../../Objects/IUCN_Distribution/Mammals/st_simplify/%s.rda", gsub(" ", "_", bi)))){
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/Mammals/st_simplify/%s.rda", gsub(" ", "_", bi)))
  }else{
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/Mammals/RAW/%s.rda", gsub(" ", "_", bi)))
  }
  
  if ((class(tmp_sf$Shape)[1]=="sfc_GEOMETRY")|(class(tmp_sf$Shape)[1]=="sfc_MULTISURFACE")){
    next()
  }
  
  tmp_sf<-tmp_sf[!st_is_empty(tmp_sf),]
  if (nrow(tmp_sf)==0){
    next()
  }
  
  
  print("calculating niche breadth")
  mask_nb<-crop(mask_10km, tmp_sf)
  mask_nb<-mask(mask_nb, tmp_sf)
  
  
  nb_p<-data.table(rasterToPoints(mask_nb))
  if (nrow(nb_p)==0){
    next()
  }
  nb_v<-start_env_layers_se[mask_10km %in% nb_p$mask_10km]
  initial_disp<-nb_v[year==2020]
  cols<-c("x", "y", "mask_10km")
  first_disp<-unique(initial_disp[, ..cols])
  mean_bio1<-mean(nb_v$bio1, na.rm=T)
  sd_bio1<-sd(nb_v$bio1, na.rm=T)
  quantile_bio1<-quantile(nb_v$bio1, c(0.25, 0.75), na.rm=T)
  IQR_bio1<-quantile_bio1[2]-quantile_bio1[1]
  range_bio1_sd_min<-mean_bio1-3*sd_bio1
  range_bio1_sd_max<-mean_bio1+3*sd_bio1
  range_bio1_IQR_min<-mean_bio1-1.5*IQR_bio1
  range_bio1_IQR_max<-mean_bio1+1.5*IQR_bio1
  
  mean_bio5<-mean(nb_v$bio5, na.rm=T)
  sd_bio5<-sd(nb_v$bio5, na.rm=T)
  quantile_bio5<-quantile(nb_v$bio5, c(0.25, 0.75), na.rm=T)
  IQR_bio5<-quantile_bio5[2]-quantile_bio5[1]
  range_bio5_sd_min<-mean_bio5-3*sd_bio5
  range_bio5_sd_max<-mean_bio5+3*sd_bio5
  range_bio5_IQR_min<-mean_bio5-1.5*IQR_bio5
  range_bio5_IQR_max<-mean_bio5+1.5*IQR_bio5
  
  mean_bio6<-mean(nb_v$bio6, na.rm=T)
  sd_bio6<-sd(nb_v$bio6, na.rm=T)
  quantile_bio6<-quantile(nb_v$bio6, c(0.25, 0.75), na.rm=T)
  IQR_bio6<-quantile_bio6[2]-quantile_bio6[1]
  range_bio6_sd_min<-mean_bio6-3*sd_bio6
  range_bio6_sd_max<-mean_bio6+3*sd_bio6
  range_bio6_IQR_min<-mean_bio6-1.5*IQR_bio6
  range_bio6_IQR_max<-mean_bio6+1.5*IQR_bio6
  
  mean_bio12<-mean(nb_v$bio12, na.rm=T)
  sd_bio12<-sd(nb_v$bio12, na.rm=T)
  quantile_bio12<-quantile(nb_v$bio12, c(0.25, 0.75), na.rm=T)
  IQR_bio12<-quantile_bio12[2]-quantile_bio12[1]
  range_bio12_sd_min<-mean_bio12-3*sd_bio12
  range_bio12_sd_max<-mean_bio12+3*sd_bio12
  range_bio12_IQR_min<-mean_bio12-1.5*IQR_bio12
  range_bio12_IQR_max<-mean_bio12+1.5*IQR_bio12
  
  mean_bio13<-mean(nb_v$bio13, na.rm=T)
  sd_bio13<-sd(nb_v$bio13, na.rm=T)
  quantile_bio13<-quantile(nb_v$bio13, c(0.25, 0.75), na.rm=T)
  IQR_bio13<-quantile_bio13[2]-quantile_bio13[1]
  range_bio13_sd_min<-mean_bio13-3*sd_bio13
  range_bio13_sd_max<-mean_bio13+3*sd_bio13
  range_bio13_IQR_min<-mean_bio13-1.5*IQR_bio13
  range_bio13_IQR_max<-mean_bio13+1.5*IQR_bio13
  
  mean_bio14<-mean(nb_v$bio14, na.rm=T)
  sd_bio14<-sd(nb_v$bio14, na.rm=T)
  quantile_bio14<-quantile(nb_v$bio14, c(0.25, 0.75), na.rm=T)
  IQR_bio14<-quantile_bio14[2]-quantile_bio14[1]
  range_bio14_sd_min<-mean_bio14-3*sd_bio14
  range_bio14_sd_max<-mean_bio14+3*sd_bio14
  range_bio14_IQR_min<-mean_bio14-1.5*IQR_bio14
  range_bio14_IQR_max<-mean_bio14+1.5*IQR_bio14
  
  min_x<-min(nb_v$x, na.rm=T)
  max_x<-max(nb_v$x, na.rm=T)
  mean_x<-mean(nb_v$x, na.rm=T)
  
  min_y<-min(nb_v$y, na.rm=T)
  max_y<-max(nb_v$y, na.rm=T)
  max_abs_y<-max(abs(nb_v$y), na.rm=T)
  mean_y<-mean(nb_v$y, na.rm=T)
  
  
  N_CELL<-nrow(nb_v[year==1850])
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
  saveRDS(fit, sprintf("%s/fit.rda", target_folder))
  
  
  
  #saveRDS(p_buffer, sprintf("%s/p_buffer.rda", target_folder))
}

if (F){
  #plot(st_geometry(tmp_sf), col="red")
  plot(st_geometry(tmp_sf_buff))
  plot(st_geometry(tmp_sf), add=T, col="red")
  plot(mask_buffer, add=T)
  
  points(env_item$x, env_item$y, pch=".", col="blue")
  points(prev_dis$x, prev_dis$y, pch=".", col="red")
  plot(mask_nb, col="red", add=T)
  
  plot(dispersal_log$`2021`$x, dispersal_log$`2021`$y, pch=".")
  points(dispersal_log$`2022`$x, dispersal_log$`2022`$y, col="red", pch=".")
  points(dispersal_log$`2021`$x, dispersal_log$`2021`$y, pch=".")
  
}
