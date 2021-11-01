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
#for BIRDS
if (F){
  r<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/EC-Earth3-Veg/SSP119/1850/bio12_eck4_1km.tif")
  plot(r)
  continent<-readOGR("../../Shape/continents", "continent")
  crs(continent)<-CRS("+proj=longlat +datum=WGS84")
  continent_eck4<-spTransform(continent, crs(r))
  continent_eck4<-continent_eck4[which(continent_eck4$CONTINENT!="Antarctica"),]
  writeOGR(continent_eck4, "../../Shape/continents", "continent_eck4", driver="ESRI Shapefile")
  r_c<-crop(r, extent(continent_eck4))
  r_c<-mask(r_c, continent_eck4)
  plot(r_c)
  v<-values(r_c)
  
  values(r_c)[!is.na(v)]<-c(1:length(v[!is.na(v)]))
  plot(r_c)
  writeRaster(r_c, "../../Raster/mask_1km.tif", datatype="INT4U")
}

if (F){
  mask_bak<-raster("../../Raster/mask_10km.tif")
  mask_high<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/mask_high_res.tif")
  vessel <- sf::st_read(dsn = "/home/huijieqiao/Experiments/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
  #vessel_simpl <- st_simplify(vessel,
  #                            preserveTopology = TRUE, 
  #                            dTolerance = 10000)
  #object.size(vessel_simpl)
  #object.size(vessel)
  
  
  #plot_map(mv_simpl)
  
  sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_bak))
  #st_write(sp_df_eck4, "../../Data/Raw/IUCN/BIRDS/BIRDS_ECK4.gdb")
  bird_df<-data.frame(sp_df_eck4)
  bird_df$Shape<-NULL
  
  colnames(bird_df)
  bird_df<-data.table(bird_df)
  saveRDS(bird_df, "../../Data/Birds/bird_df.rda")
}
mask_10km<-raster("../../Raster/mask_10km.tif")
#mask_1km<-raster("../../Raster/mask_1km.tif")

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2))
}

PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)

if (F){
  vessel <- sf::st_read(dsn = "/home/huijieqiao/Experiments/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
  sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_10km))
  
  bird_df<-readRDS("../../Data/Birds/bird_df.rda")
  unique <- unique(bird_df$SCINAME)
  i=10
  for (i in 1:length(unique)) {
    bi<-unique[i]
    print(paste(i, length(unique), bi))
    tmp_sf <- sp_df_eck4[which((bird_df$SCINAME == bi)&(bird_df$PRESENCE %in% PRESENCE)&
                                 (bird_df$ORIGIN %in% ORIGIN)&(bird_df$SEASONAL %in% SEASONAL)), ]
    saveRDS(tmp_sf, sprintf("../../Objects/IUCN_Distribution/Birds/RAW/%s.rda", gsub(" ", "_", bi)))
    tryCatch(
      {
        tmp_sf_sim<-st_simplify(tmp_sf, dTolerance = 5000)
        saveRDS(tmp_sf_sim, sprintf("../../Objects/IUCN_Distribution/Birds/st_simplify/%s.rda", gsub(" ", "_", bi)))
        tmp_sf_rm_sim<-rmapshaper::ms_simplify(tmp_sf, keep=0.01)
        saveRDS(tmp_sf_rm_sim, sprintf("../../Objects/IUCN_Distribution/Birds/ms_simplify/%s.rda", gsub(" ", "_", bi)))
      },
      error=function(cond) {
        
      }
    )
    
  }
}
bird_df<-readRDS("../../Data/Birds/bird_df.rda")
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
bird_full<-merge(bird_df, bird_disp, by.x="SCINAME", by.y="iucn_name", all=F)

unique <- unique(bird_full$SCINAME)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
i=511

if (F){
  stat<-NULL
  for (i in 1:length(unique)) {
    bi<-unique[i]
    print(paste(i, length(unique), bi))
    tmp_sf <- sp_df_eck4[which((bird_df$SCINAME == bi)&(bird_df$PRESENCE %in% PRESENCE)&
                                 (bird_df$ORIGIN %in% ORIGIN)&(bird_df$SEASONAL %in% SEASONAL)), ]
    if (class(tmp_sf$Shape)[1]=="sfc_MULTISURFACE"){
      asdfasdf
    }
    plot(st_geometry(tmp_sf))
    
    tmp_sf %>% st_collection_extract("MULTIPOLYGON")
    st_write(tmp_sf, sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)), append=FALSE)
    poly<-st_read(sprintf("../../Data/Raw/IUCN/BIRDS/Weird_Polygons/%s.shp", gsub(" ", "_", bi)))
    head(st_as_text(tmp_sf$Shape[1]))
    st_cast(st_sfc(tmp_sf), "MULTIPOLYGON")
    st_cast(st_sf(a = 1, st_sfc(tmp_sf$Shape)), "MULTIPOLYGON")
    
    sum_area<-sum(tmp_sf$Shape_Area)
    n<-nrow(tmp_sf)
    type<-class(tmp_sf$Shape)[1]
    item<-data.frame(sp=bi, sum_area=sum_area, n=n, type=type)
    if (is.null(stat)){
      stat<-item
    }else{
      stat<-bind_rows(stat, item)
    }
    unique(stat$type)
    stat%>%filter(type=="sfc_GEOMETRY")
    stat%>%filter(type=="sfc_MULTISURFACE")
    stat%>%filter(type=="sfc_MULTIPOLYGON")
    mask<-raster(extent(tmp_sf), res=res(mask_bak), crs=crs(mask_bak))
    st_crs(tmp_sf)<-crs(mask_bak)
    extract(mask_bak, tmp_sf)
  }
}

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
bi<-"Tangara rufigenis"
bird_full_sum_area<-bird_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
bird_full_sum_area<-bird_full_sum_area[order(sum_are),]
bi<-bird_full_sum_area[bird_full_sum_area$sum_are<=1.5*min(bird_full_sum_area$sum_are)]$SCINAME[1]
bird_full_sum_area<-bird_full_sum_area[sample(nrow(bird_full_sum_area), nrow(bird_full_sum_area))]
for (i in 1:length(bird_full_sum_area$SCINAME)) {
  
  bi<-bird_full_sum_area$SCINAME[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  target_folder<-sprintf("../../Objects/Dispersal/Birds/%s", gsub(" ", "_", bi))
  if (dir.exists(target_folder)){
    next()
  }
  dir.create(target_folder)
  tmp_sf<-NULL
  print("extracting the matched polygons")
  if (file.exists(sprintf("../../Objects/IUCN_Distribution/Birds/st_simplify/%s.rda", gsub(" ", "_", bi)))){
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/Birds/st_simplify/%s.rda", gsub(" ", "_", bi)))
  }else{
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/Birds/RAW/%s.rda", gsub(" ", "_", bi)))
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
