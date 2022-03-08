library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Birds"
}


env<-NULL
esm_i<-as.numeric(args[2])
if (is.na(esm_i)){
  esm_i<-floor(runif(1) * 9 + 1)
}
esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
esm_ssp_str<-esm_ssp[esm_i]
esm_ssp_str_split<-strsplit(esm_ssp_str, "_")[[1]]
esm<-esm_ssp_str_split[1]
ssp<-esm_ssp_str_split[2]
print(group)

if (group=="Birds"){
  group_df<-readRDS("../../Data/Birds/bird_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  #group_disp$estimated_disp
  #group_disp2<-readRDS("../../Objects_PNAS/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  #dd<-merge(data.frame(group_disp[, c("iucn_name", "estimated_disp")]), 
  #          data.frame(group_disp2[, c("iucn_name", "estimated_disp")]), 
  #          by.x="iucn_name", by.y="iucn_name", all.x=T, all.y=F)
  group_full<-merge(group_df, group_disp, by.x="SCINAME", by.y="iucn_name", all=F)
  unique <- unique(group_full$SCINAME)
  unique<-as.character(unique)
  
}else{
  group_df<-readRDS("../../Data/Mammals/mammal_df.rda")
  group_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  group_full<-merge(group_df, group_disp, by.x="binomial", by.y="Scientific", all=F)
  
  unique <- unique(group_full$binomial)
  unique<-as.character(unique)
  colnames(group_full)[1]<-"SCINAME"
  colnames(group_full)[27]<-"Shape_Area"
  colnames(group_disp)[1]<-"iucn_name"
  
}
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
print("reading mask layers")
mask_10km<-raster("../../Raster/mask_10km.tif")
mask_points<-data.frame(rasterToPoints(mask_10km))
mask_points = st_as_sf(mask_points, coords = c("x", "y"), crs = st_crs(mask_10km))



is_edge<-function(index, all_index, xsize){
  #print(index)
  neighbors<-c(index-1, index+1, index-xsize, index+xsize)
  if (sum(neighbors %in% all_index * 1)==4){
    F
  }else{
    T
  }
}
density<-rexp(n = 10000, rate = 0.1)
density<-density/max(density)

#hist(density)
get_disp_dist<-function(n, max_disp, density){
  v<-density * max_disp
  v[sample(length(density), n, replace = T)]
}

mean(get_disp_dist(10, 97640, density))
predict_range<-c(2021:2100)

i=2
unique<-unique[sample(length(unique), length(unique))]
x_size<-dim(mask_10km)[2]
bi<-"Lophornis pavoninus"

sp_list<-readRDS("../../Objects/extincted_species_list.rda")
sp_list<-sp_list[group.x==group&GCM==esm&SSP==ssp]
sp_list<-sp_list[sample(nrow(sp_list), nrow(sp_list))]

for (i in 1:nrow(sp_list)) {
  start_time<-Sys.time()
  print(paste("Start time:", start_time))
  item<-sp_list[i]
  dispersal<-item$M
  exposure_threshold<-item$exposure
  bi<-item$Scientific
  
  print(paste("20220220", i, length(sp_list$Scientific), bi, 
              "exposure", exposure_threshold, "dispersal", dispersal,
              "Dispersal distance", sp_list[i]$estimated_disp,
              "Distribution area", sp_list[i]$N))
  
  target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
  fit_str<-sprintf("%s/fit.rda", target_folder)
  if (!file.exists(fit_str)){
    next()
  }
  
  target<-sprintf("%s/%s_%d_dispersal_%d_10km.rda", target_folder, esm_ssp_str,
                  exposure_threshold, dispersal)
  if (file.exists(target)){
    print("skip")
    next()
  }
  saveRDS(NULL, target)
  check_point<- sprintf("%s/initial_disp_10km_exposure_%d_dispersal_%d.rda", 
                        target_folder, exposure_threshold, dispersal)
  initial_disp<-NULL
  if (file.exists(check_point)){
    initial_disp<-readRDS(check_point)
    
  }
  fit<-readRDS(fit_str)
  
  tmp_sf<-NULL
  print("extracting the matched polygons")
  if (file.exists(sprintf("../../Objects/IUCN_Distribution/%s/st_simplify/%s.rda", 
                          group, gsub(" ", "_", bi)))){
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/st_simplify/%s.rda", 
                            group, gsub(" ", "_", bi)))
  }else{
    tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", 
                            group, gsub(" ", "_", bi)))
  }
  
  if ((class(tmp_sf$Shape)[1]=="sfc_GEOMETRY")|(class(tmp_sf$Shape)[1]=="sfc_MULTISURFACE")){
    next()
  }
  
  tmp_sf<-tmp_sf[!st_is_empty(tmp_sf),]
  if (nrow(tmp_sf)==0){
    next()
  }
  
  if (dispersal!=0){
    print("calculating buffer")
    disp_dist<-min(group_disp[iucn_name==bi]$estimated_disp)
    buffer<-disp_dist * 80
    tmp_sf_buff<-st_buffer(tmp_sf, ceiling(buffer) * 1000)
    print("cut by buffer")
    mask_buffer<-crop(mask_10km, tmp_sf_buff)
    mask_buffer<-mask(mask_buffer, tmp_sf_buff)
    print("buffer to point")
    p_buffer<-data.table(rasterToPoints(mask_buffer))
  }else{
    disp_dist<-0
  }
  tmp_sf_b<-as_Spatial(st_buffer(tmp_sf, 5000))
  mask_nb<-crop(mask_10km, tmp_sf_b)
  mask_nb<-mask(mask_nb, tmp_sf_b)
  
  
  initial_disp<-data.table(rasterToPoints(mask_nb))
  print("saving potential distribution result")
  saveRDS(initial_disp, check_point)
  if (F){
    plot(mask_nb)
    dddd<-readRDS( sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda",
                           target_folder, exposure_threshold, dispersal))
    points(dddd$x, dddd$y)
    plot(dddd$x, dddd$y)
  }
  
  
  first_disp<-initial_disp
  if (is.null(env)){
    print(paste("Reading", esm_ssp_str))
    env<-readRDS(sprintf("../../Objects/stacked_layers_2021_2100_list_10km_%s.rda", esm_ssp_str))
  }
  item<-env
  if (dispersal!=0){
    item<-item[mask_10km %in% p_buffer$mask_10km]
  }
  item<-item[between(bio1, fit$range_bio1_sd_min, fit$range_bio1_sd_max)&
               between(bio5, fit$range_bio5_sd_min, fit$range_bio5_sd_max)&
               between(bio6, fit$range_bio6_sd_min, fit$range_bio6_sd_max)&
               between(bio12, fit$range_bio12_sd_min, fit$range_bio12_sd_max)&
               between(bio13, fit$range_bio13_sd_min, fit$range_bio13_sd_max)&
               between(bio14, fit$range_bio14_sd_min, fit$range_bio14_sd_max)]
  ###Simulation start here
  prev_dis<-first_disp
  dispersal_log<-list()
  year_i<-predict_range[1]
  max_dispersal<-disp_dist
  prev_dis$exposure<-0
  prev_dis$suitable<-1
  prev_dis$accumulative_disp<-0
  prev_dis$disp<-0
  year_i = 2021
  for (year_i in predict_range){
    if (dispersal==1){
      current_time<-Sys.time()
      print(paste("20220222 (start time", start_time, ", current time", 
                  current_time,
                  "):", i, length(sp_list$Scientific), 
                  bi, "exposure", exposure_threshold, 
                  "dispersal", dispersal, year_i, esm_ssp_str))
    }
    if (nrow(prev_dis)==0){
      next()
    }
    
    range_x<-range(prev_dis$x)
    range_x<-c(range_x[1]-1500*max_dispersal, range_x[2]+1500*max_dispersal)
    range_y<-range(prev_dis$y)
    range_y<-c(range_y[1]-1500*max_dispersal, range_y[2]+1500*max_dispersal)
    
    
    env_item<-item[year==year_i]
    env_item<-env_item[(x %between% range_x)&(y %between% range_y)]
    
    if (dispersal==0){
      prev_dis$accumulative_disp<-0
      prev_dis$disp<-0
      prev_dis$suitable<-0
      prev_dis[mask_10km %in% env_item$mask_10km]$suitable<-1
      prev_dis[suitable==1]$exposure<-0
      prev_dis[suitable==0]$exposure<-prev_dis[suitable==0]$exposure + 1
      prev_dis<-prev_dis[exposure<=exposure_threshold]
      if (nrow(prev_dis)>0){
        prev_dis$YEAR<-year_i
        dispersal_log[[as.character(year_i)]]<-prev_dis
        selected_cols<-c("x", "y", "mask_10km", "exposure", "suitable", "accumulative_disp", "disp")
        prev_dis<-unique(prev_dis[, ..selected_cols])
      }
    }else{
      prev_dis$suitable<-0
      prev_dis[mask_10km %in% env_item$mask_10km]$suitable<-1
      prev_dis[suitable==1]$exposure<-0
      prev_dis$disp<-get_disp_dist(nrow(prev_dis), max_dispersal * 1000, density)
      prev_dis$accumulative_disp<-prev_dis$disp + prev_dis$accumulative_disp
      
      moveable_dis<-prev_dis[suitable==1]
      if (nrow(moveable_dis)>0){
        edge_points_list<-moveable_dis[, is_edge(mask_10km, moveable_dis$mask_10km, x_size), 
                                       by = 1:nrow(moveable_dis)]
        edge_points<-moveable_dis[which(edge_points_list$V1),]
        moveable_dis$is_edge<-F
        edge_points$is_edge<-T
      }else{
        edge_points<-data.frame()
      }
      edge_points<-rbindlist(list(edge_points, 
                                  moveable_dis[!(mask_10km %in% edge_points$mask_10km)]))
      
      if (nrow(edge_points)!=0){
        edge_points[(is_edge==F)&(accumulative_disp<5000*sqrt(2))]$accumulative_disp<-5000*sqrt(2)
        pts     <- sf::st_as_sf(edge_points, coords = c("x", "y"), remove = F, crs=crs(mask_buffer))
        pts_buf <- sf::st_buffer(pts, edge_points$accumulative_disp)
        pts_buf_union<-st_cast(st_union(pts_buf), "POLYGON")
        potential_area<-mask_points[Reduce(c, st_intersects(pts_buf_union, mask_points)),]
        new_item<-data.table(st_coordinates(potential_area))
        colnames(new_item)<-c("x", "y")
        new_item$mask_10km<-potential_area$mask_10km
        new_item<-new_item[!(mask_10km %in% prev_dis$mask_10km)]
        new_item$exposure<-0
        new_item$suitable<-1
        new_item$accumulative_disp<-0
        new_item$disp<-0
        new_item<-new_item[mask_10km %in% env_item$mask_10km]
      }else{
        new_item<-NULL
      }
      
      prev_dis[suitable==0]$exposure<-prev_dis[suitable==0]$exposure + 1
      prev_dis<-prev_dis[exposure<=exposure_threshold]
      if (!is.null(new_item)){
        if (nrow(new_item)>0){
          prev_dis<-rbindlist(list(prev_dis, new_item))
        }
      }
      if (nrow(prev_dis)>0){
        prev_dis$YEAR<-year_i
        dispersal_log[[as.character(year_i)]]<-prev_dis
        selected_cols<-c("x", "y", "mask_10km", "exposure", "suitable", "accumulative_disp", "disp")
        prev_dis<-unique(prev_dis[, ..selected_cols])
      }
      prev_dis[prev_dis$accumulative_disp>10000]$accumulative_disp<-10000
    }
    
  }
  
  print("Writing result")
  saveRDS(dispersal_log, target)
  print("Done! Writing result")
  
  end_time<-Sys.time()
  print(paste("end time", end_time))
}

