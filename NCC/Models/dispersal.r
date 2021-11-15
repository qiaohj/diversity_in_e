library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
library(rmapshaper)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Birds"
}
exposure_threshold<-0
dispersal<-1
exposure_threshold<-as.numeric(args[2])
if (is.na(exposure_threshold)){
  exposure_threshold<-5
}
dispersal<-as.numeric(args[3])
if (is.na(dispersal)){
  dispersal<-1
}
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
print("reading env layers")
mask_100km<-raster("../../Raster/mask_100km.tif")
mask_points<-data.frame(rasterToPoints(mask_100km))
mask_points = st_as_sf(mask_points, coords = c("x", "y"), crs = st_crs(mask_100km))

start_env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df_100km.rda")
future_env_layers<-readRDS("../../Objects/stacked_layers_2021_2100_list_100km.rda")

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
density<-density[density<=0.5]
density<-density*2
#hist(density)
get_disp_dist<-function(n, max_disp, density){
  v<-density * max_disp
  v[sample(length(density), n, replace = T)]
}

predict_range<-c(2021:2100)

i=2
unique<-unique[sample(length(unique), length(unique))]
x_size<-dim(mask_100km)[2]
bi<-"Lophornis pavoninus"

group_full_sum_area<-group_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
group_full_sum_area<-group_full_sum_area[order(-1*sum_are),]


bi<-group_full_sum_area[group_full_sum_area$sum_are<=1.5*min(group_full_sum_area$sum_are)]$SCINAME[1]
group_full_sum_area<-group_full_sum_area[sample(nrow(group_full_sum_area), nrow(group_full_sum_area))]
i=1
bi="Poicephalus rufiventris"
for (i in 1:length(group_full_sum_area$SCINAME)) {
  bi<-group_full_sum_area$SCINAME[i]
  if (bi=="Aratinga maculata"){
    asdf
  }
  print(paste(i, length(group_full_sum_area$SCINAME), bi, 
              "exposure", exposure_threshold, "dispersal", dispersal))
  
  target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
  fit_str<-sprintf("%s/fit.rda", target_folder)
  if (!file.exists(fit_str)){
    next()
  }
  check_point<- sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda", 
                        target_folder, exposure_threshold, dispersal)
  if (file.exists(check_point)){
    next()
  }
  saveRDS(NULL, check_point)
  
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
    disp_dist<-group_disp[iucn_name==bi]$estimated_disp
    buffer<-disp_dist * 80
    tmp_sf_buff<-st_buffer(tmp_sf, ceiling(buffer) * 1000)
    print("cut by buffer")
    mask_buffer<-crop(mask_100km, tmp_sf_buff)
    mask_buffer<-mask(mask_buffer, tmp_sf_buff)
    print("buffer to point")
    p_buffer<-data.table(rasterToPoints(mask_buffer))
  }else{
    disp_dist<-0
  }
  item_str<-"UKESM1_SSP245"
  tmp_sf_b<-as_Spatial(st_buffer(tmp_sf, 50000))
  mask_nb<-crop(mask_100km, tmp_sf_b)
  mask_nb<-mask(mask_nb, tmp_sf_b)
  
  
  initial_disp<-data.table(rasterToPoints(mask_nb))
  print("saving potential distribution result")
  saveRDS(initial_disp, check_point)
  
  
  
  first_disp<-initial_disp
  
  for (item_str in names(future_env_layers)){
    print(item_str)
    item<-future_env_layers[[item_str]]
    if (dispersal!=0){
      item<-item[mask_100km %in% p_buffer$mask_100km]
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
        print(paste(i, length(unique), 
                  bi, "exposure", exposure_threshold, 
                  "dispersal", dispersal, year_i, item_str))
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
        prev_dis[mask_100km %in% env_item$mask_100km]$suitable<-1
        prev_dis[suitable==1]$exposure<-0
        prev_dis[suitable==0]$exposure<-prev_dis[suitable==0]$exposure + 1
        prev_dis<-prev_dis[exposure<=exposure_threshold]
        if (nrow(prev_dis)>0){
          prev_dis$YEAR<-year_i
          dispersal_log[[as.character(year_i)]]<-prev_dis
          selected_cols<-c("x", "y", "mask_100km", "exposure", "suitable", "accumulative_disp", "disp")
          prev_dis<-unique(prev_dis[, ..selected_cols])
        }
      }else{
        prev_dis$suitable<-0
        prev_dis[mask_100km %in% env_item$mask_100km]$suitable<-1
        prev_dis[suitable==1]$exposure<-0
        prev_dis$disp<-get_disp_dist(nrow(prev_dis), max_dispersal * 1000, density)
        prev_dis$accumulative_disp<-prev_dis$disp + prev_dis$accumulative_disp
        
        moveable_dis<-prev_dis[suitable==1]
        if (nrow(moveable_dis)>0){
          edge_points_list<-moveable_dis[, is_edge(mask_100km, moveable_dis$mask_100km, x_size), 
                                         by = 1:nrow(moveable_dis)]
          edge_points<-moveable_dis[which(edge_points_list$V1),]
          moveable_dis$is_edge<-F
          edge_points$is_edge<-T
        }else{
          edge_points<-data.frame()
        }
        edge_points<-rbindlist(list(edge_points, 
                                    moveable_dis[!(mask_100km %in% edge_points$mask_100km)]))
        
        if (nrow(edge_points)!=0){
          edge_points[(is_edge==F)&(accumulative_disp<50000*sqrt(2))]$accumulative_disp<-50000*sqrt(2)
          pts     <- sf::st_as_sf(edge_points, coords = c("x", "y"), remove = F, crs=crs(mask_buffer))
          pts_buf <- sf::st_buffer(pts, edge_points$accumulative_disp)
          pts_buf_union<-st_cast(st_union(pts_buf), "POLYGON")
          potential_area<-mask_points[Reduce(c, st_intersects(pts_buf_union, mask_points)),]
          new_item<-data.table(st_coordinates(potential_area))
          colnames(new_item)<-c("x", "y")
          new_item$mask_100km<-potential_area$mask_100km
          new_item<-new_item[!(mask_100km %in% prev_dis$mask_100km)]
          new_item$exposure<-0
          new_item$suitable<-1
          new_item$accumulative_disp<-0
          new_item$disp<-0
          new_item<-new_item[mask_100km %in% env_item$mask_100km]
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
          selected_cols<-c("x", "y", "mask_100km", "exposure", "suitable", "accumulative_disp", "disp")
          prev_dis<-unique(prev_dis[, ..selected_cols])
        }
        prev_dis[prev_dis$accumulative_disp>100000]$accumulative_disp<-100000
      }
      
    }
    
    print("Writing result")
    saveRDS(dispersal_log, sprintf("%s/%s_%d_dispersal_%d.rda", target_folder, item_str,
                                   exposure_threshold, dispersal))
    print("Done! Writing result")
  }
}

