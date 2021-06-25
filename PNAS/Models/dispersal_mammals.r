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

mammal_df<-readRDS("../../Data/Mammals/mammal_df.rda")
mammal_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
mammal_full<-merge(mammal_df, mammal_disp, by.x="binomial", by.y="Scientific", all=F)

unique <- unique(mammal_full$binomial)
unique<-as.character(unique)
PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)
print("reading env layers")

mask_100km<-raster("../../Raster/mask_100km.tif")
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
get_disp_dist<-function(n, max_disp){
  disp_seed<-rexp(n = n, rate = .1)
  disp_seed/max(disp_seed) * max_disp
}

predict_range<-c(2021:2100)
exposure_threshold<-5
i=2
unique<-unique[sample(length(unique), length(unique))]
x_size<-dim(mask_100km)[2]
bi<-"Poicephalus rufiventris"
mammal_full_sum_area<-mammal_full[, .(sum_are=sum(SHAPE_Area)), by="binomial"]
mammal_full_sum_area<-mammal_full_sum_area[order(sum_are),]
if (F){
  mammal_full_sum_area<-mammal_full[, .(sum_are=sum(SHAPE_Area)), 
                                    by=list(binomial, ForStrat, log_body_mass, Diet, estimated_disp)]
  i=1
  final_df<-NULL
  for (i in 1:length(mammal_full_sum_area$binomial)) {
    bi<-mammal_full_sum_area$binomial[i]
    print(paste(i, length(mammal_full_sum_area$binomial), bi))
    target_folder<-sprintf("../../Objects/Mammals/%s", gsub(" ", "_", bi))
    fit_str<-sprintf("%s/fit.rda", target_folder)
    if (!file.exists(fit_str)){
      next()
    }
    N_file<-length(list.files(target_folder))
    if (N_file==11){
      item<-mammal_full_sum_area[i,]
      
      fit<-readRDS(fit_str)
      
      item<-cbind(item, fit)
      colnames(item)[1]<-"SP"
      if (is.null(final_df)){
        final_df<-item
      }else{
        final_df<-rbindlist(list(final_df, item))
      }
    }
  }
  saveRDS(final_df, "../../Objects/IUCN_List/Mammals_df.rda")
  
}
bi<-mammal_full_sum_area[mammal_full_sum_area$sum_are<=1.5*min(mammal_full_sum_area$sum_are)]$binomial[1]
mammal_full_sum_area<-mammal_full_sum_area[sample(nrow(mammal_full_sum_area), nrow(mammal_full_sum_area))]

for (i in 1:length(mammal_full_sum_area$binomial)) {
  bi<-mammal_full_sum_area$binomial[i]
  print(paste(i, length(unique), bi))
  target_folder<-sprintf("../../Objects/Mammals/%s", gsub(" ", "_", bi))
  fit_str<-sprintf("%s/fit.rda", target_folder)
  if (!file.exists(fit_str)){
    next()
  }
  
  check_point<- sprintf("%s/initial_disp_exposure_%d.rda", target_folder, exposure_threshold)
  if (file.exists(check_point)){
    next()
  }
  saveRDS(NULL, check_point)
  
  fit<-readRDS(fit_str)
  
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
  
  
  print("calculating buffer")
  disp_dist<-mammal_disp[Scientific==bi]$estimated_disp
  buffer<-disp_dist * 80
  tmp_sf_buff<-st_buffer(tmp_sf, ceiling(buffer) * 1000)
  print("cut by buffer")
  mask_buffer<-crop(mask_100km, tmp_sf_buff)
  mask_buffer<-mask(mask_buffer, tmp_sf_buff)
  print("buffer to point")
  p_buffer<-data.table(rasterToPoints(mask_buffer))
  item_str<- names(future_env_layers)[9]
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
    item<-item[mask_100km %in% p_buffer$mask_100km]
    item<-item[between(bio1, fit$range_bio1_sd_min, fit$range_bio1_sd_max)&
                 between(bio5, fit$range_bio5_sd_min, fit$range_bio5_sd_max)&
                 between(bio6, fit$range_bio6_sd_min, fit$range_bio6_sd_max)&
                 between(bio12, fit$range_bio12_sd_min, fit$range_bio12_sd_max)&
                 between(bio13, fit$range_bio13_sd_min, fit$range_bio13_sd_max)&
                 between(bio14, fit$range_bio14_sd_min, fit$range_bio14_sd_max)]
    
    prev_dis<-first_disp
    dispersal_log<-list()
    year_i<-predict_range[1]
    max_dispersal<-disp_dist
    prev_dis$exposure<-0
    prev_dis$suitable<-1
    prev_dis$accumulative_disp<-0
    year_i = 2021
    for (year_i in predict_range){
      
      print(paste("", i, length(unique), bi, year_i, item_str))
      if (nrow(prev_dis)==0){
        next()
      }
      
      prev_dis$accumulative_disp<-get_disp_dist(nrow(prev_dis), max_dispersal * 1000) + prev_dis$accumulative_disp
      prev_dis[suitable==0]$accumulative_disp<-0
      moveable_dis<-prev_dis[suitable==1]
      if (nrow(moveable_dis)>0){
        edge_points_list<-moveable_dis[, is_edge(mask_100km, moveable_dis$mask_100km, x_size), 
                                       by = 1:nrow(moveable_dis)]
        edge_points<-moveable_dis[which(edge_points_list$V1),]
      }else{
        edge_points<-data.frame()
      }
      
      if (nrow(edge_points)!=0){
        if (F){
          plot(prev_dis$x, prev_dis$y)
          points(edge_points$x, edge_points$y, col="red")
        }
        
        range_x<-range(prev_dis$x)
        range_x<-c(range_x[1]-1500*max_dispersal, range_x[2]+1500*max_dispersal)
        range_y<-range(prev_dis$y)
        range_y<-c(range_y[1]-1500*max_dispersal, range_y[2]+1500*max_dispersal)
        
        
        env_item<-item[year==year_i]
        env_item<-env_item[(x %between% range_x)&(y %between% range_y)]
        
        pts     <- sf::st_as_sf(edge_points, coords = c("x", "y"), remove = F, crs=crs(mask_buffer))
        pts_buf <- sf::st_buffer(pts, edge_points$accumulative_disp)
        
        pts_buf_union<-st_cast(st_union(pts_buf), "POLYGON")
        pts_buf_union_f<-as_Spatial(st_buffer(pts_buf_union, 50000))
        potential_area<-crop(mask_100km, pts_buf_union_f)
        potential_area<-mask(potential_area, pts_buf_union_f)
        
        
        if (F){
          st_write(pts, "../../Figures/Example/point.shp")
          st_write(pts_buf, "../../Figures/Example/point_buffer.shp")
          st_write(pts_buf_union, "../../Figures/Example/point_buffer_union.shp")
          plot(st_geometry(pts))
          plot(st_geometry(pts_buf), add=T)
          plot(st_geometry(pts_buf_union), add=T, col="blue")
          plot(potential_area, add=T, col="red")
          plot(st_geometry(pts_buf_union), add=T, col="blue")
        }
        
        new_item<-data.table(rasterToPoints(potential_area))
        new_item<-new_item[!(mask_100km %in% prev_dis$mask_100km)]
        new_item$exposure<-0
        new_item$suitable<-1
        new_item$accumulative_disp<-0
        new_item<-new_item[mask_100km %in% env_item$mask_100km]
      }else{
        new_item<-NULL
      }
      prev_dis$suitable<-0
      prev_dis[mask_100km %in% env_item$mask_100km]$suitable<-1
      prev_dis[suitable==1]$exposure<-0
      prev_dis[suitable==0]$exposure<-prev_dis[suitable==0]$exposure + 1
      prev_dis<-prev_dis[exposure<=exposure_threshold]
      if (!is.null(new_item)){
        prev_dis<-rbindlist(list(prev_dis, new_item))
      }
      if (F){
        plot(prev_dis$x, prev_dis$y, col=prev_dis$suitable, pch=".")
        points(new_item$x, new_item$y, col="red")
      }
      
      
      if (nrow(prev_dis)>0){
        prev_dis$YEAR<-year_i
        dispersal_log[[as.character(year_i)]]<-prev_dis
        selected_cols<-c("x", "y", "mask_100km", "exposure", "suitable", "accumulative_disp")
        prev_dis<-unique(prev_dis[, ..selected_cols])
      }
      #prev_dis<-prev_dis[exposure==0]
      prev_dis[prev_dis$accumulative_disp>50000]$accumulative_disp<-0
    }
    if (F){
      xxxx=2021
      library(ggplot2)
      all<-rbindlist(dispersal_log)
      x_size<-range(all$x)
      y_size<-range(all$y)
      
      for (xxxx in predict_range){
        fff<-dispersal_log[[as.character(xxxx)]]
        p<-ggplot(fff)+geom_tile(aes(x=x, y=y))+xlim(x_size)+ylim(y_size)+ggtitle(xxxx)+theme_bw()
        ggsave(p, filename=sprintf("../../Figures/Example/distributions/%d_xxx.png", xxxx))
      }
    }
    print("Writing result")
    saveRDS(dispersal_log, sprintf("%s/%s_%d.rda", target_folder, item_str, exposure_threshold))
    print("Done! Writing result")
  }
}
