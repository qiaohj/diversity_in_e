library(batchtools)
library(rgl)
library(ceramic)
library(anglr)
library(ggnewscale)
library(ggplot2)
library(raster)
library(data.table)
library(polylabelr)
library(data.tree)
library(igraph)
library(dbscan)

plot_path <- function(x, y, ...) {
  plot.new()
  plot.window(range(x, na.rm = TRUE), range(y, na.rm = TRUE))
  polypath(x, y, ...)
}
x <- c(5, 10, 10, 5, 5, 6, 6, 7, 7, 6, 8, 8, 9, 9, 8)
y <- c(5, 5, 10, 10, 5, 6, 7, 7, 6, 6, 8, 9, 9, 8, 8)
plot_path(x, y, col = "grey", border = NA)
points(poi(x, y))
## Not run:
# Find visual centers for North Carolina counties
library(sf)
nc <- st_read(system.file("shape/nc.shp", package="sf"))
locations = do.call(rbind, poi(nc, precision=0.01))
plot(st_geometry(nc))
points(locations)
## End(Not run)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
alt<-raster("../../Raster/ALT/alt_eck4_high_res.tif")
mask<-raster("../../Raster/mask.tif")
mask_p<-data.frame(rasterToPoints(mask))
source("commonFuns/colors.r")
source("commonFuns/functions.r")
#p2<-data.frame(rasterToPoints(alt))
#p2[which(p2$x<=-12103059), "alt_eck4_high_res"]<-NA
#p2[which((p2$x>12912000)&(p2$y>5000000)), "alt_eck4_high_res"]<-NA
#values(alt)[!is.na(values(alt))]<-p2$alt_eck4_high_res
if (F){
  slope = terrain(alt, opt='slope')
  aspect = terrain(alt, opt='aspect')
  hill = hillShade(slope, aspect)
  dem_spdf <- as(alt, "SpatialPixelsDataFrame")
  dem_spdf <- as.data.frame(dem_spdf)
  colnames(dem_spdf) <- c("value", "x", "y")
  
  hill_spdf <- as(hill, "SpatialPixelsDataFrame")
  hill_spdf <- as.data.frame(hill_spdf)
  colnames(hill_spdf) <- c("value", "x", "y")
}
p_bak<-ggplot() + 
  geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
  map_theme

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1
df_sp_list<-list()
for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list)

ttt<-2
df_sp_list<-df_sp_list[area>ttt]
df_sp_list$sp<-gsub(" ", "_", df_sp_list$sp)
j=9
euc.dist <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2) ^ 2+(y1-y2)^2)
}
threshold=1
start_range<-c(600000, 2200000)
end_range<-c(4000000, 7000000)
layer_item<-layer_df[j,]
smooth_path<-readRDS(sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s.rda", 
                             threshold, layer_item$LABEL))
print(sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s.rda", threshold, layer_item$LABEL))

smooth_path<-smooth_path[sp %in% df_sp_list$sp]
smooth_path$line_group<-paste(smooth_path$sp, smooth_path$continent_i)
smooth_path$index<-raster::extract(mask, smooth_path[, c("x", "y")])
smooth_path<-smooth_path[!is.na(index)]
smooth_path$pre_x<--1
smooth_path[2:nrow(smooth_path), "pre_x"]<-smooth_path[1:(nrow(smooth_path)-1), "x"]
smooth_path$pre_y<--1
smooth_path[2:nrow(smooth_path), "pre_y"]<-smooth_path[1:(nrow(smooth_path)-1), "y"]
smooth_path$dist<-euc.dist(smooth_path$x, smooth_path$y, smooth_path$pre_x, smooth_path$pre_y)
smooth_path[1, "dist"]<-0
smooth_path[YEAR==2020]$dist<-0
smooth_path<-smooth_path[dist<100000]

smooth_path$alpha<-((smooth_path$YEAR-2020)/80)^5
smooth_path_start<-smooth_path[(between(x, start_range[1], start_range[2]))&(alpha==0)&(survive=="SURVIVE")]
smooth_path_start<-smooth_path_start[y>=5000000]
smooth_path_end<-smooth_path[(between(x, end_range[1], end_range[2]))&(alpha==1)&(survive=="SURVIVE")]
smooth_path_end<-smooth_path_end[y>=5000000]
target_sp<-smooth_path_start[sp %in% smooth_path_end$sp]
target_sp<-target_sp[group=="Mammals"]
smooth_path_target<-smooth_path[sp %in% target_sp$sp]
smooth_path_target<-smooth_path_target[continent_i==2]
p<-p_bak+geom_path(data=smooth_path_target, aes(x=x, y=y, alpha=alpha, color=group,
                                                group=line_group))+
  scale_alpha_continuous()+
  scale_color_manual(values = color_groups)
p
#Miniopterus_schreibersii
sp_test<-target_sp[2, ]
path_test<-smooth_path_target[sp==sp_test$sp]

target_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", sp_test$group, sp_test$sp)
start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
p<-p_bak+geom_tile(data=start_dis, aes(x=x, y=y), fill="grey50")+
  geom_path(data=path_test, aes(x=x, y=y, alpha=alpha, color=group,
                                         group=line_group))+
  scale_alpha_continuous()+
  scale_color_manual(values = color_groups)
p


#burn in
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")

dispersal<-1
mask<-raster("../../Raster/mask_index.tif")
start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
start_dis$mask_index<-extract(mask, data.frame(start_dis[, c("x", "y")]))
selected_cols<-c("x", "y", "mask_index")
start_dis<-unique(start_dis[, ..selected_cols])
start_dis$exposure<-0

predict_range<-c(2021:2100)
prev_size<-nrow(prev_dis)
new_size=0

#potential_dis

model<-readRDS(sprintf("%s/fit.rda", target_folder))
env_item_start<-data.table(env_layers[[paste(layer_item$LABEL, 2021, sep="_")]])
potential_dis_start<-env_item_start[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                      (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                      (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]

env_item_end<-data.table(env_layers[[paste(layer_item$LABEL, 2100, sep="_")]])
potential_dis_end<-env_item_end[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                                 (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                                 (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]

dis_details<-readRDS(sprintf("%s/dispersal_%d/%s_%d.rda", target_folder, threshold, layer_item$LABEL, dispersal))

dispersal_based_end<-dis_details[YEAR==2100]
start_point_mean<-data.table(x=mean(start_dis$x), y=mean(start_dis$y))
start_point_poi<-poi(start_dis$x, start_dis$y, precision=0.01)
start_point_poi<-data.table(x=start_point_poi$x, y=start_point_poi$y)

start_point_dbscan<-dbscan(start_dis[, c("x", "y")], minPts = 3, eps=110000)
start_dis$cluster<-start_point_dbscan$cluster
start_dis_filter<-start_dis[cluster>0]
start_point_poi<-start_dis_filter%>%dplyr::group_by(cluster)%>%
  dplyr::summarise(c_x=poi(x, y)$x,
                   c_y=poi(x, y)$y)
start_point_pam<-start_dis_filter%>%dplyr::group_by(cluster)%>%
  dplyr::summarise(c_x=cluster::pam(data.frame(x, y), k=1)$medoids[,1],
                   c_y=cluster::pam(data.frame(x, y), k=1)$medoids[,2])



if (F){
  
  ggplot(start_dis)+geom_tile(aes(x=x, y=y, fill=factor(cluster)))+
    geom_point(data=start_point_pam, aes(x=c_x, y=c_y))
}
end_point_mean<-data.table(x=mean(dispersal_based_end$x), y=mean(dispersal_based_end$y))
end_point_poi<-poi(dispersal_based_end$x, dispersal_based_end$y, precision=0.01)
end_point_poi<-data.table(x=end_point_poi$x, y=end_point_poi$y)
end_point_pam<-data.frame(cluster::pam(dispersal_based_end[, c("x", "y")], k=5)$medoids)



p<-p_bak+
  #geom_tile(data=potential_dis_start, aes(x=x, y=y), fill="lightblue", alpha=0.6)+
  #geom_tile(data=potential_dis_end, aes(x=x, y=y), fill="red")+
  geom_tile(data=dispersal_based_end, aes(x=x, y=y), fill="grey50")+
  geom_tile(data=start_dis, aes(x=x, y=y), fill="black", alpha=0.4)+
  geom_path(data=path_test, aes(x=x, y=y, alpha=alpha, color=group,
                                group=line_group))+
  geom_point(data=start_point_mean, aes(x=x, y=y), color="blue", shape=3)+
  geom_point(data=start_point_poi, aes(x=x, y=y), color="blue")+
  geom_point(data=start_point_pam, aes(x=x, y=y), color="blue", shape=2)+
  
  geom_point(data=end_point_mean, aes(x=x, y=y), color="red", shape=3)+
  geom_point(data=end_point_poi, aes(x=x, y=y), color="red")+
  geom_point(data=end_point_pam, aes(x=x, y=y), color="red", shape=2)+
  
  scale_alpha_continuous()+
  scale_color_manual(values = color_groups)
p
dispersal_log<-NULL

year<-2021
i=1

g_df<-list()
prev_dis<-start_dis
prev_dis$YEAR<-2020
for (year in predict_range){
  #for (year in c(2021, 2022, 2023)){
  print(year)
  current_dis<-dis_details[YEAR==year]
  new_dis<-current_dis[!(mask_index %in% prev_dis$mask_index)]
  g_df_item<-NULL
  for (i in c(1:nrow(new_dis))){
    item<-new_dis[i,]
    distance<-sqrt((item$x-prev_dis$x)^2+(item$y-prev_dis$y)^2)/100000
    prev_dis_cells<-prev_dis[round(distance-1)==0]
    prev_dis_cells<-prev_dis_cells[YEAR==max(prev_dis_cells$YEAR)]
    if (nrow(prev_dis_cells)>0){
      
      item_df<-data.frame(from=sprintf("%d_%d", year-1, prev_dis_cells$mask_index),
                          to=sprintf("%d_%d", year, item$mask_index), 
                          from_x=prev_dis_cells$x, from_y=prev_dis_cells$y,
                          to_x=item$x, to_y=item$y, 
                          YEAR=year)
      g_df_item<-bind(g_df_item, item_df)
    }
  }
  g_df[[as.character(year)]]<-g_df_item
  prev_dis<-dis_details[YEAR<=year]
}
g_df<-rbindlist(g_df)
g_df$TYPE<-"MIDDLE"
g_df[YEAR==2021]$TYPE<-"START"
g_df[!(to %in% g_df$from)]$TYPE<-"END"
vertices<-unique(g_df$from, g_df$to)
g<-graph_from_data_frame(g_df, directed = T)
g

plot(g)
E(g)
p<-p_bak+
  #geom_tile(data=potential_dis_start, aes(x=x, y=y), fill="lightblue", alpha=0.6)+
  #geom_tile(data=potential_dis_end, aes(x=x, y=y), fill="red")+
  geom_tile(data=dispersal_based_end, aes(x=x, y=y), fill="black", alpha=0.4)+
  geom_tile(data=start_dis, aes(x=x, y=y), fill="blue", alpha=0.4)+
  #geom_path(data=path_test, aes(x=x, y=y, alpha=alpha, color=group,
  #                              group=line_group))+
  #geom_point(data=start_point_mean, aes(x=x, y=y), color="blue", shape=3)+
  #geom_point(data=start_point_poi, aes(x=x, y=y), color="blue")+
  #geom_point(data=start_point_pam, aes(x=x, y=y), color="blue", shape=2)+
  
  #geom_point(data=end_point_mean, aes(x=x, y=y), color="red", shape=3)+
  #geom_point(data=end_point_poi, aes(x=x, y=y), color="red")+
  #geom_point(data=end_point_pam, aes(x=x, y=y), color="red", shape=2)+
  
  scale_alpha_continuous()+
  scale_color_manual(values = color_groups)
p

from_labels<-unique(g_df[TYPE=="START"]$from)
to_labels<-unique(g_df[TYPE=="END"]$to)
to_labels_end<-unique(g_df[YEAR==2100]$to)

sm_path_full<-list()
path_group<-1
for (f_i in c(1:length(from_labels))){

  from_label<-rev(from_labels[f_i])
  print(paste(f_i, length(from_labels), from_label))
  for (to_label in to_labels){
    #for (to_label in to_labels[1:100]){
    x<-get.shortest.paths(g, from_label, to_label)$vpath[[1]]
    if (length(x)>0){
      path<-g_df[(from %in% names(x))&(to %in% names(x))]
      
      
      
      keyyears<-unique(path$YEAR)
      keyyears<-unique(round(quantile(keyyears, seq(0, 1, 0.2))))
      all_c_few<-path%>%dplyr::filter(YEAR %in% keyyears)
      if (nrow(all_c_few)<=1){
        next()
      }
      full_path<-rbind(data.frame(x=all_c_few$from_x, y=all_c_few$from_y), 
                       data.frame(x=all_c_few[nrow(all_c_few), ]$to_x, y=all_c_few[nrow(all_c_few), ]$to_y))
      
      sm_path<-data.frame(xspline(full_path$x, full_path$y, shape=1, draw=F, repEnds=T))
      sm_path$YEAR<-seq(keyyears[1], keyyears[length(keyyears)], 
                        by=(keyyears[length(keyyears)]-keyyears[1])/(nrow(sm_path)-1))
      sm_path$path_group<-path_group
      sm_path_full[[sprintf("%s_%s_%d", from_label, to_label, path_group)]]<-sm_path
      path_group<-path_group+1
      
        
    }
    
  }
  
}
sm_path_full_df<-rbindlist(sm_path_full)
    sm_path_full_df$alpha<-((sm_path_full_df$YEAR-2020)/80)^5
p<-p_bak+
  #geom_tile(data=potential_dis_start, aes(x=x, y=y), fill="lightblue", alpha=0.6)+
  #geom_tile(data=potential_dis_end, aes(x=x, y=y), fill="red")+
  geom_tile(data=dispersal_based_end, aes(x=x, y=y), fill="black", alpha=0.4)+
  geom_path(data=sm_path_full_df, aes(x=x, y=y, group=path_group, alpha=alpha), color="red")+
  geom_tile(data=start_dis, aes(x=x, y=y), fill="blue", alpha=0.4)+
  #geom_path(data=path_test, aes(x=x, y=y, alpha=alpha, color=group,
  #                              group=line_group))+
  #geom_point(data=start_point_mean, aes(x=x, y=y), color="blue", shape=3)+
  #geom_point(data=start_point_poi, aes(x=x, y=y), color="blue")+
  #geom_point(data=start_point_pam, aes(x=x, y=y), color="blue", shape=2)+
  
  #geom_point(data=end_point_mean, aes(x=x, y=y), color="red", shape=3)+
  #geom_point(data=end_point_poi, aes(x=x, y=y), color="red")+
  #geom_point(data=end_point_pam, aes(x=x, y=y), color="red", shape=2)+
  
  scale_alpha_continuous()+
  scale_color_manual(values = color_groups)
p


#p_t<-p+geom_path(data=sm_path_full_df, aes(x=x, y=y, group=path_group), color="red")
ggsave(p, filename="../../Figures/test.png", width=15, height=10)
