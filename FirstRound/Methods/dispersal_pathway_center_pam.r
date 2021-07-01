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
library(tidyr)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)


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
sp_test<-df_sp_list[1,]
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
r_continent<-raster("../../Raster/Continent_ect4.tif")
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

dis_details<-readRDS(sprintf("%s/dispersal_%d/%s_%d.rda", 
                             target_folder, threshold, layer_item$LABEL, dispersal))

dispersal_based_end<-dis_details[YEAR==2100]
start_point_mean<-data.table(x=mean(start_dis$x), y=mean(start_dis$y))
start_point_poi<-poi(start_dis$x, start_dis$y, precision=0.01)
start_point_poi<-data.table(x=start_point_poi$x, y=start_point_poi$y)

start_point_dbscan<-dbscan(start_dis[, c("x", "y")], minPts = 3, eps=110000)
start_dis$cluster<-start_point_dbscan$cluster
start_dis_filter<-start_dis[cluster>0]
start_point_poi<-start_dis_filter%>%dplyr::group_by(cluster)%>%
  dplyr::summarise(x=poi(x, y)$x,
                   y=poi(x, y)$y)
start_point_pam<-start_dis_filter%>%dplyr::group_by(cluster)%>%
  dplyr::summarise(x=cluster::pam(data.frame(x, y), k=1)$medoids[,1],
                   y=cluster::pam(data.frame(x, y), k=1)$medoids[,2])



if (F){
  
  ggplot(start_dis)+geom_tile(aes(x=x, y=y, fill=factor(cluster)))+
    geom_point(data=start_point_pam, aes(x=x, y=y))
}
end_point_mean<-data.table(x=mean(dispersal_based_end$x), y=mean(dispersal_based_end$y))
end_point_poi<-poi(dispersal_based_end$x, dispersal_based_end$y, precision=0.01)
end_point_poi<-data.table(x=end_point_poi$x, y=end_point_poi$y)
end_point_pam<-data.frame(cluster::pam(dispersal_based_end[, c("x", "y")], k=5)$medoids)


dispersal_log<-NULL

year<-2021
i=1
get_min_cluster<-function(N){
  Ncluster<-round(N/100)
  if (Ncluster<1){
    Ncluster<-2
  }
  if (Ncluster>5){
    Ncluster<-5
  }
  Ncluster
}
prev_dis<-start_dis
prev_dis$YEAR<-2020
cluster_threshold<-9
dis_details$continent<-raster::extract(r_continent, dis_details[, c("x", "y")])
dis_details<-dis_details[!is.na(continent)]

keyyears<-unique(dis_details$YEAR)
keyyears<-keyyears[1:length(keyyears)]
if ((length(keyyears)>6)){
  keyyears<-unique(round(quantile(keyyears, seq(0, 1, 0.2))))
}
g_df<-list()
for (year in keyyears){
  #for (year in c(2021, 2022, 2023)){
  print(year)
  current_dis<-dis_details[between(YEAR, year-1, year+1)]
  current_dis<-unique(current_dis[, c("x", "y", "mask_index", "continent")])
  #point_dbscan<-dbscan(current_dis[, c("x", "y")], minPts = 3, eps=110000)
  #fviz_nbclust(current_dis[, c("x", "y")], pam, method = "gap_stat")
  #fviz_gap_stat(gap_stat)
  
  #gap_stat <- clusGap(current_dis[, c("x", "y")], FUN = kmeans, nstart = 25,
  #                    K.max = 10, B = 50)
  c_i=1
  current_dis$cluster<-0
  for(c_i in unique(current_dis$continent)){
    current_dis_item<-current_dis[continent==c_i]
    if (nrow(current_dis_item)>1){
      pamk.best <- pamk(current_dis_item[, c("x", "y")], 
                        krange=get_min_cluster(nrow(current_dis_item)):10)
      #n_cluster<-NbClust(current_dis_item[, c("x", "y")], min.nc=4, max.nc=10, 
      #        method = "centroid", index = "silhouette")
      current_dis[continent==c_i]$cluster<-pamk.best$pamobject$clustering
      centers_pam<-data.frame(pamk.best$pamobject$medoids)
    }else{
      current_dis[continent==c_i]$cluster<-1
      centers_pam<-current_dis_item[, c("x", "y")]
    }
    
    centers_pam$YEAR<-year
    centers_pam$continent<-c_i
    
    g_df[[sprintf("%d_%d", year, c_i)]]<-centers_pam
  }
  #print(n_cluster$Best.nc)
  #current_dis$cluster<-point_dbscan$cluster
  #current_dis_filter<-current_dis[cluster>0]
  #cluster_count<-current_dis_filter[, .N, by=cluster]
  #cluster_count_filter<-cluster_count[N>=cluster_threshold]
  #if (nrow(cluster_count_filter)==0){
  #  cluster_count_filter<-cluster_count[N==max(cluster_count$N)]
  #}
  #current_dis_filter<-current_dis_filter[cluster %in% cluster_count_filter$cluster]
  
  #centers_pam<-current_dis_filter%>%dplyr::group_by(cluster)%>%
  #  dplyr::summarise(x=cluster::pam(data.frame(x, y), k=1)$medoids[,1],
  #                   y=cluster::pam(data.frame(x, y), k=1)$medoids[,2])
  
  
  
  if (F){
    clusplot(pam(current_dis_item[,c("x", "y")], pamk.best$nc))
    
    
    print(ggplot(current_dis)+geom_tile(aes(x=x, y=y, fill=factor(cluster)))+
      geom_point(data=rbindlist(g_df), aes(x=x, y=y), color="black"))
    x<-readline(prompt="Enter name: ")
    if (x=="x"){
      break()
    }
  }
  
}
g_df<-rbindlist(g_df)
g_df$mask_index<-raster::extract(mask, g_df[, c("x", "y")])
if (F){
  p_bak+geom_point(data=g_df, aes(x=x, y=y, color=YEAR))
}
year<-2022


path_df<-NULL
for (c_i in unique(g_df$continent)){
  g_df_sub<-g_df[continent==c_i]
  year_i<-1
  for (year_i in c(1:(length(keyyears)-1))){
    centers_1<-g_df_sub[YEAR==keyyears[year_i]]
    centers_2<-g_df_sub[YEAR==keyyears[year_i+1]]
    centers<-g_df_sub[YEAR %in% keyyears[year_i:(year_i+1)]]
    i<-1
    for (i in c(1:nrow(centers))){
      center<-centers[i,]
      if (center$YEAR==keyyears[year_i]){
        from_label<-sprintf("%d_%d_%d", keyyears[year_i], center$mask_index, center$continent)
        
        if (from_label %in% path_df$from){
          next()
        }
        distances<-sqrt((center$x-centers_2$x)^2+(center$y-centers_2$y)^2)/100000
        next_centers<-centers_2[distances==min(distances)]
        to_label<-sprintf("%d_%d_%d", keyyears[year_i+1], next_centers$mask_index, next_centers$continent)
        item_df<-data.frame(from=from_label,
                            to=to_label, 
                            from_x=center$x, from_y=center$y,
                            to_x=next_centers$x, to_y=next_centers$y, 
                            YEAR=keyyears[year_i+1], min_dist=min(distances),
                            continent=c_i)
        path_df<-bind(path_df, item_df)
      }else{
        to_label<-sprintf("%d_%d_%d", keyyears[year_i+1], center$mask_index, center$continent)
        
        if (to_label %in% path_df$to){
          next()
        }
        distances<-sqrt((center$x-centers_1$x)^2+(center$y-centers_1$y)^2)/100000
        next_centers<-centers_1[distances==min(distances)]
        from_label<-sprintf("%d_%d_%d", keyyears[year_i], next_centers$mask_index, next_centers$continent)
        item_df<-data.frame(from=from_label,
                            to=to_label, 
                            from_x=next_centers$x, from_y=next_centers$y, 
                            to_x=center$x, to_y=center$y,
                            YEAR=keyyears[year_i+1], min_dist=min(distances),
                            continent=c_i)
        path_df<-bind(path_df, item_df)
      }
      
      
      
    }
    
    
  }
}
if (F){
  p_bak+geom_path(data=path_df, aes(x=from_x, y=from_y, color=continent))
}

hist(path_df$min_dist)
path_df$TYPE<-"MIDDLE"
path_df[YEAR==min(path_df$YEAR)]$TYPE<-"START"
path_df[!(to %in% path_df$from)]$TYPE<-"END"
from_v<-path_df[, c("from", "from_x", "from_y", "continent")]
colnames(from_v)<-c("label", "x", "y", "continent")
to_v<-path_df[, c("to", "to_x", "to_y", "continent")]
colnames(to_v)<-c("label", "x", "y", "continent")


keycenters<-unique(rbindlist(list(from_v, to_v)))
keycenters$name<-keycenters$label
keycenters<-keycenters %>%
  separate(label, c("YEAR", "mask_index", "c_i"), "_")
col_order<-rev(colnames(keycenters))
keycenters<-keycenters[,..col_order]
keycenters[name=="2037_9619_1"]
x<-data.frame(table(keycenters$name))
g<-graph_from_data_frame(path_df, directed = T, vertices=keycenters)

g
ve<-V(g)
ve[[1]]
plot(g, layout=layout_with_lgl)
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

from_labels<-unique(path_df[TYPE=="START"]$from)
to_labels<-unique(path_df[TYPE=="END"]$to)
to_labels_end<-unique(path_df[YEAR==2100]$to)

sm_path_full<-list()
path_group<-1
path_df$path_group<-0
for (f_i in c(1:length(from_labels))){
  
  from_label<-rev(from_labels[f_i])
  print(paste(f_i, length(from_labels), from_label))
  for (to_label in to_labels){
    #for (to_label in to_labels[1:100]){
    x<-get.shortest.paths(g, from_label, to_label)$vpath[[1]]
    print(x)
    if (length(x)>0){
      
      path<-path_df[(from %in% names(x))&(to %in% names(x))]
      
      
      
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
      path_df[(from %in% names(x))&(to %in% names(x))]$path_group<-path_group
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
  #geom_path(data=path_df, aes(x=from_x, y=from_y, color=factor(path_group),
  #                              group=path_group))
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

