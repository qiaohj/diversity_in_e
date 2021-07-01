
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
p_t<-p

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
      
      
      p_t<-p_t+#geom_path(data=path, aes(x=from_x, y=from_y), color="pink")+
        geom_path(data=sm_path, aes(x=x, y=y), color="red")
      
    }
    
  }
  
}

p_t

ggsave(p_t, filename="../../Figures/test.png", width=15, height=10)