library(ceramic)
library(raster)
library(data.table)
library(igraph)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(tidyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
rm(list=ls())
source("commonFuns/colors.r")
source("commonFuns/functions.r")
plot.new()
setDTthreads(threads=1)
print(sprintf("%d CPUs are using", getDTthreads()))
get_min_cluster<-function(N){
  Ncluster<-round(N/100)
  if (Ncluster<1){
    Ncluster<-1
  }
  if (Ncluster>5){
    Ncluster<-5
  }
  Ncluster
}
get_max_cluster<-function(N){
  Ncluster<-round(N/10)
  if (Ncluster<2){
    Ncluster<-2
  }
  if (Ncluster>10){
    Ncluster<-10
  }
  Ncluster
}

GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1
df_sp_list<-list()
group<-"Mammals"
for (group in c("Mammals", "Birds")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
  df_list$group<-group
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list, fill=T)

ttt<-2
#df_sp_list<-df_sp_list[area>ttt]
df_sp_list$SP<-gsub(" ", "_", df_sp_list$SP)
df_sp_list<-df_sp_list[sample(nrow(df_sp_list), nrow(df_sp_list)),]
j=1
mask<-raster("../../Raster/mask_100km.tif")
r_continent<-raster("../../Raster/Continent_ect4.tif")
sp_i<-412
exposure<-5
l_i<-1

points_10km<-readRDS("../../Raster/points_10km.rda")
colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-data.table(rasterToPoints(mask_100km))
colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
points_10km<-merge(points_10km, points_100km, by="mask_100km")

for (l_i in c(1:nrow(layer_df))){
  layer_item<-layer_df[l_i,]
  for (exposure in c(0, 5)){
    for (sp_i in c(1:nrow(df_sp_list))){
      
      sp_test<-df_sp_list[sp_i,]
      print(sprintf("combination:%d/%d %s, exposure:%d, sp:%d/%d, %s, %s", 
                    l_i, nrow(layer_df), layer_item$LABEL, exposure, sp_i, nrow(df_sp_list),
                    sp_test$SP, sp_test$group))
      target_folder<-sprintf("../../Objects/cluster_based_pathway_10km/%s/%s/%s/exposure_%d", 
                      sp_test$group, sp_test$SP, layer_item$LABEL, exposure)
      
      if (dir.exists(target_folder)){
        print("skip")
        next()
      }
      dir.create(target_folder, showWarnings = F, recursive = T)
      source_folder<-sprintf("../../Objects/Dispersal/%s/%s", sp_test$group, sp_test$SP)
      if (!file.exists(sprintf("%s/%s_%d_dispersal_1_10km.rda", 
                               source_folder, layer_item$LABEL, exposure))){
        next()
      }
      dis_details_10km<-readRDS(sprintf("%s/%s_%d_dispersal_1_10km.rda", 
                                   source_folder, layer_item$LABEL, exposure))
      dis_details_10km<-rbindlist(dis_details_10km)
      if (is.null(dis_details_10km)){
        next()
      }
      if (nrow(dis_details_10km)==0){
        next()
      }
      dis_details<-merge(dis_details_10km, points_10km, by.x="mask_10km", by.y="mask_10km")
      if (nrow(dis_details)==0){
        next()
      }
      dis_details<-dis_details[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR", "exposure", "suitable")]
      
      colnames(dis_details)[c(1:2)]<-c("x", "y")
      
      dis_details$continent<-raster::extract(r_continent, dis_details[, c("x", "y")])
      dis_details<-dis_details[!is.na(continent)]
      if (nrow(dis_details)==0){
        next()
      }
      keyyears<-unique(dis_details$YEAR)
      keyyears<-keyyears[1:length(keyyears)]
      if ((length(keyyears)>6)){
        keyyears<-unique(round(quantile(keyyears, seq(0, 1, 0.2))))
      }
      if (length(keyyears)==1){
        next()
      }
      g_df<-list()
      for (year in keyyears){
        current_dis<-dis_details[between(YEAR, year-1, year+1)]
        current_dis<-unique(current_dis[, c("x", "y", "mask_100km", "continent")])

        c_i=1
        current_dis$cluster<-0
        for(c_i in unique(current_dis$continent)){
          current_dis_item<-current_dis[continent==c_i]
          N_row<-nrow(current_dis_item)
          if (N_row>2){
            pamk.best <- pamk(current_dis_item[, c("x", "y")], 
                              krange=get_min_cluster(N_row):get_max_cluster(N_row))
            current_dis[continent==c_i]$cluster<-pamk.best$pamobject$clustering
            centers_pam<-data.frame(pamk.best$pamobject$medoids)
          }else{
            current_dis[continent==c_i]$cluster<-1
            centers_pam<-current_dis_item[1, c("x", "y")]
          }
          
          centers_pam$YEAR<-year
          centers_pam$continent<-c_i
          
          g_df[[sprintf("%d_%d", year, c_i)]]<-centers_pam
        }
      }
      g_df<-rbindlist(g_df)
      g_df$mask_100km<-raster::extract(mask, g_df[, c("x", "y")])
      saveRDS(g_df, sprintf("%s/key_centers.rda", target_folder))
      
      
      
      path_df<-NULL
      for (c_i in unique(g_df$continent)){
        g_df_sub<-g_df[continent==c_i]
        for (year_i in c(1:(length(keyyears)-1))){
          centers_1<-g_df_sub[YEAR==keyyears[year_i]]
          centers_2<-g_df_sub[YEAR==keyyears[year_i+1]]
          if ((nrow(centers_1)==0)|(nrow(centers_2)==0)){
            next()
          }
          centers<-g_df_sub[YEAR %in% keyyears[year_i:(year_i+1)]]
          i<-1
          for (i in c(1:nrow(centers))){
            center<-centers[i,]
            if (center$YEAR==keyyears[year_i]){
              from_label<-sprintf("%d_%d_%d", keyyears[year_i], center$mask_100km, center$continent)
              
              if (from_label %in% path_df$from){
                next()
              }
              distances<-sqrt((center$x-centers_2$x)^2+(center$y-centers_2$y)^2)/100000
              next_centers<-centers_2[distances==min(distances)]
              to_label<-sprintf("%d_%d_%d", keyyears[year_i+1], next_centers$mask_100km, next_centers$continent)
              item_df<-data.table(from=from_label,
                                  to=to_label, 
                                  from_x=center$x, from_y=center$y,
                                  to_x=next_centers$x, to_y=next_centers$y, 
                                  YEAR=keyyears[year_i+1], min_dist=min(distances),
                                  continent=c_i)
              path_df<-bind(path_df, item_df)
            }else{
              to_label<-sprintf("%d_%d_%d", keyyears[year_i+1], center$mask_100km, center$continent)
              
              if (to_label %in% path_df$to){
                next()
              }
              distances<-sqrt((center$x-centers_1$x)^2+(center$y-centers_1$y)^2)/100000
              next_centers<-centers_1[distances==min(distances)]
              from_label<-sprintf("%d_%d_%d", keyyears[year_i], next_centers$mask_100km, next_centers$continent)
              item_df<-data.table(from=from_label,
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
      if (nrow(path_df)==1){
        next()
      }
      path_df$TYPE<-"MIDDLE"
      path_df[YEAR==min(path_df$YEAR)]$TYPE<-"START"
      path_df[(TYPE!="START")&(!(to %in% path_df$from))]$TYPE<-"END"
      saveRDS(path_df, sprintf("%s/pathways.rda", target_folder))
      
      from_v<-path_df[, c("from", "from_x", "from_y", "continent")]
      colnames(from_v)<-c("label", "x", "y", "continent")
      to_v<-path_df[, c("to", "to_x", "to_y", "continent")]
      colnames(to_v)<-c("label", "x", "y", "continent")
      
      
      keycenters<-unique(rbindlist(list(from_v, to_v)))
      keycenters$name<-keycenters$label
      keycenters<-keycenters %>%
        separate(label, c("YEAR", "mask_100km", "c_i"), "_")
      col_order<-rev(colnames(keycenters))
      keycenters<-keycenters[,..col_order]
      
      g<-graph_from_data_frame(path_df, directed = T, vertices=keycenters)
      saveRDS(g, sprintf("%s/pathway_graph.rda", target_folder))
      
      from_labels<-unique(path_df[TYPE=="START"]$from)
      to_labels<-unique(path_df[TYPE=="END"]$to)
      if ((length(from_labels)==0)|(length(to_labels)==0)){
        next()
      }
      sm_path_full<-list()
      path_group<-1
      path_df$path_group<-0
      for (f_i in c(1:length(from_labels))){
        from_label<-rev(from_labels[f_i])
        for (to_label in to_labels){
          x<-get.shortest.paths(g, from_label, to_label)$vpath[[1]]
          if (length(x)>0){
            path<-path_df[(from %in% names(x))&(to %in% names(x))]
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
      if (nrow(sm_path_full_df)==0){
        next()
      }
      sm_path_full_df$alpha<-((sm_path_full_df$YEAR-2020)/80)^5
      saveRDS(sm_path_full_df, sprintf("%s/smoothed_pathway.rda", target_folder))
      if (F){
        mask_black<-raster("../../Raster/mask.tif")
        mask_p<-data.frame(rasterToPoints(mask_black))
        start_dis<-readRDS(sprintf("%s/occ_with_env.rda", source_folder))
        dispersal_based_end<-dis_details[YEAR==2100]
        p_bak<-ggplot() + 
          geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
          map_theme
        p<-p_bak+
          geom_tile(data=dispersal_based_end, aes(x=x, y=y), fill="black", alpha=0.4)+
          geom_path(data=sm_path_full_df, aes(x=x, y=y, group=path_group, alpha=alpha), color="red")+
          geom_tile(data=start_dis, aes(x=x, y=y), fill="blue", alpha=0.4)+
          #geom_path(data=path_df, aes(x=from_x, y=from_y, color=YEAR,group=path_group))+
          
          scale_alpha_continuous()
        p
        ggsave(p, filename=sprintf("%s/fig.pdf", target_folder))
      }
    }
  }
}

