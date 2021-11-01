library(ceramic)
library(raster)
library(data.table)
library(igraph)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(tidyr)
library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
rm(list=ls())
source("commonFuns/colors.r")
source("commonFuns/functions.r")
plot.new()
setDTthreads(threads=1)
print(sprintf("%d CPUs are using", getDTthreads()))
GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
if (is.na(group)){
  group<-"Mammals"
}
exposure<-as.numeric(args[2])
if (is.na(exposure)){
  exposure<-5
}

index<-as.numeric(args[3])
if (is.na(index)){
  index<-1
}


layer_df<-expand.grid(SSP=SSPs, GCM=GCMs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1



j=1
mask<-raster("../../Raster/mask_index.tif")
r_continent<-raster("../../Raster/Continent_ect4.tif")
sp_i<-1
#exposure<-5
l_i<-1
#group<-"Mammals"
rm<-c()
#for (l_i in c(1:nrow(layer_df))){
for (l_i in c(index)){
  layer_item<-layer_df[l_i,]
  for (exposure in c(exposure)){
    
    for (group in c(group)){
      target_rda<-sprintf("../../Objects/cluster_based_pathway/merged/%s_%s_exposure_%d.rda",
                          group, layer_item$LABEL, exposure)
      if (file.exists(target_rda)){
        next()
      }
      df_sp_list<-data.table(readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group)))
      
      
      ttt<-2
      #df_sp_list<-df_sp_list[sum_are>ttt]
      df_sp_list$SP<-gsub(" ", "_", df_sp_list$SP)
      smooth_path_full<-list()
      for (sp_i in c(1:nrow(df_sp_list))){
        
        sp_test<-df_sp_list[sp_i,]
        print(sprintf("combination:%d/%d %s, exposure:%d, sp:%d/%d, %s, %s", 
                      l_i, nrow(layer_df), layer_item$LABEL, exposure, sp_i, nrow(df_sp_list),
                      sp_test$SP, group))
        target_folder<-sprintf("../../Objects/cluster_based_pathway/%s/%s/%s/exposure_%d", 
                               group, sp_test$SP, layer_item$LABEL, exposure)
        #smooth_path_rds<-sprintf("%s/smoothed_pathway.rda", target_folder)
        path_rds<-sprintf("%s/pathways.rda", target_folder)
        #if ((file.exists(smooth_path_rds))&(file.exists(path_rds))){
        if (file.exists(path_rds)){
          #next()
          path<-readRDS(path_rds)
          path$from_label<-path$from
          path$to_label<-path$to
          path<-path %>%
            separate(from_label, c("from_YEAR", "from_mask_index", "from_c_i"), "_")
          path<-path %>%
            separate(to_label, c("to_YEAR", "to_mask_index", "to_c_i"), "_")
          path$from_YEAR<-as.numeric(path$from_YEAR)
          path$to_YEAR<-as.numeric(path$to_YEAR)
          path$max_dist<-path$to_YEAR-path$from_YEAR
          path<-path[max_dist>=min_dist]
          path_df<-path
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
          keyyears<-unique(c(path_df$from_YEAR, path_df$to_YEAR))
          
          g<-graph_from_data_frame(path_df, directed = T, vertices=keycenters)
          
          
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
          smooth_path<-rbindlist(sm_path_full)
          if (nrow(smooth_path)==0){
            next()
          }
          smooth_path$alpha<-((smooth_path$YEAR-2020)/80)^5
          #saveRDS(smooth_path, sprintf("%s/smoothed_pathway.rda", target_folder))
          
          smooth_path$group<-group
          smooth_path$sp<-sp_test$SP
          smooth_path$GCM<-layer_item$GCM
          smooth_path$SSP<-layer_item$SSP
          smooth_path_full[[sp_test$SP]]<-smooth_path
        }
        else{
          rm<-c(rm, sprintf("rm -rf %s", target_folder))
          
        }
      }
      smooth_path_full<-rbindlist(smooth_path_full)
      saveRDS(smooth_path_full, target_rda)
    }
  }
}

#write.table(rm, file="~/xxx.sh", quote=F, row.names = F, col.names = F)
