library(raster)
library(dplyr)
#library(alphahull)
library(concaveman)
library(sf)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=2
#dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
dispersals<-data.frame(M=1, N=1)
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
final_df<-NULL
beginCluster()
for (i in c(1:nrow(df_list))){
  print(paste(i, nrow(df_list)))
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp)
  
  target<-sprintf("%s/dispersal", target_folder)
  model<-"Mean"
  j=1
  for (j in c(1:nrow(layer_df))){
    layer_item<-layer_df[j,]
    k=1
    for (k in c(1:nrow(dispersals))){
      dispersal<-dispersals[k,]
      
      dispersal_log<-readRDS(sprintf("%s/%s_%d_%d.rda", target, layer_item$LABEL, dispersal$M, dispersal$N))
      if (is.null(dispersal_log)){
        next()
      }
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      #start_dis<-start_dis%>%dplyr::filter(in_out==1)
      start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
      colnames(start_dis)<-c("x", "y")
      start_dis$mask_index<-extract(mask, start_dis)
      
      dispersal_log_end<-dispersal_log%>%dplyr::filter(YEAR==2100)
      dispersal_log_others<-dispersal_log%>%dplyr::filter((YEAR!=2100)&!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))%>%
        dplyr::distinct(x, y, mask_index, YEAR)
      
      #points<-st_sfc(st_point(as.matrix(start_dis[, c("x", "y")])))
      points<-as.matrix(bind_rows(start_dis[, c("x", "y")], dispersal_log_end[, c("x", "y")]))
      alpha_hull<-concaveman(points)
      dispersal_log_others$inout<-point.in.polygon(dispersal_log_others$x, dispersal_log_others$y,
                       alpha_hull[,1], alpha_hull[,2])
      dispersal_log_others<-dispersal_log_others%>%
        dplyr::filter(!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))
      
      if (F){
        ggplot()+
          geom_path(data=as.data.frame(alpha_hull), aes(x=V1, y=V2), color="grey")+
          geom_point(data=dispersal_log_others, aes(x=x, y=y, color=factor(inout)), shape=3)+
          geom_point(data=start_dis, aes(x=x, y=y), color="red", shape=8)+
          geom_point(data=dispersal_log_end, aes(x=x, y=y), color="black", shape=8)
      }
      dispersal_log_others<-dispersal_log_others%>%dplyr::filter(inout==1)
      dispersal_log_others$GCM<-layer_item$GCM
      dispersal_log_others$SSP<-layer_item$SSP
      dispersal_log_others$N_SP<-1
      if (nrow(dispersal_log_others)==0){
        next()
      }
      if (is.null(final_df)){
        final_df<-dispersal_log_others
      }else{
        final_df<-left_join(final_df, dispersal_log_others, by=c("mask_index", "YEAR", "GCM", "SSP"))
        final_df[is.na(final_df)]<-0
        final_df$N_SP<-final_df$N_SP.x+final_df$N_SP.y
        final_df<-final_df[, c("mask_index", "YEAR", "GCM", "SSP", "N_SP")]
      }
      
    }
  }
}
saveRDS(final_df, sprintf("../../Figures/dispersal_usage/%s.rda", group))
endCluster()