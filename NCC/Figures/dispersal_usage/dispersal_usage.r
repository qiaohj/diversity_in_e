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

source("commonFuns/functions.r")
predict_range<-c(2021:2100)
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-5
}
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

ttt<-2
df_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", group))
df_list<-df_list[which(df_list$area>ttt),]
i=3369
#dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
dispersals<-c(1)
#df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
final_df<-NULL
#beginCluster()
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  #if ((i<10)|(item$sp=="Physalaemus_evangelistai")){
    
  #}else{
  #  next()
  #}
  if (item$area<=0){
    next()
  }
  
  target_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", group, item$sp)
  print(paste(i, nrow(df_list), item$sp, target_folder))
  start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
  model<-readRDS(sprintf("%s/fit.rda", target_folder))
  #start_dis<-start_dis%>%dplyr::filter(in_out==1)
  start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y, mask_index)
  if (is.null(start_dis)){
    next()
  }
  if (nrow(start_dis)==0){
    next()
  }
  #colnames(start_dis)<-c("x", "y")
  
  target<-sprintf("%s/dispersal_%d", target_folder, threshold)
  
  j=7
  for (j in c(1:nrow(layer_df))){
    layer_item<-layer_df[j,]
    k=1
    for (k in c(1:length(dispersals))){
      dispersal<-dispersals[k]
      print(sprintf("%s/%s_%d.rda", target, layer_item$LABEL, dispersal, threshold))
      dispersal_log<-readRDS(sprintf("%s/%s_%d.rda", target, layer_item$LABEL, dispersal))
      if (is.null(dispersal_log)){
        next()
      }
      
      
      dispersal_log_end<-dispersal_log%>%dplyr::filter(YEAR==2100)
      if (F){
        plot(dispersal_log_end$x, dispersal_log_end$y)
        range(dispersal_log_end$PR)
        range(dispersal_log_end[which(dispersal_log_end$exposure==0),]$PR)
        range(dispersal_log_end[which(dispersal_log_end$exposure==0),]$TEMP)
        model[, c("range_PR_sd_min", "range_PR_sd_max", "range_TEMP_sd_min", "range_TEMP_sd_max")]
        table(dispersal_log_end$exposure)
        points(start_dis$x, start_dis$y, col="red")
      }
      if (is.null(dispersal_log_end)){
        next()
      }
      if (nrow(dispersal_log_end)==0){
        next()
      }
      
      dist_ma<-start_dis%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, dispersal_log_end)/100000)
      min_d<-min(dist_ma$dist)
      if (min_d<=1){
        next()
      }

      if (!("mask_index" %in% names(start_dis))){
        print("EXTINCT INDEX")
        start_dis$mask_index<-extract(mask, start_dis)
      }
      dispersal_log_others<-dispersal_log%>%
        dplyr::filter((YEAR!=2100)&!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))%>%
        dplyr::distinct(x, y, mask_index, YEAR)
      if (is.null(dispersal_log_others)){
        next()
      }
      if (nrow(dispersal_log_others)==0){
        next()
      }
      
      
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
      dispersal_log_others$dispersal<-dispersal
      if (nrow(dispersal_log_others)==0){
        next()
      }
      if (is.null(final_df)){
        final_df<-dispersal_log_others
      }else{
        final_df<-full_join(final_df, dispersal_log_others, by=c("mask_index", "YEAR", "GCM", "SSP", "dispersal"))
        
        final_df[is.na(final_df)]<-0
        final_df$N_SP<-final_df$N_SP.x+final_df$N_SP.y
        final_df<-final_df[, c("mask_index", "YEAR", "GCM", "SSP", "N_SP", "dispersal")]
        print(nrow(final_df))
      }
      
    }
  }
}
saveRDS(final_df, sprintf("../../Figures_Full_species/dispersal_usage_%d/%s_ttt_%d.rda", 
                          threshold, group, ttt))
#endCluster()

if (F){
  df1<-readRDS("../../Objects/Niche_Models/Amphibians/Leptopelis_zebra/dispersal/EC-Earth3-Veg_SSP585_1.rda")
  df2<-readRDS("../../Objects/Niche_Models/Amphibians/Leptopelis_zebra/dispersal/UKESM1_SSP245_2.rda")
  }
