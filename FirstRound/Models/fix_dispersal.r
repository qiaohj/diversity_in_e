library(raster)
#library(rgdal)
#library(rgeos)
#library(MASS)
#library(cluster)
library(data.table)
library(batchtools)
#library(ggplot2)
rm(list=ls())
#in_Ellipsoid <- stats::qchisq(0.95, 2)
setDTthreads(threads=1)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=2
#dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
dispersals<-c(0:2)
exposure_threshold<-1
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
#group<-"Birds"
#item<-data.table(sp="Cinnyris_jugularis", area=1)
print(paste(getDTthreads(), "CPUs are using in data.table"))
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  
  target_folders<-c(sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/dispersal", target_folder)
    
    
    j=9
    for (j in c(1:nrow(layer_df))){
      
      
      layer_item<-layer_df[j,]
      print(paste(i, nrow(df_list), item$sp, layer_item$LABEL))
      year<-predict_range[1]
      model<-readRDS(sprintf("%s/fit.rda", target_folder))
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      selected_cols<-c("x", "y", "mask_index")
      start_dis<-unique(start_dis[, ..selected_cols])
      start_dis$exposure<-0
      k=2
      
      for (k in c(1:length(dispersals))){
        dispersal<-dispersals[k]
        prev_dis<-start_dis
        dispersal_log<-NULL
        tda_target<-sprintf("%s/%s_%d.rda", target, layer_item$LABEL, dispersal)
        if (file.exists(tda_target)){
          #next()
          print(paste("Testing", tda_target))
          
          temp_df<-NULL
          tryCatch({
            temp_df<-readRDS(tda_target)
          }, warning = function(w) {
            #
          }, error = function(e) {
            print("Error")
          }, finally = {
            #
          })
          if (!is.null(temp_df)){
            next()
          }
        }
        
        if (!dir.exists(target)){
          dir.create(target, showWarnings = F)
        }
        saveRDS(NULL, tda_target)
        for (year in predict_range){
          
          #print(year)
          if (nrow(prev_dis)==0){
            next()
          }
          range_x<-range(prev_dis$x)
          range_x<-c(range_x[1]-150000*dispersal, range_x[2]+150000*dispersal)
          range_y<-range(prev_dis$y)
          range_y<-c(range_y[1]-150000*dispersal, range_y[2]+150000*dispersal)
          
          env_item<-data.table(env_layers[[paste(layer_item$LABEL, year, sep="_")]])
          env_item<-env_item[(x %between% range_x)&(y %between% range_y)]
          env_item$dist<-env_item[, min_dist(x, y, prev_dis), by = 1:nrow(env_item)]$V1/100000
          
          if (dispersal==0){
            env_item<-env_item[dist<1]
          }else{
            env_item<-env_item[dist<=dispersal]
          }
          if (F){
            plot(env_item$x, env_item$y)
            points(env_item$x, env_item$y, col="blue")
            points(start_dis$x, start_dis$y, col="red")
          }
          env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
          env_item[is.na(exposure)]$exposure<-0
          if (F){
            env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                        (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-0
            exposure<-env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                                   (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure
            if (length(exposure)>0){
              env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                           (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-exposure+1
            }
            env_item<-env_item[exposure<exposure_threshold]
          }
          env_item<-env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                                (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]
          prev_dis<-env_item
          
          if (nrow(prev_dis)>0){
            prev_dis$M<-dispersal
            prev_dis$YEAR<-year
            dispersal_log<-bind(dispersal_log, prev_dis)
            selected_cols<-c("x", "y", "mask_index", "exposure")
            prev_dis<-unique(prev_dis[, ..selected_cols])
          }
        }
        print("Writing result")
        
        saveRDS(dispersal_log, tda_target)
        print("Done! Writing result")
        if (F){
          yyear=2021
          for (yyear in c(2021:2100)){
            item<-dispersal_log[YEAR==yyear]
            p<-ggplot(item)+geom_point(aes(x=x, y=y, color=exposure))+
              xlim(range(dispersal_log$x))+ylim(range(dispersal_log$y))+ggtitle(yyear)
            print(p)
            x<-readline(prompt="X=exit: ")
            if (toupper(x)=="X"){
              break()
            }
          }
        }
      }
      
    }
  }
  
}
