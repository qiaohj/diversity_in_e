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
  group<-"Mammals"
}
mask<-raster("../../Raster/mask_index.tif")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")

#mask<-raster("../../Raster/mask_index.tif")
#mask_p<-data.frame(rasterToPoints(mask))
predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=2
#dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
dispersals<-c(0:1)
exposure_threshold<-5
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
  
  target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
  
  target<-sprintf("%s/dispersal_5", target_folder)
  
  if (dir.exists(target)){
    next()
  }
  
  dir.create(target, showWarnings = F)
  j=1
  for (j in c(1:nrow(layer_df))){
    
    
    layer_item<-layer_df[j,]
    print(paste(i, nrow(df_list), item$sp, layer_item$LABEL))
    year<-predict_range[1]
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
    start_dis$mask_index<-extract(mask, data.frame(start_dis[, c("x", "y")]))
    saveRDS(start_dis, sprintf("%s/occ_with_env.rda", target_folder))
    selected_cols<-c("x", "y", "mask_index")
    start_dis<-unique(start_dis[, ..selected_cols])
    start_dis$exposure<-0
    start_dis$is_new<-F
    k=2
    
    for (k in c(1:length(dispersals))){
      dispersal<-dispersals[k]
      prev_dis<-start_dis
      dispersal_log<-NULL
      for (year in predict_range){
        
        #print(year)
        if (nrow(prev_dis)==0){
          next()
        }
        disperable<-prev_dis[exposure==0]
        undisperable<-prev_dis[exposure!=0]
        #p_item<-p_item+geom_point(data=prev_dis, aes(x=x, y=y))
        env_item_ori<-data.table(env_layers[[paste(layer_item$LABEL, year, sep="_")]])
        env_item<-env_item_ori
        if (nrow(disperable)>0){
          range_x<-range(disperable$x)
          range_x<-c(range_x[1]-150000*dispersal, range_x[2]+150000*dispersal)
          range_y<-range(disperable$y)
          range_y<-c(range_y[1]-150000*dispersal, range_y[2]+150000*dispersal)
          
          
          env_item<-env_item[(x %between% range_x)&(y %between% range_y)]
          env_item$dist<-env_item[, min_dist(x, y, disperable), by = 1:nrow(env_item)]$V1/100000
          
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
          undisperable<-undisperable[!(mask_index %in% env_item$mask_index)]
          undisperable<-ljoin(undisperable, env_item_ori, by= c("x", "y", "mask_index"))
          undisperable$dist<-0
          env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
          env_item[is.na(exposure)]$exposure<-0
          env_item[is.na(is_new)]$is_new<-T
          env_item<-bind(env_item, undisperable)
          #p_item<-p_item+geom_point(data=env_item, aes(x=x, y=y, color=factor(is_new)))
        }else{
          env_item<-ijoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
          env_item$dist<-0
        }
        
        #Step 1. Remove the new/unsuitable pixels
        env_item<-env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                              (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                              (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))|
                             (!is_new)]
        #Step 2. reset the exposure of suitable area to 0
        env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                    (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                    (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-0
        #p_item+geom_point(data=env_item, aes(x=x, y=y, color=factor(exposure)))
        
        #Step 3. increase exposure of unsuitable areas
        exposure<-env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                               (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                               (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure
        if (length(exposure)>0){
          env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                       (TEMP_MAX %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max))&
                       (TEMP_MIN %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-exposure+1
        }
        
        #Step 4. Remove the exposure >=5
        env_item<-env_item[exposure<exposure_threshold]
        #Step 5. set is_new
        env_item$is_new<-F
        
        
        prev_dis<-env_item
        
        if (nrow(prev_dis)>0){
          prev_dis$M<-dispersal
          prev_dis$YEAR<-year
          dispersal_log<-bind(dispersal_log, prev_dis)
          selected_cols<-c("x", "y", "mask_index", "exposure", "is_new")
          prev_dis<-unique(prev_dis[, ..selected_cols])
          if (F){
            p_item<-ggplot(prev_dis, aes(x=x, y=y, fill=factor(exposure)))+geom_tile()+
              xlim(c(996781.7, 3696781.7))+ylim(c(-2174772.1, 925227.9))+
              ggtitle(year)
            print(p_item)
            x<-readline(prompt="X=exit: ")
            if (toupper(x)=="X"){
              break()
            }
          }
          #prev_dis<-prev_dis[exposure==0]
        }
      }
      print("Writing result")
      saveRDS(dispersal_log, sprintf("%s/%s_%d.rda", target, layer_item$LABEL, dispersal))
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
