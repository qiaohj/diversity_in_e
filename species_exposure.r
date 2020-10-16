library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
source("addEllipse.R")
source("genCircle.R")
source("mve_box.r")
source("functions.r")

NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")
var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"
predict_range<-c(2015:2100)
predict_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=predict_range)
predict_layer_df$VAR2<-"sum"
predict_layer_df[which(predict_layer_df$VAR=="tasmax"), "VAR2"]<-"max"
predict_layer_df[which(predict_layer_df$VAR=="tasmin"), "VAR2"]<-"min"
predict_layer_df$names<-sprintf("%s_%s_%d", predict_layer_df$GCM, predict_layer_df$VAR,
                              predict_layer_df$Y)
var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"
predict_layer_files<-sprintf(var_tamplate, predict_layer_df$GCM, predict_layer_df$SSP, predict_layer_df$VAR,
                           predict_layer_df$Y, predict_layer_df$VAR2)
predict_layers<-stack(predict_layer_files)
names(predict_layers)<-predict_layer_df$names
if (F){
  mask<-raster("../../Raster/mask.tif")
  v_no_na<-values(mask)[which(!is.na(values(mask)))]
  v_no_na<-c(1:length(v_no_na))
  values(mask)[which(!is.na(values(mask)))]<-v_no_na
  plot(mask)
  writeRaster(mask, "../../Raster/mask_index.tif")
}
if (F){
  mask<-data.frame(rasterToPoints(raster("../../Raster/mask_index.tif")))
  env_layers<-NULL
  for (i in c(1:nrow(layer_df))){
    print(paste("Init layer list:", i, nrow(layer_df)))
    item<-layer_df[i,]
    pr<-sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum")
    tasmax<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max")
    tasmin<-sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min")
    layers<-stack(c(pr, tasmax, tasmin))
    names(layers)<-c("PR", "TEMP_MAX", "TEMP_MIN")
    layers_table<-mask
    layers_table[,c("PR", "TEMP_MAX", "TEMP_MIN")]<-raster::extract(layers, mask[, c("x", "y")])
    layers_table$year<-item$Y
    layers_table$GCM<-item$GCM
    layers_table$SSP<-item$SSP
    
    env_layers<-bind(env_layers, layers_table)
  }
  saveRDS(env_layers, "../../Objects/stacked_layers_2015_2100_df.rda")
}
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))

layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=predict_range)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, layer_df$Y, sep="_")

print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2015_2100_df.rda")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=100
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folders<-c(sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folder<-target_folders[1]
  for (target_folder in target_folders){
    target<-sprintf("%s/exposure", target_folder)
    if (dir.exists(target)){
      next()
    }
    dir.create(target, showWarnings = F)
    print(paste(group, i, nrow(df_list), item$sp, target))
    
    model<-readRDS(sprintf("%s/fit.rda", target_folder))
    xy<-mve_box(model$center, model$cov, col="black", p.interval=0.95)
    x_range<-range(xy[,1])
    y_range<-range(xy[,2])
    
    occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
    v<-raster::extract(mask, occ[, c("x", "y")])
    #occ_predict<-mask_p%>%dplyr::filter(mask_index %in% v)
    
    env_item<-env_layers%>%dplyr::filter(mask_index %in% v)
    env_item$dist2 <- stats::mahalanobis(env_item[, c("PR", "TEMP_MAX")], center = model$center, 
                                         cov = model$cov)
    env_item$in_out_2<-F
    env_item[which(env_item$dist2<in_Ellipsoid), "in_out_2"]<-T
    env_item$in_mve<-env_item$in_out_2
    
    env_item$temp_max_range<-between(env_item$TEMP_MAX, y_range[1], y_range[2])
    env_item$prec_range<-between(env_item$PR, x_range[1], x_range[2])
    env_item$in_range<-env_item$temp_max_range&env_item$prec_range
    env_item$range_type<-0 # in mve
    env_item[which((env_item$in_range)&(!env_item$in_mve)), "range_type"]<-1 #in range, but out mve
    env_item[which((!env_item$in_range)&(!env_item$in_mve)&(!env_item$temp_max_range)), "range_type"]<-2 #out of range because of the max temp
    env_item[which((!env_item$in_range)&(!env_item$in_mve)&(!env_item$prec_range)), "range_type"]<-3 #out of range because of the prec
    env_item[which((!env_item$in_range)&(!env_item$in_mve)&(!env_item$temp_max_range)&(!env_item$prec_range)), "range_type"]<-4 #out of range because of the prec and temp
    #table(env_item$range_typ)
    #table(env_item$in_mve)
    
    saveRDS(env_item, sprintf("%s/exposure.rda", target))
    
    if (F){
      df_ori<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      colors<-rainbow(5)
      plot(env_item$PR, env_item$TEMP_MAX, col=colors[env_item$range_type+1], 
           xlim=range(env_item$PR), ylim=range(env_item$TEMP_MAX), pch=".")
      
      
      addEllipse(model$center, model$cov, col="black", p.interval=0.95)
      points(df_ori$PR, df_ori$TEMP, col="purple")
      
      box<-matrix(c(x_range[1], y_range[1], x_range[1], y_range[2], x_range[2], y_range[2], x_range[2], y_range[1], x_range[1], y_range[1]), ncol=2, byrow=T)
      lines(box)
      
    }
    
  }
  
}