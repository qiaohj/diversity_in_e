library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
source("addEllipse.R")
source("genCircle.R")
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
start_range<-c(2000:2014)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)
start_layer_df$VAR2<-"sum"
start_layer_df[which(start_layer_df$VAR=="tasmax"), "VAR2"]<-"max"
start_layer_df[which(start_layer_df$VAR=="tasmin"), "VAR2"]<-"min"
start_layer_df$names<-sprintf("%s_%s_%d", start_layer_df$GCM, start_layer_df$VAR,
                              start_layer_df$Y)
var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"
start_layer_files<-sprintf(var_tamplate, start_layer_df$GCM, start_layer_df$SSP, start_layer_df$VAR,
                           start_layer_df$Y, start_layer_df$VAR2)
start_layers<-stack(start_layer_files)
names(start_layers)<-start_layer_df$names
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))

i=1
var_pair<-start_layer_df[which(start_layer_df$VAR=="pr"),]
var_pair1<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmax"),], by=c("GCM", "SSP", "Y"))
var_pair2<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmin"),], by=c("GCM", "SSP", "Y"))
var_pair1<-var_pair1%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair1)<-c("GCM", "PR", "Y", "TEMP")
var_pair2<-var_pair2%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair2)<-c("GCM", "PR", "Y", "TEMP")
var_pair<-bind_rows(var_pair1, var_pair2)
var_pair$PR_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$PR, var_pair$Y, sep="_")
var_pair$TEMP_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$TEMP, var_pair$Y, sep="_")
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
  if (dir.exists(target_folder)){
    next()
  }
  dir.create(target_folder, showWarnings = F)
  print(paste(i, nrow(df_list), item$sp))
  occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
  v<-extract(start_layers, occ[, c("x", "y")])
  v<-data.frame(v)
  col<-colnames(v)
  j=1
  all_v<-NULL
  for (j in c(1:nrow(var_pair))){
    var_item<-var_pair[j,]
    v_item<-data.frame(PR=v[, var_item$PR_NAME], TEMP=v[, var_item$TEMP_NAME], X=occ$x, Y=occ$y)
    if (is.null(all_v)){
      all_v<-v_item
    }else{
      all_v<-bind_rows(all_v, v_item)
    }
  }
  
  fit <- cov.rob(all_v[, c("PR", "TEMP")], quantile.used=NDquntil(nrow(all_v), 0.95),  method = "mve")
  saveRDS(fit, sprintf("%s/fit.rda", target_folder))
  
  all_v$dist <- stats::mahalanobis(all_v[, c("PR", "TEMP")], center = fit$center, 
                                 cov = fit$cov)
  all_v$in_out<-0
  all_v[which(all_v$dist<in_Ellipsoid), "in_out"]<-1
  if (F){
    colors<-c("red", "blue")
    plot(all_v$PR, all_v$TEMP, col=colors[all_v$in_out+1])
    addEllipse(fit$center, fit$cov, col="red", p.interval=0.95)
  }
  saveRDS(all_v, sprintf("%s/occ_with_env.rda", target_folder))
}
