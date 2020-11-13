library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("addEllipse.R")
source("genCircle.R")
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)
group<-"Amphibians"
sp<-"Abavorana_luctuosa"
occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, sp))

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

v<-extract(start_layers, occ[, c("x", "y")])
v<-data.frame(v)
col<-colnames(v)

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

j=1
all_v<-NULL
for (j in c(1:nrow(var_pair))){
  var_item<-var_pair[j,]
  v_item<-data.frame(PR=v[, var_item$PR_NAME], TEMP=v[, var_item$TEMP_NAME], 
                     X=occ$x, Y=occ$y, Year=var_item$Y, VAR=var_item$TEMP)
  if (is.null(all_v)){
    all_v<-v_item
  }else{
    all_v<-bind_rows(all_v, v_item)
  }
}
all_v_mean<-all_v%>%dplyr::group_by(Year, X, Y, VAR)%>%dplyr::summarise(PR=mean(PR),
                                                                   TEMP=mean(TEMP))

fit <- cov.rob(all_v[, c("PR", "TEMP")], quantile.used=NDquntil(nrow(all_v), 0.95),  method = "mve")
all_v$dist <- stats::mahalanobis(all_v[, c("PR", "TEMP")], center = fit$center, 
                                 cov = fit$cov)
all_v$in_out<-0
all_v[which(all_v$dist<in_Ellipsoid), "in_out"]<-1

fit_mean <- cov.rob(all_v_mean[, c("PR", "TEMP")], quantile.used=NDquntil(nrow(all_v_mean), 0.95),  method = "mve")
all_v_mean$dist <- stats::mahalanobis(all_v_mean[, c("PR", "TEMP")], center = fit_mean$center, 
                                 cov = fit_mean$cov)
all_v_mean$in_out<-0
all_v_mean[which(all_v_mean$dist<in_Ellipsoid), "in_out"]<-1
if (F){
  p_start_layers<-data.frame(rasterToPoints(start_layers))
  colors<-c("red", "blue")
  plot(p_start_layers$UKESM1_pr_2014, p_start_layers$UKESM1_tasmax_2014, pch=".", col="grey")
  points(p_start_layers$UKESM1_pr_2014, p_start_layers$UKESM1_tasmin_2014, pch=".", col="grey")
  points(all_v[which(all_v$in_out==1),]$PR, all_v[which(all_v$in_out==1),]$TEMP, col=colors[1])
  addEllipse(fit$center, fit$cov, col=colors[1], p.interval=0.95)
  points(all_v_mean[which(all_v_mean$in_out==1),]$PR, all_v_mean[which(all_v_mean$in_out==1),]$TEMP, col=colors[2])
  addEllipse(fit_mean$center, fit_mean$cov, col=colors[2], p.interval=0.95)
}