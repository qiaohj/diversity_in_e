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
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(2000:2014)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

if (F){
  mask<-data.frame(rasterToPoints(raster("../../Raster/mask_index.tif")))
  start_env_layers<-NULL
  for (i in c(1:nrow(start_layer_df))){
    print(paste("Init layer list:", i, nrow(start_layer_df)))
    item<-start_layer_df[i,]
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
    
    start_env_layers<-bind(start_env_layers, layers_table)
  }
  saveRDS(start_env_layers, "../../Objects/stacked_layers_2000_2014_df.rda")
}
mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_2000_2014_df.rda")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))


var_pair<-start_layer_df[which(start_layer_df$VAR=="pr"),]
var_pair1<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmax"),], by=c("GCM", "SSP", "Y"))
#var_pair2<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmin"),], by=c("GCM", "SSP", "Y"))
var_pair1<-var_pair1%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair1)<-c("GCM", "PR", "Y", "TEMP")
#var_pair2<-var_pair2%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
#colnames(var_pair2)<-c("GCM", "PR", "Y", "TEMP")
#var_pair<-bind_rows(var_pair1, var_pair2)
var_pair<-var_pair1
var_pair$PR_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$PR, var_pair$Y, sep="_")
var_pair$TEMP_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$TEMP, var_pair$Y, sep="_")
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
i=1
if (F){
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area>=100){
      asdf
    }
  }
}
for (i in c(1:nrow(df_list))){
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp)
  if (dir.exists(target_folder)){
    next()
  }
  dir.create(target_folder, showWarnings = F, recursive=T)
  print(paste(i, nrow(df_list), item$sp))
  occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
  v<-raster::extract(mask, occ[, c("x", "y")])
  v<-start_env_layers %>%dplyr::filter(mask_index %in% v)
  j=1
  
  all_v<-v%>%dplyr::group_by(year, x, y)%>%dplyr::summarise(PR=mean(PR),
                                                                TEMP=mean(TEMP_MAX))
  fit <- cov.rob(all_v[, c("PR", "TEMP")], quantile.used=NDquntil(nrow(all_v), 0.95),  method = "mve")
  saveRDS(fit, sprintf("%s/fit.rda", target_folder))
  
  all_v$dist <- stats::mahalanobis(all_v[, c("PR", "TEMP")], center = fit$center, 
                                   cov = fit$cov)
  all_v$in_out<-0
  all_v[which(all_v$dist<in_Ellipsoid), "in_out"]<-1
  if (F){
    table(all_v$in_out)
    colors<-c("red", "blue")
    plot(all_v$PR, all_v$TEMP, col=colors[all_v$in_out+1])
    addEllipse(fit$center, fit$cov, col="red", p.interval=0.95)
  }
  saveRDS(all_v, sprintf("%s/occ_with_env.rda", target_folder))
}

if (F){
  readRDS(sprintf("%s/occ_with_env.rda", target_folder))
  readRDS(sprintf("%s/%s.rda", target, layer_item$LABEL))
}