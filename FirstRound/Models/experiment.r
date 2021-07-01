library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(batchtools)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
group<-"Amphibians"
source("commonFuns/addEllipse.R")
source("commonFuns/genCircle.R")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(1850:2020)
mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
#start_env_layers<-start_env_layers[,list(PR=mean(PR),TEMP=mean(TEMP)),by=.(x, y, mask_index, year)]

#item<-data.table(sp="Limnodynastes_interioris", area=1)
group<-"Amphibians"
item<-data.table(sp="Euphlyctis_hexadactylus", area=1)
group<-"Birds"
item<-data.table(sp="Cinnyris_jugularis", area=1)
occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
v<-raster::extract(mask, occ[, c("x", "y")])
all_v<-start_env_layers[mask_index %in% v]

p<-ggplot(start_env_layers[sample(nrow(start_env_layers), 5000),], aes(x=TEMP, y=PR))+geom_point(color="black", alpha=0.2, size=0.1)+theme_bw()

# To address this, we estimated species-realized niche limits using the climate projections from the
# historical run of each climate model (1850–2020).
# To prevent estimates of niche being inflated by either extreme outliers in the time series or from 
# the overestimation of species ranges, we excluded outlier temperature values within each
# grid cell, defined as those more than three standard deviations from the mean.
sampled_all_v<-all_v

mean_PR<-mean(sampled_all_v$PR, na.rm=T)
sd_PR<-sd(sampled_all_v$PR, na.rm=T)
offset_PR<-sd_PR*3
range_PR<-c(mean_PR-offset_PR, mean_PR+offset_PR)

mean_TEMP<-mean(sampled_all_v$TEMP, na.rm=T)
sd_TEMP<-sd(sampled_all_v$TEMP, na.rm=T)
offset_TEMP<-sd_TEMP*3
range_TEMP<-c(mean_TEMP-offset_TEMP, mean_TEMP+offset_TEMP)

fit_3sd<-data.frame(mean_PR=mean_PR, sd_PR=sd_PR, range_PR_min=range_PR[1], range_PR_max=range_PR[2],
                mean_TEMP=mean_TEMP, sd_TEMP=sd_TEMP, range_TEMP_min=range_TEMP[1], range_TEMP_max=range_TEMP[2])


colors<-c("blue", "red")
p<-p+geom_point(data=all_v, aes(x=TEMP, y=PR), size=0.1, color="purple", alpha=0.5)
box<-data.frame(x=c(fit_3sd$range_TEMP_min, fit_3sd$range_TEMP_min, fit_3sd$range_TEMP_max, fit_3sd$range_TEMP_max, fit_3sd$range_TEMP_min),
       y=c(fit_3sd$range_PR_max, fit_3sd$range_PR_min, fit_3sd$range_PR_min, fit_3sd$range_PR_max, fit_3sd$range_PR_max))
p<-p+geom_path(data=box, aes(x=x, y=y), color="red")

quantile_PR<-quantile(sampled_all_v$PR, c(0.25, 0.75), na.rm=T)
IQR_PR<-quantile_PR[2]-quantile_PR[1]
offset_PR<-IQR_PR*1.5
range_PR<-c(mean_PR-offset_PR, mean_PR+offset_PR)

quantile_TEMP<-quantile(sampled_all_v$TEMP, c(0.25, 0.75), na.rm=T)
IQR_TEMP<-quantile_TEMP[2]-quantile_TEMP[1]
offset_TEMP<-IQR_TEMP*1.5
range_TEMP<-c(mean_TEMP-offset_TEMP, mean_TEMP+offset_TEMP)

fit_IQR<-data.frame(mean_PR=mean_PR, sd_PR=sd_PR, range_PR_min=range_PR[1], range_PR_max=range_PR[2],
                mean_TEMP=mean_TEMP, sd_TEMP=sd_TEMP, range_TEMP_min=range_TEMP[1], range_TEMP_max=range_TEMP[2])
box<-data.frame(x=c(fit_IQR$range_TEMP_min, fit_IQR$range_TEMP_min, fit_IQR$range_TEMP_max, fit_IQR$range_TEMP_max, fit_IQR$range_TEMP_min),
                y=c(fit_IQR$range_PR_max, fit_IQR$range_PR_min, fit_IQR$range_PR_min, fit_IQR$range_PR_max, fit_IQR$range_PR_max))
p<-p+geom_path(data=box, aes(x=x, y=y), color="blue")

NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)


fit_MVE <- cov.rob(all_v[, c("PR", "TEMP")], quantile.used=NDquntil(nrow(all_v), 0.95),  method = "mve")



lines<-data.frame(addEllipse(fit_MVE$center, fit_MVE$cov, col="red", p.interval=0.95))
p<-p+geom_path(data=lines, aes(x=X2, y=X1), color="red")
ggsave(p, filename=sprintf("../../Figures/Experiment/%s_%s_E.png", group, item$sp))

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
predict_range<-c(2021:2100)

layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
print("Loading ENV DATA ...")
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")
dispersal<-1
layer_item<-layer_df[9,]
target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
selected_cols<-c("x", "y", "mask_index")
start_dis<-unique(start_dis[, ..selected_cols])
start_dis$exposure<-0

prev_dis<-start_dis
dispersal_log_3sd<-NULL
year=2021
exposure_threshold=1
for (year in predict_range){
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
  
  env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
  env_item[is.na(exposure)]$exposure<-0
  
  model<-fit_3sd
  env_item[((PR %between% c(model$range_PR_min, model$range_PR_max))&
              (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-0
  exposure<-env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                         (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure
  if (length(exposure)>0){
    env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                 (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-exposure+1
  }
  env_item<-env_item[exposure<exposure_threshold]
  prev_dis<-env_item
  selected_cols<-c("x", "y", "mask_index", "exposure")
  prev_dis<-unique(prev_dis[, ..selected_cols])
  if (nrow(prev_dis)>0){
    prev_dis$M<-dispersal
    prev_dis$YEAR<-year
    dispersal_log_3sd<-bind(dispersal_log_3sd, prev_dis)
  }
}
p<-ggplot()+
  geom_point(data=dispersal_log_3sd[YEAR==2100], aes(x=x, y=y, color=factor(exposure)))+
  geom_point(data=dispersal_log_3sd[YEAR==2021], aes(x=x, y=y), color="green")+
  geom_point(data=start_dis, aes(x=x, y=y), color="black")
ggsave(p, filename=sprintf("../../Figures/Experiment/%s_%s_G_3sd_1.png", group, item$sp))



prev_dis<-start_dis
dispersal_log_3sd<-NULL
year=2021
exposure_threshold=5
for (year in predict_range){
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
  
  env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
  env_item[is.na(exposure)]$exposure<-0
  
  model<-fit_3sd
  env_item[((PR %between% c(model$range_PR_min, model$range_PR_max))&
              (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-0
  exposure<-env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                         (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure
  if (length(exposure)>0){
    env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                 (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-exposure+1
  }
  env_item<-env_item[exposure<exposure_threshold]
  prev_dis<-env_item
  selected_cols<-c("x", "y", "mask_index", "exposure")
  prev_dis<-unique(prev_dis[, ..selected_cols])
  if (nrow(prev_dis)>0){
    prev_dis$M<-dispersal
    prev_dis$YEAR<-year
    dispersal_log_3sd<-bind(dispersal_log_3sd, prev_dis)
  }
}
p<-ggplot()+
  geom_point(data=dispersal_log_3sd[YEAR==2100], aes(x=x, y=y, color=exposure))+
  geom_point(data=dispersal_log_3sd[YEAR==2021], aes(x=x, y=y), color="green")+
  geom_point(data=start_dis, aes(x=x, y=y), color="black")
ggsave(p, filename="../../Figures/Experiment/G_3sd_5.png")



prev_dis<-start_dis

dispersal_log_IQR<-NULL

year=2021
exposure_threshold=1
for (year in predict_range){
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
  
  env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
  env_item[is.na(exposure)]$exposure<-0
  
  model<-fit_IQR
  env_item[((PR %between% c(model$range_PR_min, model$range_PR_max))&
              (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-0
  exposure<-env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                         (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure
  if (length(exposure)>0){
    env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                 (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-exposure+1
  }
  env_item<-env_item[exposure<exposure_threshold]
  prev_dis<-env_item
  selected_cols<-c("x", "y", "mask_index", "exposure")
  prev_dis<-unique(prev_dis[, ..selected_cols])
  if (nrow(prev_dis)>0){
    prev_dis$M<-dispersal
    prev_dis$YEAR<-year
    dispersal_log_IQR<-bind(dispersal_log_IQR, prev_dis)
  }
}
p<-ggplot()+
  geom_point(data=dispersal_log_IQR[YEAR==2100], aes(x=x, y=y, color=factor(exposure)))+
  geom_point(data=dispersal_log_IQR[YEAR==2021], aes(x=x, y=y), color="green")+
  geom_point(data=start_dis, aes(x=x, y=y), color="black")
ggsave(p, filename=sprintf("../../Figures/Experiment/%s_%s_G_IQR_1.png", group, item$sp))



prev_dis<-start_dis
dispersal_log_3sd<-NULL
dispersal_log_IQR<-NULL
dispersal_log_MVE<-NULL
year=2021
exposure_threshold=5
for (year in predict_range){
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
  
  env_item<-ljoin(env_item, prev_dis, by= c("x", "y", "mask_index"))
  env_item[is.na(exposure)]$exposure<-0
  
  model<-fit_IQR
  env_item[((PR %between% c(model$range_PR_min, model$range_PR_max))&
              (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-0
  exposure<-env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                         (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure
  if (length(exposure)>0){
    env_item[!((PR %between% c(model$range_PR_min, model$range_PR_max))&
                 (TEMP %between% c(model$range_TEMP_min, model$range_TEMP_max)))]$exposure<-exposure+1
  }
  env_item<-env_item[exposure<exposure_threshold]
  prev_dis<-env_item
  selected_cols<-c("x", "y", "mask_index", "exposure")
  prev_dis<-unique(prev_dis[, ..selected_cols])
  if (nrow(prev_dis)>0){
    prev_dis$M<-dispersal
    prev_dis$YEAR<-year
    dispersal_log_IQR<-bind(dispersal_log_IQR, prev_dis)
  }
}
p<-ggplot()+
  geom_point(data=dispersal_log_IQR[YEAR==2100], aes(x=x, y=y, color=exposure))+
  geom_point(data=dispersal_log_IQR[YEAR==2021], aes(x=x, y=y), color="green")+
  geom_point(data=start_dis, aes(x=x, y=y), color="black")
ggsave(p, filename="../../Figures/Experiment/G_IQR_5.png")


testsp<-"Euphlyctis_hexadactylus/Amphinians Cinnyris_jugularis/Birds"
dispersal_log<-readRDS("../../Objects/Niche_Models/Amphibians/Euphlyctis_hexadactylus/dispersal/EC-Earth3-Veg_SSP119_2.rda")
model<-readRDS("../../Objects/Niche_Models/Amphibians/Euphlyctis_hexadactylus/fit.rda")
yyear=2100
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
