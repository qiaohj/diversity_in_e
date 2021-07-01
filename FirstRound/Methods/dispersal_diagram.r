library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(ggpubr)
library(data.table)
library(batchtools)

setwd("Y:/Script/diversity_in_e")

source("commonFuns/colors.r")
source("commonFuns/functions.r")

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}

group<-"Amphibians"

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(1850:2020)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
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
  j=1
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area<=0){
      next()
    }
    
    
    print(paste(i, nrow(df_list), item$sp))
    occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
    if ((nrow(occ)>20)&(nrow(occ)<25)){
      if (j==1){
        asdf
      }
      j=j+1
    }else{
      next()
    }
  }
}
example_sp<-"Dendropsophus_walfordi"
target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
all_v<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
#plot(mask)
target<-sprintf("%s/dispersal_5", target_folder)
dispersal_df<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1_1"))

year=2021
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

layer_item<-layer_df[9,]

start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
selected_cols<-c("x", "y", "mask_index")
start_dis<-unique(start_dis[, ..selected_cols])
start_dis$exposure<-0
start_dis$is_new<-F
dispersal<-1
predict_range<-c(2021:2100)
env_layers<-readRDS("../../Objects/stacked_layers_2021_2100.rda")
model<-readRDS(sprintf("%s/fit.rda", target_folder))
exposure_threshold<-5
xrange<-c(-8000000, -3000000)
yrange<-c(-3300000, 1500000)
  
  
prev_dis<-start_dis
dispersal_log<-NULL
year<-2021
image_index<-1
for (year in predict_range){
  print(year)
  current_env<-data.table(env_layers[[paste(layer_item$LABEL, year, sep="_")]])
  
  #Suitable area
  suitable<-current_env[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                        (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]
  p<-ggplot()+geom_tile(data=suitable, 
                        aes(x=x, y=y), fill=colors_black[2])+
    xlim(xrange)+
    ylim(yrange)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none")+
    ggtitle(year-1)
  #print(year)
  if (nrow(prev_dis)==0){
    next()
  }
  disperable<-prev_dis[exposure==0]
  undisperable<-prev_dis[exposure!=0]
  if (nrow(disperable)>0){
    p1<-p+geom_tile(data=disperable, 
                aes(x=x, y=y), fill=colors_red[3])
  }
  if (nrow(undisperable)>0){
    p1<-p1+geom_tile(data=undisperable, 
                aes(x=x, y=y), fill=colors_blue[3])+
      geom_text(data=undisperable, aes(x=x, y=y, label=exposure), size=2)
  }
  ggsave(p1, filename=sprintf("../../Figures/Methods/dispersal_diagram/png/%d.png", image_index),
         width=6, height=6)
  image_index<-image_index+1
  
  ggsave(p1, filename=sprintf("../../Figures/Methods/dispersal_diagram/pdf/p1_%d.pdf", year),
         width=6, height=6)
  
  #p_item<-p_item+geom_point(data=prev_dis, aes(x=x, y=y))
  env_item_ori<-current_env
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
  
  
  new_pixels<-env_item[!(mask_index %in% c(prev_dis$mask_index))]
  if (nrow(new_pixels)>0){
    p2<-p1+geom_tile(data=new_pixels, aes(x=x, y=y), fill=colors_green[4])
  }else{
    p2<-p1
  }
  ggsave(p2, filename=sprintf("../../Figures/Methods/dispersal_diagram/png/%d.png", image_index),
         width=6, height=6)
  image_index<-image_index+1
  ggsave(p2, filename=sprintf("../../Figures/Methods/dispersal_diagram/pdf/p2_%d.pdf", year),
         width=6, height=6)
  #Step 1. Remove the new/unsuitable pixels
  env_item<-env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                        (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))|
                       (!is_new)]
  #Step 2. reset the exposure of suitable area to 0
  env_item[((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
              (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-0
  #p_item+geom_point(data=env_item, aes(x=x, y=y, color=factor(exposure)))
  
  #Step 3. increase exposure of unsuitable areas
  exposure<-env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                         (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure
  if (length(exposure)>0){
    env_item[!((PR %between% c(model$range_PR_sd_min, model$range_PR_sd_max))&
                 (TEMP %between% c(model$range_TEMP_sd_min, model$range_TEMP_sd_max)))]$exposure<-exposure+1
  }
  exposure_pixels<-env_item[exposure>0]
  p3<-p+geom_tile(data=env_item, aes(x=x, y=y, fill=factor(exposure)))+
    geom_text(data=env_item[exposure!=0], aes(x=x, y=y, label=exposure), size=2)+
    scale_fill_manual(breaks=c(0:5), values=c(colors_red[3], rep(colors_blue[3], 4),
                                              colors_blue[7]))
  #ggsave(p3, filename=sprintf("../../Figures/Methods/dispersal_diagram/png/p3_%d.png", year),
  #       width=6, height=6)
  ggsave(p3, filename=sprintf("../../Figures/Methods/dispersal_diagram/pdf/p3_%d.pdf", year),
         width=6, height=6)
  #Step 4. Remove the exposure >=5
  is_p4<-(max(env_item$exposure)==5)
  env_item<-env_item[exposure<exposure_threshold]
  if (is_p4){
    p4<-p+geom_tile(data=env_item, aes(x=x, y=y, fill=factor(exposure)))+
      geom_text(data=env_item[exposure!=0], aes(x=x, y=y, label=exposure), size=2)+
      scale_fill_manual(breaks=c(0:5), values=c(colors_red[3], rep(colors_blue[3], 4),
                                                colors_blue[7]))
    
    #ggsave(p4, filename=sprintf("../../Figures/Methods/dispersal_diagram/png/p4_%d.png", year),
    #       width=6, height=6)
    ggsave(p4, filename=sprintf("../../Figures/Methods/dispersal_diagram/pdf/p4_%d.pdf", year),
           width=6, height=6)
  }
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

#ffmpeg -r 1 -start_number 2015 -i %04d.png -y ../dispersal_diagram.gif
#ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../dispersal_diagram.mp4
