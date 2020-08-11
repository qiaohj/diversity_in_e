library(raster)
#library(rgdal)
#library(rgeos)
#library(MASS)
#library(cluster)
library(dplyr)
#library(ggplot2)
source("addEllipse.R")
source("genCircle.R")
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}
bind<-function(df1, df2){
  if (is.null(df1)){
    df1<-df2
  }else{
    df1<-dplyr::bind_rows(df1, df2)
  }
  return(df1)
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

predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
for (i in c(1:nrow(df_list))){
  
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  #target_folders<-c(sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp),
  #                  sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folders<-c(sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, item$sp))
  target_folder<-target_folders[2]
  for (target_folder in target_folders){
    target<-sprintf("%s/dispersal", target_folder)
    model<-"Normal"
    if (grepl("Mean_GCM", target)){
      model<-"Mean"
    }
    if (dir.exists(target)){
      next()
    }
    
    dir.create(target, showWarnings = F)
    j=1
    for (j in c(1:nrow(layer_df))){
      layer_item<-layer_df[j,]
      year<-predict_range[1]
      enm_folder<-sprintf("%s/predict", target_folder)
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(X, Y)
      colnames(start_dis)<-c("x", "y")
      k=1
      
      print("Init all potential distributions")
      distributoins<-list()
      for (year in predict_range){
        #print(year)
        env_item<-readRDS(sprintf("%s/%s_%d.rda", enm_folder, layer_item$LABEL, year))
        distributoins[[as.character(year)]]<-env_item
      }
      for (k in c(1:nrow(dispersals))){
        dispersal<-dispersals[k,]
        prev_dis<-start_dis
        dispersal_log<-NULL
        for (year in predict_range){
          print(paste(item$sp, layer_item$LABEL, year, 
                      j, ":", nrow(layer_df), "/",
                      i, ":", nrow(df_list), "/",
                      k, ":", nrow(dispersals), "/",
                      model))
          env_item<-distributoins[[as.character(year)]]
          for (mt in c(1:dispersal$N)){
            if (nrow(prev_dis)==0){
              next()
            }
            env_item<-env_item%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, prev_dis)/100000)
            if (dispersal$M==0){
              prev_dis<-env_item%>%dplyr::filter(dist<1)
            }else{
              prev_dis<-env_item%>%dplyr::filter(dist<=dispersal$M)
            }
            if (nrow(prev_dis)>0){
              prev_dis$M<-dispersal$M
              prev_dis$N<-dispersal$N
              prev_dis$N_Step<-mt
              prev_dis$YEAR<-year
              dispersal_log<-bind(dispersal_log, prev_dis)
            }
            if (F){
              ggplot(env_item, aes(x=x, y=y, color=dist)) + geom_point() +
                geom_point(data=env_item%>%dplyr::filter(dist<=dispersal$M), color="purple")+
                geom_point(data=prev_dis, aes(x=x, y=y), color="red")
            }
          }
        }
        saveRDS(dispersal_log, sprintf("%s/%s_%d_%d.rda", target, layer_item$LABEL, dispersal$M, dispersal$N))
        if (F){
          for (year in predict_range){
            env_item<-distributoins[[as.character(year)]]
            p<-ggplot(env_item, aes(x=x, y=y), color="grey") + geom_point() +
              geom_point(data=dispersal_log%>%dplyr::filter(YEAR==year), aes(x=x, y=y), color="red")+
              theme_bw()+ggtitle(year)+xlim(5096782, 15096782)+ylim(-1525228, 2525228)
            #print(p)
            ggsave(p, file=sprintf("../../Figures/Dispersal_Example/1_1/%d.png", year))
            
          }
        }
      }
      
    }
  }
  
}
