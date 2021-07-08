
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

if (F){
  source("commonFuns/functions.r")
  
  keyspots<-read.csv("../../Objects/keyspots/keyspots.csv", stringsAsFactors = F)
  mask<-raster("../../Raster/mask_100km.tif")
  plot(mask)
  cord.dec = SpatialPoints(keyspots[, c("lon", "lat")], proj4string=CRS("+proj=longlat"))
  cord.eck4 <- spTransform(cord.dec, crs(mask))
  plot(cord.eck4, add=T)
  keyspots$x<-cord.eck4@coords[,1]
  keyspots$y<-cord.eck4@coords[,2]
  
  pc100km <- gBuffer( cord.eck4, width=100*10e3, byid=TRUE )
  plot(pc100km, add=T)
  pc100km <- SpatialPolygonsDataFrame(pc100km, data=keyspots)
  writeOGR( pc100km, "../../Shape/pc100km", "pc100km", driver="ESRI Shapefile" ) 
  #name<-pc100km$name[1]
  #plot(pc100km[which(pc100km$name==name),], col="red", add=T)
  
  
  #points(over_points_result$x, over_points_result$y, add=T, col="red")
  #head(over_points_result)
  
  all_df<-NULL
  for (group in c("Amphibians", "Reptiles", "Mammals", "Birds")){
    print(group)
    df<-readRDS(sprintf("../../Objects/Species_exposure/%s_exposure_se_acc.rda", group))
    over_points_result_se<-df%>%dplyr::group_by(SSP, year)%>%
      dplyr::summarise(mean=mean(exposure_ratio_mean),
                       sd=sd(exposure_ratio_mean),
                       CI=CI(exposure_ratio_mean)[2]-CI(exposure_ratio_mean)[3],
                       mean_type_0_ratio=mean(range_type_0_ratio_mean),
                       sd_type_0_ratio=sd(range_type_0_ratio_mean),
                       CI_type_0_ratio=CI(range_type_0_ratio_mean)[2]-CI(range_type_0_ratio_mean)[3],
                       mean_type_1_ratio=mean(range_type_1_ratio_mean),
                       sd_type_1_ratio=sd(range_type_1_ratio_mean),
                       CI_type_1_ratio=CI(range_type_1_ratio_mean)[2]-CI(range_type_1_ratio_mean)[3],
                       mean_type_2_ratio=mean(range_type_2_ratio_mean),
                       sd_type_2_ratio=sd(range_type_2_ratio_mean),
                       CI_type_2_ratio=CI(range_type_2_ratio_mean)[2]-CI(range_type_2_ratio_mean)[3],
                       mean_type_3_ratio=mean(range_type_3_ratio_mean),
                       sd_type_3_ratio=sd(range_type_3_ratio_mean),
                       CI_type_3_ratio=CI(range_type_3_ratio_mean)[2]-CI(range_type_3_ratio_mean)[3],
                       mean_type_4_ratio=mean(range_type_4_ratio_mean),
                       sd_type_4_ratio=sd(range_type_4_ratio_mean),
                       CI_type_4_ratio=CI(range_type_4_ratio_mean)[2]-CI(range_type_4_ratio_mean)[3])
    over_points_result_se$group<-group
    over_points_result_se$name<-"Global"
    over_points_result_se$lon<-0
    over_points_result_se$lat<-0
    all_df<-bind(all_df, over_points_result_se)
    sp_exposure<-SpatialPoints(df[, c("x", "y")], proj4string=crs(mask))
    
    for (i in c(1:nrow(keyspots))){
      keyspot<-keyspots[i,]
      over_points<-over(sp_exposure, pc100km[which(pc100km$name==keyspot$name),])
      over_points_result<-df[!is.na(over_points$name),]
      over_points_result_se<-over_points_result%>%dplyr::group_by(SSP, year)%>%
        dplyr::summarise(mean=mean(exposure_ratio_mean),
                         sd=sd(exposure_ratio_mean),
                         CI=CI(exposure_ratio_mean)[2]-CI(exposure_ratio_mean)[3],
                         mean_type_0_ratio=mean(range_type_0_ratio_mean),
                         sd_type_0_ratio=sd(range_type_0_ratio_mean),
                         CI_type_0_ratio=CI(range_type_0_ratio_mean)[2]-CI(range_type_0_ratio_mean)[3],
                         mean_type_1_ratio=mean(range_type_1_ratio_mean),
                         sd_type_1_ratio=sd(range_type_1_ratio_mean),
                         CI_type_1_ratio=CI(range_type_1_ratio_mean)[2]-CI(range_type_1_ratio_mean)[3],
                         mean_type_2_ratio=mean(range_type_2_ratio_mean),
                         sd_type_2_ratio=sd(range_type_2_ratio_mean),
                         CI_type_2_ratio=CI(range_type_2_ratio_mean)[2]-CI(range_type_2_ratio_mean)[3],
                         mean_type_3_ratio=mean(range_type_3_ratio_mean),
                         sd_type_3_ratio=sd(range_type_3_ratio_mean),
                         CI_type_3_ratio=CI(range_type_3_ratio_mean)[2]-CI(range_type_3_ratio_mean)[3],
                         mean_type_4_ratio=mean(range_type_4_ratio_mean),
                         sd_type_4_ratio=sd(range_type_4_ratio_mean),
                         CI_type_4_ratio=CI(range_type_4_ratio_mean)[2]-CI(range_type_4_ratio_mean)[3])
      over_points_result_se$group<-group
      over_points_result_se$name<-keyspot$name
      over_points_result_se$lon<-keyspot$lon
      over_points_result_se$lat<-keyspot$lat
      all_df<-bind(all_df, over_points_result_se)
    }
  }
  saveRDS(all_df, "../../Figures/Species_Exposure/Species_Exposure.rda")
  
  envs<-readRDS("../../Objects/stacked_layers_2015_2100.rda")
  label<- names(envs)[1]
  env_list<-list()
  for (label in names(envs)){
    labeles<-strsplit(label, "_")[[1]]
    GCM<-labeles[1]
    SSP<-labeles[2]
    year<-as.numeric(labeles[3])
    label_name<-paste(SSP, year, sep="_")
    envs[[label]]$GCM<-GCM
    envs[[label]]$SSP<-SSP
    envs[[label]]$year<-year
    if (label_name %in% names(env_list)){
      env_list[[label_name]]<-bind_rows(env_list[[label_name]], envs[[label]])
    }else{
      env_list[[label_name]]<-envs[[label]]
    }
    
  }
  label<-names(env_list)[1]
  all_env_df<-NULL
  for (label in names(env_list)){
    labeles<-strsplit(label, "_")[[1]]
    SSP<-labeles[1]
    year<-as.numeric(labeles[2])
    
    env_df<-env_list[[label]]%>%dplyr::group_by(SSP)%>%
      dplyr::summarise(mean_PR=mean(PR),
                       mean_TEMP_MAX=mean(TEMP_MAX),
                       mean_TEMP_MIN=mean(TEMP_MIN))
    env_df$SSP<-SSP
    env_df$year<-year
    env_df$name<-"Global"
    all_env_df<-bind(all_env_df, env_df)
    points<-SpatialPoints(env_list[[label]][, c("x", "y")], proj4string=crs(mask))
    for (i in c(1:nrow(keyspots))){
      keyspot<-keyspots[i,]
      over_points<-over(points, pc100km[which(pc100km$name==keyspot$name),])
      over_points_result<-env_list[[label]][!is.na(over_points$name),]
      env_df<-over_points_result%>%dplyr::group_by(SSP, year)%>%
        dplyr::summarise(mean_PR=mean(PR),
                         mean_TEMP_MAX=mean(TEMP_MAX),
                         mean_TEMP_MIN=mean(TEMP_MIN))
      env_df$SSP<-SSP
      env_df$year<-year
      env_df$name<-keyspot$name
      all_env_df<-bind(all_env_df, env_df)
    }
  }
  saveRDS(all_env_df, "../../Figures/Species_Exposure/all_env_df.rda")
  
  
}

#1: in range, but out mve
#2: out of range because of the max temp
#3: out of range because of the prec
#4: out of range because of the prec and temp

source("colors.r")
all_df<-readRDS("../../Figures/Species_Exposure/Species_Exposure.rda")


p<-ggplot(all_df, aes(x=year, y=mean, color=factor(group)))+
  geom_line(aes(linetype=factor(SSP)))+
  theme_bw()+
  facet_wrap(vars(name))
ggsave(p, file="../../Figures/Species_Exposure/species_exposure_overall.png", width=12, height=6)
names(all_df)[21]<-"Group"
p<-ggplot(all_df%>%dplyr::filter(name=="Global"), aes(x=year, y=mean, color=Group))+
  geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI, color=Group), alpha=0.7, width=0.25)+
  geom_line(aes(linetype=SSP))+
  scale_color_manual(values=color_groups)+
  scale_linetype_manual(values=linetype_ssp)+
  xlab("Year")+
  ylab("Species exposure proportion")+
  theme_bw()
p
ggsave(p, file="../../Figures/Species_Exposure/species_exposure_global_only.png", width=12, height=6)


all_env_df<-readRDS("../../Figures/Species_Exposure/all_env_df.rda")
p<-ggplot(all_env_df, aes(x=year, y=mean_TEMP_MAX))+
  geom_line(aes(linetype=factor(SSP)))+
  theme_bw()+
  facet_wrap(vars(name))
ggsave(p, file="../../Figures/Species_Exposure/species_exposure_temp.png", width=12, height=6)
p<-ggplot(all_env_df, aes(x=year, y=mean_PR))+
  geom_line(aes(linetype=factor(SSP)))+
  theme_bw()+
  facet_wrap(vars(name))
ggsave(p, file="../../Figures/Species_Exposure/species_exposure_prec.png", width=12, height=6)

df_type_all<-NULL
for (type in c(0:4)){
  df_item<-all_df[, c("SSP", "year", "name", "group", sprintf("mean_type_%d_ratio", type))]
  colnames(df_item)<-c("SSP", "year", "name", "group", "ratio")
  if (type==0){
    df_item$type<-"IN NICHE"
  }
  if (type==1){
    df_item$type<-"IN RANGE BUT OUT NICHE"
  }
  if (type==2){
    df_item$type<-"OUT OF RANGE BECAUSE OF TEMPERATURE"
  }
  if (type==3){
    df_item$type<-"OUT OF RANGE BECAUSE OF PRECIPITATION"
  }
  if (type==4){
    df_item$type<-"OUT OF RANGE BECAUSE OF BOTH"
  }
  df_type_all<-bind(df_type_all, df_item)
}

p<-ggplot(df_type_all[which(df_type_all$type!="IN NICHE"),], aes(x=year, y=ratio, color=factor(type)))+
  geom_line(aes(linetype=factor(SSP)))+
  theme_bw()+
  facet_wrap(~name+group, ncol=4, scale="free")
ggsave(p, file="../../Figures/Species_Exposure/species_exposure_env.png", width=12, height=18)

g<-"Amphibians"
year<-2100
SSP<-"SSP119"
pc100km<-readOGR("../../Shape/pc100km", "pc100km") 
png("../../Figures/Species_Exposure/Species_Exposure_map.png", width=1000, height=1000)
par(mfrow = c(4, 3))
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  
  for (SSP in c("SSP119", "SSP245", "SSP585")){
    r<-raster(sprintf("../../Objects/Species_exposure/Exposure_by_year/%s/TIF/%s_%d.tif", g, SSP, year))
    plot(r, main=paste(g, SSP, year), cex=2)
    plot(pc100km, add=T)
  }
  
}

dev.off()  
