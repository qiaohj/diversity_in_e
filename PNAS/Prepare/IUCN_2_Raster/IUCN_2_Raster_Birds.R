library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
library(rmapshaper)
 
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(3)
print(sprintf("Current core number is %d", getDTthreads()))
#for BIRDS
if (F){
  r<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/EC-Earth3-Veg/SSP119/1850/bio12_eck4_1km.tif")
  plot(r)
  continent<-readOGR("../../Shape/continents", "continent")
  crs(continent)<-CRS("+proj=longlat +datum=WGS84")
  continent_eck4<-spTransform(continent, crs(r))
  continent_eck4<-continent_eck4[which(continent_eck4$CONTINENT!="Antarctica"),]
  writeOGR(continent_eck4, "../../Shape/continents", "continent_eck4", driver="ESRI Shapefile")
  r_c<-crop(r, extent(continent_eck4))
  r_c<-mask(r_c, continent_eck4)
  plot(r_c)
  v<-values(r_c)
  
  values(r_c)[!is.na(v)]<-c(1:length(v[!is.na(v)]))
  plot(r_c)
  writeRaster(r_c, "../../Raster/mask_1km.tif", datatype="INT4U")
}

if (F){
  mask_bak<-raster("../../Raster/mask_10km.tif")
  mask_high<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/mask_high_res.tif")
  vessel <- sf::st_read(dsn = "/home/huijieqiao/Experiments/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
  #vessel_simpl <- st_simplify(vessel,
  #                            preserveTopology = TRUE, 
  #                            dTolerance = 10000)
  #object.size(vessel_simpl)
  #object.size(vessel)
  
  
  #plot_map(mv_simpl)
  
  sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_bak))
  #st_write(sp_df_eck4, "../../Data/Raw/IUCN/BIRDS/BIRDS_ECK4.gdb")
  bird_df<-data.frame(sp_df_eck4)
  bird_df$Shape<-NULL
  
  colnames(bird_df)
  bird_df<-data.table(bird_df)
  saveRDS(bird_df, "../../Data/Birds/bird_df.rda")
}
mask_10km<-raster("../../Raster/mask_10km.tif")
#mask_1km<-raster("../../Raster/mask_1km.tif")

min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2))
}

PRESENCE<-c(1,2,3,4,5)
ORIGIN<-c(1,2,3,5,6)
SEASONAL<-c(1,2)

if (F){
  vessel <- sf::st_read(dsn = "/home/huijieqiao/Experiments/IUCN/BIRDS/BOTW.gdb", layer = "All_Species")
  sp_df_eck4<-st_transform(vessel, crs = st_crs(mask_10km))
  
  bird_df<-readRDS("../../Data/Birds/bird_df.rda")
  unique <- unique(bird_df$SCINAME)
  i=10
  for (i in 1:length(unique)) {
    bi<-unique[i]
    print(paste(i, length(unique), bi))
    tmp_sf <- sp_df_eck4[which((bird_df$SCINAME == bi)&(bird_df$PRESENCE %in% PRESENCE)&
                                 (bird_df$ORIGIN %in% ORIGIN)&(bird_df$SEASONAL %in% SEASONAL)), ]
    saveRDS(tmp_sf, sprintf("../../Objects/IUCN_Distribution/Birds/RAW/%s.rda", gsub(" ", "_", bi)))
    tryCatch(
      {
        tmp_sf_sim<-st_simplify(tmp_sf, dTolerance = 5000)
        saveRDS(tmp_sf_sim, sprintf("../../Objects/IUCN_Distribution/Birds/st_simplify/%s.rda", gsub(" ", "_", bi)))
        tmp_sf_rm_sim<-rmapshaper::ms_simplify(tmp_sf, keep=0.01)
        saveRDS(tmp_sf_rm_sim, sprintf("../../Objects/IUCN_Distribution/Birds/ms_simplify/%s.rda", gsub(" ", "_", bi)))
      },
      error=function(cond) {
        
      }
    )
    
  }
}
