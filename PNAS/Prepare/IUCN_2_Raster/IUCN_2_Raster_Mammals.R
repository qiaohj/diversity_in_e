library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#for MAMMALS

if (F){
  PRESENCE<-c(1,2,3,4,5)
  ORIGIN<-c(1,2,3,5,6)
  SEASONAL<-c(1,2)
  sp_df<-sf::st_read(dsn="/home/huijieqiao/Experiments/IUCN/MAMMALS", layer="MAMMALS_TERRESTRIAL_ONLY") 
  mask_bak<-raster("../../Raster/mask_10km.tif")
  sp_df_eck4<-st_transform(sp_df, crs = st_crs(mask_bak))
  #writeOGR(sp_df_eck4, "../../Data/Raw/IUCN/MAMMALS", "MAMMALS_TERRESTRIAL_ONLY_ECK4", driver="ESRI Shapefile")
  
  mammals_df<-data.frame(sp_df_eck4)
  
  mammals_df$geometry<-NULL
  
  colnames(mammals_df)
  mammals_df<-data.table(mammals_df)
  
  saveRDS(mammals_df, "../../Data/Mammals/mammals_df.rda")
  
  unique <- unique(mammals_df$binomial)
  i=10
  for (i in 1:length(unique)) {
    bi<-unique[i]
    print(paste(i, length(unique), bi))
    tmp_sf <- sp_df_eck4[which((mammals_df$binomial == bi)&(mammals_df$presence %in% PRESENCE)&
                                 (mammals_df$origin %in% ORIGIN)&(mammals_df$seasonal %in% SEASONAL)), ]
    saveRDS(tmp_sf, sprintf("../../Objects/IUCN_Distribution/Mammals/RAW/%s.rda", gsub(" ", "_", bi)))
    tryCatch(
      {
        tmp_sf_sim<-st_simplify(tmp_sf, dTolerance = 5000)
        saveRDS(tmp_sf_sim, sprintf("../../Objects/IUCN_Distribution/Mammals/st_simplify/%s.rda", gsub(" ", "_", bi)))
        tmp_sf_rm_sim<-rmapshaper::ms_simplify(tmp_sf, keep=0.01)
        saveRDS(tmp_sf_rm_sim, sprintf("../../Objects/IUCN_Distribution/Mammals/ms_simplify/%s.rda", gsub(" ", "_", bi)))
      },
      error=function(cond) {
        
      }
    )
    
  }
}
