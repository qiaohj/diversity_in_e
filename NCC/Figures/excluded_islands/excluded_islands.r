library(sf)
library(raster)
library(dplyr)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#the island layer was downloaded via https://rmgsc.cr.usgs.gov/gie/
if (F){
  islands<-st_read("../../Shape/USGSEsriWCMC_GlobalIslands_v3/v10/globalislandsfix.gdb")
  
  mask_100km<-raster("../../Raster/mask_100km.tif")
  
  mask_10km<-raster("../../Raster/mask_10km.tif")
  islands<-st_transform(islands, crs=st_crs(mask_10km))
  i=4435
  extend_100km<-extent(mask_100km)
  extend_10km<-extent(mask_10km)
  all<-list()
  
  for (i in c(1:nrow(islands))){
    print(paste(i, nrow(islands)))
    f<-islands[i,]
    extend<-st_bbox(f)
    extend<-c(extend[1], extend[3], extend[2], extend[4])
    
    if (between(extend[1], extend_100km[1], extend_100km[2])&
        between(extend[2], extend_100km[1], extend_100km[2])&
        between(extend[3], extend_100km[3], extend_100km[4])&
        between(extend[4], extend_100km[3], extend_100km[4])){
      m_100km<-crop(mask_100km, extend)
      v_100km<-values(m_100km)
      v_100km<-v_100km[!is.na(v_100km)]
      N_100km<-length(v_100km)
      
      if (N_100km>1){
        m_100km<-mask(m_100km, f)
        v_100km<-values(m_100km)
        v_100km<-v_100km[!is.na(v_100km)]
        N_100km<-length(v_100km)
      }
    }else{
      print("no overlap in 100km")
      N_100km<-0
    }
    
    if (between(extend[1], extend_10km[1], extend_10km[2])&
        between(extend[2], extend_10km[1], extend_10km[2])&
        between(extend[3], extend_10km[3], extend_10km[4])&
        between(extend[4], extend_10km[3], extend_10km[4])){
      m_10km<-crop(mask_10km, extend)
      v_10km<-values(m_10km)
      v_10km<-v_10km[!is.na(v_10km)]
      N_10km<-length(v_10km)
      if (N_10km>1){
        m_10km<-mask(m_10km, f)
        v_10km<-values(m_10km)
        v_10km<-v_10km[!is.na(v_10km)]
        N_10km<-length(v_10km)
      }
    }else{
      print("no overlap in 10km")
      N_10km<-0
    }
    
    item<-data.frame(Name_USGSO=f$Name_USGSO, IslandArea=f$IslandArea, Plate=f$Plate,
                     USGS_ISID=f$USGS_ISID, NEAR_FID=f$NEAR_FID, N_100km=N_100km, N_10km=N_10km)
    all[[length(all)+1]]<-item
  }
  all2<-rbindlist(all)
  saveRDS(all2, "../../Figures/excluded_islands/excluded_islands.rda")
}

all2<-readRDS("../../Figures/excluded_islands/excluded_islands.rda")

all_10km<-all2[N_10km==0]
all_10km<-all_10km[trimws(Name_USGSO)!=""]

islands<-st_read("../../Shape/USGSEsriWCMC_GlobalIslands_v3/v10/globalislandsfix.gdb")
islands_missed<-merge(islands, all_10km, by=c("Name_USGSO", "IslandArea", "Plate", "NEAR_FID", "USGS_ISID"))
islands_missed<-islands_missed[which(islands_missed$Name_USGSO!="Greenland"),]
st_write(islands_missed, "../../Shape/excluded_island/excluded_island.shp", append=F)
quantile(islands_missed$IslandArea, c(0.99))
islands_missed_top<-islands_missed[which(islands_missed$IslandArea>=100),]
st_write(islands_missed_top, "../../Shape/excluded_island/excluded_island_top.shp", append=F)
hist(islands_missed$IslandArea)
max_area<-max(islands_missed$IslandArea)
islands_missed[which(islands_missed$IslandArea>1000),]

mask<-raster("../../Raster/mask_10km.tif")
islands_missed<-st_transform(islands_missed, crs=st_crs(mask))
plot(mask, col="lightgrey")
plot(st_geometry(islands_missed[which(islands_missed$IslandArea>100),]), add=T, col="black", border="black")
plot(st_geometry(islands_missed[which(islands_missed$IslandArea>500),]), add=T, col="blue", border="blue")
plot(st_geometry(islands_missed[which(islands_missed$IslandArea>1000),]), add=T, col="red", border="red")

length(unique(islands_missed$Name_USGSO))

unique(islands_missed[which(islands_missed$IslandArea>100),]$Name_USGSO)
cols<-c("Name_USGSO", "IslandArea")

large_islands<-unique(data.frame(islands_missed[which(islands_missed$IslandArea>100), cols]))
large_islands<-large_islands[, cols]
write.csv(large_islands, "../../Shape/excluded_island_100km.csv", row.names = F)


large_islands<-unique(data.frame(islands_missed[which(islands_missed$IslandArea>200), cols]))
large_islands<-large_islands[, cols]
write.csv(large_islands, "../../Shape/excluded_island_200km.csv", row.names = F)
unique(islands_missed[which(islands_missed$IslandArea>500),]$Name_USGSO)
unique(islands_missed[which(islands_missed$IslandArea>1000),]$Name_USGSO)
