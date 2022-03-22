library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(data.table)
library(sf)
library(fasterize)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")





group_disp_birds<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
colnames(group_disp_birds)[5]<-"sp"
colnames(group_disp_birds)[4]<-"Migration"
group_disp_birds$group<-"Birds"
group_disp_mammals<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
colnames(group_disp_mammals)[1]<-"sp"
group_disp_mammals$group<-"Mammals"


group_full<-rbindlist(list(group_disp_birds, group_disp_mammals), fill = T)
group_full$BodyMass.Value<-NULL
group_full$log_body_mass<-NULL
cols<-c("group", "sp", "HWI", "Diet", "Migration", "body_mass", "estimated_disp", "ForStrat")
group_full<-group_full[, ..cols]
table(group_full$group)
group_full$is_migratory_bird<-F

group_full<-group_full[sp!=" NA"]

unique_sp<-group_full[, .(N=.N), by=list(group, sp)]
unique_sp_t<-unique_sp[N>1]
item1<-group_full[sp=="Cuculus saturatus"][1]
item2<-group_full[sp=="Corapipo leucorrhoa"][2]
item3<-group_full[sp=="Tangara cyanoptera"][2]
group_full<-group_full[!(sp %in% unique_sp_t$sp)]
group_full<-rbindlist(list(group_full, item1, item2, item3))
group_full$in_ebird<-F
group_full$res<-""
i<-2
for (i in 1:nrow(group_full)){
  print(paste(i, nrow(group_full)))
  item<-group_full[i]
  #f1<-sprintf("../../Objects/Dispersal/%s/%s/initial_disp_10km_exposure_0_dispersal_0.rda", 
  #            item$group, gsub(" ", "_", item$sp))
  f1<-list.files(sprintf("../../Objects/Dispersal/%s/%s/", item$group, gsub(" ", "_", item$sp)), pattern="initial_disp_10km_*")
  res<-""
  mig_bird<-F
  ebird_bird<-F
  if (length(f1)>1){
    res<-"10km"
  }else{
    f1<-list.files(sprintf("../../Objects/Dispersal/%s/%s/", item$group, gsub(" ", "_", item$sp)), pattern="initial_disp_ex*")
    if (length(f1)>1){
      res<-"100km"
    }
  }
  if (res!=""){
   
    f1<-sprintf("../../Objects/Dispersal/%s/%s/fit_seasonal_2.rda", 
                item$group, gsub(" ", "_", item$sp))
    if (file.exists(f1)){
      mig_bird<- T
    }
    
    f1<-sprintf("../../Objects/Dispersal/%s/%s/fit_ebird.rda", 
                item$group, gsub(" ", "_", item$sp))
    if (file.exists(f1)){
      ebird_bird<- T
    }
  }
  group_full[i]$in_ebird<-ebird_bird
  group_full[i]$is_migratory_bird<-mig_bird
  group_full[i]$res<-res
}


sp_v_bird<-data.table(readRDS("../../Objects/Diversity_exposure_0_dispersal_0_10km_2_100km/Birds/EC-Earth3-Veg_SSP119/sp_dis.rda"))
sp_v_bird<-sp_v_bird[YEAR==2020]
sp_v_mammal<-data.table(readRDS("../../Objects/Diversity_exposure_0_dispersal_0_10km_2_100km/Mammals/EC-Earth3-Veg_SSP119/sp_dis.rda"))
sp_v_mammal<-sp_v_mammal[YEAR==2020]
dddd<-readRDS("../../Objects/Diversity_exposure_0_dispersal_0/Mammals/EC-Earth3-Veg_SSP119/sp_dis.rda")
length(unique(dddd$sp))
sp_v_mammal[!(sp %in% unique(dddd$sp))]$sp
sp_list<-c(sp_v_bird$sp, sp_v_mammal$sp)
sp_list<-gsub("_", " ", sp_list)
group_full_filter<-group_full[sp %in% sp_list]
table(group_full_filter$in_ebird)
table(group_full_filter$group)
#for resolution
group_full_filter[, .(N=.N), by=list(group, res, is_migratory_bird)]
group_full_filter[, .(N=.N), by=list(group, res)]

write.csv(group_full_filter, "../../Objects/species_list.csv", row.names = F)

group_full_filter<-data.table(read.csv("../../Objects/species_list.csv", head=T, stringsAsFactors = F))

poly<-readRDS(sprintf("../../Objects/IUCN_Distribution/Mammals/RAW/%s.rda", "Abrocoma_budini"))

xmin<-st_bbox(st_geometry(poly))[1]-200000
xmax<-st_bbox(st_geometry(poly))[3]+200000
ymin<-st_bbox(st_geometry(poly))[2]-200000
ymax<-st_bbox(st_geometry(poly))[4]+200000

plot(st_geometry(poly), xlim=c(xmin, xmax),  ## with c()
     ylim=c(ymin, ymax))

par(mfrow=c(1,2))   
mask_10km<-raster("../../Raster/mask_10km.tif")
plot(mask_10km, xlim=c(xmin, xmax),  ## with c()
     ylim=c(ymin, ymax))
plot(st_geometry(poly), add=T)

mask_100km<-raster("../../Raster/mask_100km.tif")
plot(mask_100km, xlim=c(xmin, xmax),  ## with c()
     ylim=c(ymin, ymax))
plot(st_geometry(poly), add=T)
