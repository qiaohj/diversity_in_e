setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(sf)
library(ggplot2)
r<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/bioclim/wc2.0_bio_30s_01.tif")
r2<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/bioclim/wc2.0_bio_30s_12.tif")
plot(r)
plot(r2)
bi<-"Poicephalus rufiventris"
tmp_sf<-readRDS(sprintf("../../Objects/IUCN_Distribution/Birds/st_simplify/%s.rda", gsub(" ", "_", bi)))
plot(st_geometry(tmp_sf))

no_disp_no_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_0_dispersal_0.rda")
no_disp_no_exposure<-no_disp_no_exposure[['2100']]
ggplot(no_disp_no_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[no_disp_no_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)


no_disp_with_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_5_dispersal_0.rda")
no_disp_with_exposure<-no_disp_with_exposure[['2100']]
ggplot(no_disp_with_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[no_disp_with_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)

with_disp_with_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_5_dispersal_1.rda")
with_disp_with_exposure<-with_disp_with_exposure[['2100']]
ggplot(with_disp_with_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[with_disp_with_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)


with_disp_no_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_0_dispersal_1.rda")
with_disp_no_exposure<-with_disp_no_exposure[['2100']]
ggplot(with_disp_no_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[with_disp_no_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)

with_disp_no_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_0_dispersal_1.rda")
with_disp_no_exposure<-with_disp_no_exposure[['2040']]
#ggplot(with_disp_no_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE, main=2040)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[with_disp_no_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)


with_disp_no_exposure<-readRDS("../../Objects/Dispersal/Birds/Poicephalus_rufiventris/UKESM1_SSP245_0_dispersal_1.rda")
with_disp_no_exposure<-with_disp_no_exposure[['2080']]
#ggplot(with_disp_no_exposure)+geom_tile(aes(x=x, y=y))
mask_item<-crop(mask_100km, extent(1000000,   5000000,  -3000000, 3500000))
plot(mask_item, col="lightgrey", legend=FALSE, main=2080)
plot(st_geometry(tmp_sf), add=T, border="blue")
suitable_area<-mask_100km
values(suitable_area)<-NA
values(suitable_area)[with_disp_no_exposure$mask_100km]<-1
suitable_area<-crop(suitable_area, extent(1000000,   5000000,  -3000000, 3500000))
plot(suitable_area, col=alpha("indianred1", 0.2), add=T, legend=FALSE)



