library(dggridR)
library(terra)
library(sf)
library(dplyr)

setwd("/media/huijieqiao/WD22T_50/Diversity_in_Future/diversity_in_e")
bioclim<-rast("/media/huijieqiao/WD22T_50/CMIP6/Bioclimatic_Alex/CNRM-ESM2-1/SSP119/1850/CNRM-ESM2-1_SSP119_1850_raw_based.tif")
dggs<-dgconstruct(res = 11)
dg_closest_res(dggs, col="spacing_km", val=20)
#Resolution: 11, Area (km^2): 287.933536398634, Spacing (km): 16.758963498128, CLS (km): 19.147021538141

earth<-dgearthgrid(dggs)
earth_fixed <- st_wrap_dateline(earth, options = c("WRAPDATELINE=YES"))

shpfname = "../Data/Shape/isea3h11.shp"
write_sf(earth_fixed, shpfname)
#neighbors<-find_connected_hexagon(earth)
shp<-read_sf("../Data/Shape/isea3h11.shp")
biolist<-rast(c(sprintf("../Data/Raster/wc2.1_10m_bio/wc2.1_10m_bio_%d.tif", c(1:19)),
                "../Data/Raster/wc2.1_10m_elev.tif"))