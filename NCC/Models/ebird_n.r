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
setDTthreads(1)
print(sprintf("Current core number is %d", getDTthreads()))
mask_10km<-raster("../../Raster/mask_10km.tif")
bird_df<-readRDS("../../Data/Birds/bird_df.rda")
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
bird_full<-merge(bird_df, bird_disp, by.x="SCINAME", by.y="iucn_name", all=F)

unique <- unique(bird_full$SCINAME)
unique<-as.character(unique)

bird_full_sum_area<-bird_full[, .(sum_are=sum(Shape_Area)), by="SCINAME"]
bird_full_sum_area<-bird_full_sum_area[order(sum_are),]
bi<-bird_full_sum_area[bird_full_sum_area$sum_are<=1.5*min(bird_full_sum_area$sum_are)]$SCINAME[1]
bird_full_sum_area<-bird_full_sum_area[sample(nrow(bird_full_sum_area), nrow(bird_full_sum_area))]
ebird_base<-"/media/huijieqiao/SSD_Fast/ES50_eBird/Tables/Speccies_Split_202112"
i=2
for (i in 7731:length(bird_full_sum_area$SCINAME)) {
  
  bi<-bird_full_sum_area$SCINAME[i]
  #bi="Pseudophryne occidentalis"
  print(paste(i, length(unique), bi))
  ebird_full<-sprintf("%s/%s", ebird_base, bi)
  if (!file.exists(ebird_full)){
    next()
  }
  target<-sprintf("../../Objects/Dispersal/Birds/%s/n_ebird.rda", gsub(" ", "_", bi))
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  ff<-""
  cols_ll<-c("SCIENTIFIC_NAME", "LATITUDE", "LONGITUDE", "YEAR")
  all_items<-list()
  for (ff in list.files(ebird_full, pattern="\\.rda")){
    item<-readRDS(sprintf("%s/%s/%s", ebird_base, bi, ff))
    item<-item[APPROVED==1|REVIEWED==1]
    item<-item[, ..cols_ll]
    item<-unique(item)
    if (nrow(item)>0){
      all_items[[ff]]<-item
    }
  }
  all_items<-rbindlist(all_items)
  cols_ll<-c("LONGITUDE", "LATITUDE")
  points<-st_as_sf(x=all_items, coords = cols_ll, 
                   crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  points<-st_transform(points, crs=st_crs(mask_10km))
  all_items$mask_10km<-raster::extract(mask_10km, st_coordinates(points))
  all_items_n<-all_items[, .(N=.N), by=list(mask_10km, YEAR)]
  
  print("saving result")
  saveRDS(all_items_n, target)
  
  
  
  #saveRDS(p_buffer, sprintf("%s/p_buffer.rda", target_folder))
}
