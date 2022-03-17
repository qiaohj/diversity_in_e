library(sf)
library(raster)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_10km.tif")
if (F){
  mountain<-st_read("../../Shape/GMBA mountain inventory V1.2/GMBA Mountain Inventory_v1.2-World.shp")
  mountain<-st_transform(mountain, crs=st_crs(mask))
  
  all_info<-readRDS("../../Objects/Island/species_mountain_10km.rda")
  mountain_species<-all_info[mountain>0]
  i=1
  for (i in c(1:nrow(mountain_species))){
    print(paste(i, nrow(mountain_species)))
    sp_item<-mountain_species[i]
    sp<-sp_item$sp
    target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", sp_item$group, sp)
    if (!file.exists(target)){
      next()
    }
    dis<-readRDS(target)
    if (nrow(dis)==0){
      next()
    }
    extend<-st_bbox(dis)
    extend<-c(extend[1], extend[3], extend[2], extend[4])
    
    if (between(extend[1], extant_mask[1], extant_mask[2])&
        between(extend[2], extant_mask[1], extant_mask[2])&
        between(extend[3], extant_mask[3], extant_mask[4])&
        between(extend[4], extant_mask[3], extant_mask[4])){
      dis_c<-crop(mountain, extend)
      dis_c<-mask(dis_c, dis)
      v<-values(dis_c)
      v_plain<-length(v[which(v==0)])
      v_mountain<-length(v[which(v==1)])
      if (F){
        plot(st_geometry(dis))
        plot(dis_c, add=T)
      }
      item<-data.frame(sp=sp, plain=v_plain, mountain=v_mountain, group="Mammals")
      all_info[[sp]]<-item
    }else{
      print("no overlap, skip")
    }
  }
}
