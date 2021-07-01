library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
r_hotspot<-raster("../../Raster/hotspots_2016_1/hotspots_2016_eck4.hotspot_area.tif")
r_hotspot_all<-raster("../../Raster/hotspots_2016_1/hotspots_2016_eck4.all.tif")
year<-2100
group<-"Amphibians"
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
SSP<-SSPs[3]
GCM<-GCMs[3]
M<-1
for (SSP in SSPs){
  for (M in c(0:2)){
    r_stack_all<-NULL
    for (group in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      r_stack<-NULL
      for (GCM in GCMs){
        print(paste(SSP, M, group, GCM, year))
        r<-raster(sprintf("../../Figures/Diversity/%s/species.richness/%s_%s_%d_1/%d.tif", group, GCM, SSP, M, year))
        
        if (is.null(r_stack)){
          r_stack<-r
        }else{
          r_stack<-stack(r_stack, r)
        }
      }
      r_ave<-mean(r_stack)
      if (is.null(r_stack_all)){
        r_stack_all<-r_ave
      }else{
        r_stack_all<-stack(r_stack_all, r_ave)
      }
    }   
    richness<-sum(r_stack_all)
    richness_p<-data.frame(rasterToPoints(richness))
    richness_p$is_big_hotspot<-raster::extract(r_hotspot_all, richness_p[, c("x", "y")])
    richness_p$is_key_hotspot<-raster::extract(r_hotspot, richness_p[, c("x", "y")])
    richness_p$hotspot<-"OUT"
    richness_p[!is.na(richness_p$is_big_hotspot),]$hotspot<-"IN"
    richness_p[!is.na(richness_p$is_key_hotspot),]$hotspot<-"KEY"
    richness_se<-richness_p%>%dplyr::group_by(hotspot)%>%
      dplyr::summarise(mean_richness=mean(layer),
                       sd_richness=sd(layer))
    plot(richness)
    plot(r_hotspot_all, add=T, col="black", alpha=0.3)              
  }
}