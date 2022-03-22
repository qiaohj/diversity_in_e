library(data.table)
library(raster)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
p<-data.table(rasterToPoints(raster("../../Raster/mask_100km.tif")))
for (var in c(1,5,6,12,13,14)){
  bio<-raster(sprintf("../../Raster/VoCC_mask/%s/%s/bio%d_voccMag.tif", gcm, ssp, var))
  p_bio<-data.table(rasterToPoints(bio))
  p<-merge(p, p_bio, by=c("x", "y"))
}
for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
  for (ssp in c("SSP119", "SSP245", "SSP585")){
    
    for (group in c("Birds", "Mammals")){
      for (dispersal in c(0, 1)){
        for (exposure in c(0, 5)){
          f_str<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d_10km_2_100km/%s/%s_%s/sp_dis.rda", 
                         exposure, dispersal, group, gcm, ssp)
          div_str<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d_10km_2_100km/%s/%s_%s/diversity_df.rda", 
                         exposure, dispersal, group, gcm, ssp)
          div<-readRDS(div_str)
          div_2020<-rbindlist(div[["2020"]])
          
          dis<-data.table(readRDS(f_str))
          dis_2100<-dis[YEAR==2100]
          dis_2020<-dis[YEAR==2020]
          dis_2020$is_extinct<-!(dis_2020$sp %in% dis_2100$sp)
          dis_2020_vocc<-merge(div_2020, p, by=c("mask_100km"))
          dis_2020_vocc_se<-dis_2020_vocc[, .(bio1_voccMag=mean(bio1_voccMag),
                                              sd_bio1_voccMag=sd(bio1_voccMag),
                                              bio5_voccMag=mean(bio5_voccMag),
                                              sd_bio5_voccMag=sd(bio5_voccMag),
                                              bio6_voccMag=mean(bio6_voccMag),
                                              sd_bio6_voccMag=sd(bio6_voccMag),
                                              bio12_voccMag=mean(bio12_voccMag),
                                              sd_bio12_voccMag=sd(bio12_voccMag),
                                              bio13_voccMag=mean(bio13_voccMag),
                                              sd_bio13_voccMag=sd(bio13_voccMag),
                                              bio14_voccMag=mean(bio14_voccMag),
                                              sd_bio14_voccMag=sd(bio14_voccMag)),
                                          by=list(sp, res)]
          dis_2020_vocc_se<-merge(dis_2020_vocc_se, dis_2020, by=c("sp"))
          if (group=="Mammals"){
            traits<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
            colnames(traits)[1]<-"sp"
            traits$sp<-gsub(" ", "_", traits$sp)
          }
          dis_2020_vocc_se<-merge(dis_2020_vocc_se, traits, by="sp")
          niche<-readRDS(sprintf("../../Objects/Species_property/%s_property.rda", group))
          
          
          niche$scaled_range_bio1_sd_max<-scale(c(niche$range_bio1_sd_max, 
                                          niche$range_bio1_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio1_sd_min<-scale(c(niche$range_bio1_sd_min, 
                                                   niche$range_bio1_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio1<-niche$scaled_range_bio1_sd_max-niche$scaled_range_bio1_sd_min
          
          niche$scaled_range_bio5_sd_max<-scale(c(niche$range_bio5_sd_max, 
                                                   niche$range_bio5_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio5_sd_min<-scale(c(niche$range_bio5_sd_min, 
                                                   niche$range_bio5_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio5<-niche$scaled_range_bio5_sd_max-niche$scaled_range_bio5_sd_min
          
          niche$scaled_range_bio6_sd_max<-scale(c(niche$range_bio6_sd_max, 
                                                   niche$range_bio6_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio6_sd_min<-scale(c(niche$range_bio6_sd_min, 
                                                   niche$range_bio6_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio6<-niche$scaled_range_bio6_sd_max-niche$scaled_range_bio6_sd_min
          
          niche$scaled_range_bio12_sd_max<-scale(c(niche$range_bio12_sd_max, 
                                                   niche$range_bio12_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio12_sd_min<-scale(c(niche$range_bio12_sd_min, 
                                                   niche$range_bio12_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio12<-niche$scaled_range_bio12_sd_max-niche$scaled_range_bio12_sd_min
          
          niche$scaled_range_bio13_sd_max<-scale(c(niche$range_bio13_sd_max, 
                                                   niche$range_bio13_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio13_sd_min<-scale(c(niche$range_bio13_sd_min, 
                                                   niche$range_bio13_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio13<-niche$scaled_range_bio13_sd_max-niche$scaled_range_bio13_sd_min
         
          niche$scaled_range_bio14_sd_max<-scale(c(niche$range_bio14_sd_max, 
                                                   niche$range_bio14_sd_min))[1:nrow(niche)]
          niche$scaled_range_bio14_sd_min<-scale(c(niche$range_bio14_sd_min, 
                                                   niche$range_bio14_sd_max))[1:nrow(niche)]
          niche$scaled_range_SD_bio14<-niche$scaled_range_bio14_sd_max-niche$scaled_range_bio14_sd_min
          niche$sd_nb_volume<-niche$scaled_range_SD_bio1*niche$scaled_range_SD_bio5*niche$scaled_range_SD_bio6*
            niche$scaled_range_SD_bio12*niche$scaled_range_SD_bio13*niche$scaled_range_SD_bio14
          
          dis_2020_vocc_se_niche<-merge(dis_2020_vocc_se, niche, by="sp")
          dis_2020_vocc_se_niche$nb_cuts<-as.numeric(
            Hmisc::cut2(dis_2020_vocc_se_niche$sd_nb_volume,
                        cuts=quantile(dis_2020_vocc_se_niche$sd_nb_volume, seq(0, 1, by=0.1))))
          dis_2020_vocc_se_niche$disp_cuts<-as.numeric(
            Hmisc::cut2(dis_2020_vocc_se_niche$estimated_disp,
                        cuts=quantile(dis_2020_vocc_se_niche$estimated_disp, seq(0, 1, by=0.1))))
          dis_2020_vocc_se_niche$cuts<-paste(dis_2020_vocc_se_niche$nb_cuts, dis_2020_vocc_se_niche$disp_cuts)
          
          dis_2020_vocc_se_niche$is
          
          ggplot(dis_2020_vocc_se_niche)+geom_boxplot(aes(y=bio1_voccMag, x=cuts, color=is_extinct))+
            facet_wrap(~island)
          
        }
      }
    }
  }
}


