library(data.table)
library(raster)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  pt<-data.table(rasterToPoints(raster("../../Raster/mask_100km.tif")))
  vocc_list<-list()
  for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
    for (ssp in c("SSP119", "SSP245", "SSP585")){
      p<-pt
      for (var in c(1,5,6,12,13,14)){
        bio<-raster(sprintf("../../Raster/VoCC_mask/%s/%s/bio%d_voccMag.tif", gcm, ssp, var))
        p_bio<-data.table(rasterToPoints(bio))
        p<-merge(p, p_bio, by=c("x", "y"))
      }
      p$gcm<-gcm
      p$ssp<-ssp
      vocc_list[[paste(ssp, gcm)]]<-p
    }
  }
  saveRDS(vocc_list, "../../Objects/vocc_list.rda")
}
vocc_list<-readRDS("../../Objects/vocc_list.rda")

exposure=0
g<-"Mammals"
for (exposure in c(0, 5)){
  for (g in c("Birds", "Mammals")){
    print(paste(g, exposure))
    df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s_10km.rda", exposure, g))
    df<-rbindlist(df)
    df_se<-df[, .(N=.N), by=list(x, y, mask, GCM, SSP, dispersal)]
    if (F){
      ggplot(df_se[GCM=="UKESM1"&SSP=="SSP585"&dispersal==0])+geom_tile(aes(x=x, y=y, color=N))
    }
  }
}
for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
  for (ssp in c("SSP119", "SSP245", "SSP585")){
    vocc_item<-vocc_list[[paste(ssp, gcm)]]
    
    item<-df_se[GCM==gcm&SSP==ssp&dispersal==0]
    item_with_vocc<-merge(item, vocc_item, by.x="mask", by.y="mask_100km")
    cor(item_with_vocc$N, item_with_vocc$bio1_voccMag)
    cor(item_with_vocc$N, item_with_vocc$bio12_voccMag)
  }
}

ratio_final<-data.table(readRDS("../../Figures/when_where_extinction_all/ratio_final_group_10km.rda"))
bird_property<-data.table(read.csv("../../Objects/Supp.Tables/bird_property_10km.csv", stringsAsFactors = F))
mammal_property<-data.table(read.csv("../../Objects/Supp.Tables/mammal_property_10km.csv", stringsAsFactors = F))
property<-rbindlist(list(bird_property, mammal_property), fill=T)

bird_extinct<-data.table(read.csv("../../Objects/Supp.Tables/bird_extinction_10km.csv", stringsAsFactors = F))
mammal_extinct<-data.table(read.csv("../../Objects/Supp.Tables/mammal_extinction_10km.csv", stringsAsFactors = F))
extinct<-rbindlist(list(bird_extinct, mammal_extinct), fill=T)
extinct<-extinct[!is.na(extinct_year)]
extinct$sp<-gsub(" ", "_", extinct$sp)
property$scaled_range_bio1_sd_min<-scale(c(property$range_bio1_sd_min, property$range_bio1_sd_max))[1:nrow(property)]
property$scaled_range_bio1_sd_max<-scale(c(property$range_bio1_sd_max, property$range_bio1_sd_min))[1:nrow(property)]
property$scaled_range_bio5_sd_min<-scale(c(property$range_bio5_sd_min, property$range_bio5_sd_max))[1:nrow(property)]
property$scaled_range_bio5_sd_max<-scale(c(property$range_bio5_sd_max, property$range_bio5_sd_min))[1:nrow(property)]
property$scaled_range_bio6_sd_min<-scale(c(property$range_bio6_sd_min, property$range_bio6_sd_max))[1:nrow(property)]
property$scaled_range_bio6_sd_max<-scale(c(property$range_bio6_sd_max, property$range_bio6_sd_min))[1:nrow(property)]
property$scaled_range_bio12_sd_min<-scale(c(property$range_bio12_sd_min, property$range_bio12_sd_max))[1:nrow(property)]
property$scaled_range_bio12_sd_max<-scale(c(property$range_bio12_sd_max, property$range_bio12_sd_min))[1:nrow(property)]
property$scaled_range_bio13_sd_min<-scale(c(property$range_bio13_sd_min, property$range_bio13_sd_max))[1:nrow(property)]
property$scaled_range_bio13_sd_max<-scale(c(property$range_bio13_sd_max, property$range_bio13_sd_min))[1:nrow(property)]
property$scaled_range_bio14_sd_min<-scale(c(property$range_bio14_sd_min, property$range_bio14_sd_max))[1:nrow(property)]
property$scaled_range_bio14_sd_max<-scale(c(property$range_bio14_sd_max, property$range_bio14_sd_min))[1:nrow(property)]


property$scaled_nb_bio1_sd<-property$scaled_range_bio1_sd_max - property$scaled_range_bio1_sd_min
property$scaled_nb_bio5_sd<-property$scaled_range_bio5_sd_max - property$scaled_range_bio5_sd_min
property$scaled_nb_bio6_sd<-property$scaled_range_bio6_sd_max - property$scaled_range_bio6_sd_min
property$scaled_nb_bio12_sd<-property$scaled_range_bio12_sd_max - property$scaled_range_bio12_sd_min
property$scaled_nb_bio13_sd<-property$scaled_range_bio13_sd_max - property$scaled_range_bio13_sd_min
property$scaled_nb_bio14_sd<-property$scaled_range_bio14_sd_max - property$scaled_range_bio14_sd_min

property$nb_volume<-property$scaled_nb_bio1_sd*property$scaled_nb_bio5_sd*property$scaled_nb_bio6_sd*
  property$scaled_nb_bio12_sd*property$scaled_nb_bio13_sd*property$scaled_nb_bio14_sd
disp_2020_vocc_sp_full_list<-list()
models<-list()

for (group in c("Mammals", "Birds")){
  for (dispersal in c(0, 1)){
    for (exposure in c(0, 5)){
      for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
        for (ssp in c("SSP119", "SSP245", "SSP585")){
          f<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d_10km_2_100km/%s/%s_%s/diversity_df.rda", 
                     exposure, dispersal, group, gcm, ssp)
          print(f)
          disp<-readRDS(f)
          disp_2020<-rbindlist(disp[["2020"]])
          colnames(disp_2020)[c(1,2)]<-c("x", "y")
          disp_2020_vocc<-merge(disp_2020, vocc_list[[paste(ssp, gcm)]], by=c("x", "y", "mask_100km"))
          disp_2020_vocc_sp<-disp_2020_vocc[, .(bio1_voccMag=mean(bio1_voccMag),
                                                bio5_voccMag=mean(bio5_voccMag),
                                                bio6_voccMag=mean(bio6_voccMag),
                                                bio12_voccMag=mean(bio12_voccMag),
                                                bio13_voccMag=mean(bio13_voccMag),
                                                bio14_voccMag=mean(bio14_voccMag)),
                                            by=list(sp)]
          property$sp<-gsub(" ", "_", property$sp)
          extinct_item<-extinct[SSP==ssp&ESM==gcm]
          disp_2020_vocc_sp_full<-merge(disp_2020_vocc_sp, property, by="sp")
          disp_2020_vocc_sp_full$extinct<-0
          disp_2020_vocc_sp_full[sp %in% extinct_item$sp]$extinct<-1
          disp_2020_vocc_sp_full$GCM<-gcm
          disp_2020_vocc_sp_full$SSP<-ssp
          disp_2020_vocc_sp_full$dispersal<-dispersal
          disp_2020_vocc_sp_full$exposure<-exposure
          disp_2020_vocc_sp_full_list[[length(disp_2020_vocc_sp_full_list)+1]]<-disp_2020_vocc_sp_full
          m<-lm(extinct~nb_volume+estimated_disp+bio1_voccMag, data=disp_2020_vocc_sp_full)
          summ<-summary(m)
          r2<-summ$r.squared
          p_nb_volume<-summ$coefficients[2, 4]
          p_estimated_disp<-summ$coefficients[3, 4]
          p_bio1_voccMag<-summ$coefficients[4, 4]
          models[[length(models)+1]]<-data.frame(GCM=gcm, SSP=ssp, dispersal=dispersal, exposure=exposure,
                                                 group=group, r2=r2, p_nb_volume=p_nb_volume, p_estimated_disp=p_estimated_disp,
                                                 p_bio1_voccMag=p_bio1_voccMag)
        }
      }
    }
  }
}

models<-rbindlist(models)
saveRDS(models, "../../Figures/VoCC/lm_model.rda")
disp_2020_vocc_sp_full_list_all<-rbindlist(disp_2020_vocc_sp_full_list)
saveRDS(disp_2020_vocc_sp_full_list_all, "../../Figures/VoCC/sp_vocc.rda")

models$nb_volume_sig<-ifelse(models$p_nb_volume<=0.05, T, F)
models$estimated_disp_sig<-ifelse(models$p_estimated_disp<=0.05, T, F)
models$bio1_voccMag_sig<-ifelse(models$p_bio1_voccMag<=0.05, T, F)
models_se_nb<-models[, .(N=.N), by=list(nb_volume_sig, SSP)]
models_se_disp<-models[dispersal==1, .(N=.N), by=list(estimated_disp_sig, SSP, group)]
models_se_vocc<-models[, .(N=.N), by=list(bio1_voccMag_sig, SSP, group)]
