library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  df_10km<-readRDS(sprintf("../../Objects/Diversity_exposure_0_dispersal_1_10km_2_100km/Mammals/%s/diversity_df.rda", esm_ssp[1]))
  df_100km<-readRDS(sprintf("../../Objects/Diversity_exposure_0_dispersal_1/Mammals/%s/diversity_df.rda", esm_ssp[1]))
  
  sp_10km_2100<-names(df_10km[["2100"]])
  length(sp_10km_2100)
  sp_100km_2100<-unique(rbindlist(df_100km[["2100"]])$sp)
  length(sp_100km_2100)
  
  sp_100km_2100[!(sp_100km_2100 %in% sp_10km_2100)]
  
  df_10km_disp_1<-readRDS("../../Objects/Diversity_exposure_5_dispersal_1_10km_2_100km/Mammals/EC-Earth3-Veg_SSP119/diversity_df.rda")
  sp_10km_disp_1_2100<-names(df_10km_disp_1[["2100"]])
  length(sp_10km_disp_1_2100)
}


#d1<-readRDS("../../Objects/Dispersal/Mammals/Microryzomys_altissimus/UKESM1_SSP585_5_dispersal_1_10km.rda")
#d2<-readRDS("../../Objects/Dispersal/Mammals/Microryzomys_altissimus/UKESM1_SSP585_5_dispersal_1.rda")

esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
df_100km<-readRDS(sprintf("../../Objects/Diversity_exposure_5_dispersal_1/Mammals/%s/diversity_df.rda", esm_ssp[1]))

coms<-expand.grid(eee=esm_ssp, dispersal=c(0, 1), sp=unique(rbindlist(df_100km[["2100"]])$sp))
rm_files<-sprintf("../../Objects/Dispersal/Mammals/%s/%s_0_dispersal_%d_10km.rda", coms$sp, coms$eee, coms$dispersal)

skipped<-c()
for (i in 1:length(rm_files)){
  print(i)
  f<-rm_files[i]
  if (!file.exists(f)){
    next()
  }
  if (f %in% skipped){
    next()
  }
  dd<-readRDS(f)
  if (length(dd)==5){
    print(paste("removing", f, "filesize", file.size(f)))
  
    unlink(f)
  }else{
    skipped<-c(skipped, f)
  }
  
}
length(all)
