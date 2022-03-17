library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  sp_dis_se_list<-list()
  for (exposure in c(0, 5)){
    GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
    SSPs<-c("SSP119", "SSP245", "SSP585")
    
    predict_range<-c(2021:2100)
    layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
    layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
    
    
    #df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
    i=1
    j=4
    k=1
    #dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
    dispersals<-c(0:1)
    
    
    for (group in c("Mammals", "Birds")){
      for (j in c(1:nrow(layer_df))){
        layer<-layer_df[j,]
        for (k in c(1:length(dispersals))){
          folders<-c(sprintf("Diversity_exposure_%d_dispersal_%d_10km_2_100km", exposure, dispersals[k]))
          folder<-folders[1]
          for (folder in folders){
            layer$M<-dispersals[k]
            layer$TYPE<-folder
            print(paste("READING DATA", target_folder))
            target_folder<-sprintf("../../Objects/%s/%s/%s", folder, group, layer$LABEL)
            sp_dis<-data.table(readRDS(sprintf("%s/sp_dis.rda", target_folder)))
            sp_dis_se<-sp_dis[, .(YEAR=max(YEAR)), by=list(sp)]
            sp_dis_se<-sp_dis_se[YEAR==2025]
            sp_dis_se$group<-group
            sp_dis_se$gcm<-layer$GCM
            sp_dis_se$ssp<-layer$SSP
            sp_dis_se$exposure<-exposure
            sp_dis_se$dispersal<-dispersals[k]
            sp_dis_se_list[[length(sp_dis_se_list)+1]]<-sp_dis_se
          }
        }
      }
    }
  }

  sp_dis_se_list<-rbindlist(sp_dis_se_list)
  sp_dis_se_list[group=="Mammals"&exposure==5]
  sp_dis_se_list_x<-sp_dis_se_list[exposure==0]
  rm_files<-sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d_10km.rda", 
                    sp_dis_se_list_x$group, sp_dis_se_list_x$sp, 
                    sp_dis_se_list_x$gcm, sp_dis_se_list_x$ssp,
                    sp_dis_se_list_x$exposure, sp_dis_se_list_x$dispersal)
  sp_dis_se_list_x[sp=="Brachypteracias_leptosomus"]
  unlink(rm_files)
  sp_dis_se_list[, .(N=.N), by=list(group, gcm, ssp, exposure, dispersal)]
  
  
  df_10km<-readRDS(sprintf("../../Objects/Diversity_exposure_0_dispersal_1_10km_2_100km/Birds/%s/diversity_df.rda", esm_ssp[1]))
  df_100km<-readRDS(sprintf("../../Objects/Diversity_exposure_0_dispersal_1/Birds/%s/diversity_df.rda", esm_ssp[1]))
  
  sp_10km_2100<-names(df_10km[["2100"]])
  length(sp_10km_2100)
  sp_100km_2100<-unique(rbindlist(df_100km[["2100"]])$sp)
  length(sp_100km_2100)
  
  sp_100km_2100[!(sp_100km_2100 %in% sp_10km_2100)]
  
  df_10km_disp_1<-readRDS("../../Objects/Diversity_exposure_5_dispersal_1_10km_2_100km/Birds/EC-Earth3-Veg_SSP119/diversity_df.rda")
  sp_10km_disp_1_2100<-names(df_10km_disp_1[["2100"]])
  length(sp_10km_disp_1_2100)
}


#d1<-readRDS("../../Objects/Dispersal/Mammals/Microryzomys_altissimus/UKESM1_SSP585_5_dispersal_1_10km.rda")
#d2<-readRDS("../../Objects/Dispersal/Mammals/Microryzomys_altissimus/UKESM1_SSP585_5_dispersal_1.rda")

esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
df_100km<-readRDS(sprintf("../../Objects/Diversity_exposure_5_dispersal_1/Birds/%s/diversity_df.rda", esm_ssp[1]))

coms<-expand.grid(eee=esm_ssp, dispersal=c(0, 1), sp=unique(rbindlist(df_100km[["2100"]])$sp))
rm_files<-sprintf("../../Objects/Dispersal/Birds/%s/%s_0_dispersal_%d_10km.rda", coms$sp, coms$eee, coms$dispersal)

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
  
    #unlink(f)
  }else{
    skipped<-c(skipped, f)
  }
  
}
length(all)
