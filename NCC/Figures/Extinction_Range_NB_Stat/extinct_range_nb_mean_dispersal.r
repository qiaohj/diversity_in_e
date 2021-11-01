library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")



exposure<-c(0, 5)
SSPs<-c("SSP119", "SSP245", "SSP585")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")

combs<-expand.grid(SSP=SSPs, GCM=GCMs, exposure=exposure)

if (T){
  bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  if (F){
    write.csv(bird_disp, "../../Objects/estimate_disp_dist/estimate_disp_dist_bird.csv", row.names=F)
  }
  tempelate<-"../../Objects/Dispersal/Birds/%s/%s_%s_%d_dispersal_1.rda"
  i_e=1
  i=9
  bird_disp<-bird_disp[sample(nrow(bird_disp), nrow(bird_disp)),]
  for (i in c(1:nrow(bird_disp))){
    
    sp<-bird_disp[i, ]$iucn_name
    print(paste(i, nrow(bird_disp), sp))
    target<-sprintf("../../Objects/Dispersal_distances/Birds/%s.rda", gsub(" ", "_", sp))
    target_info<-sprintf("../../Objects/Dispersal_distances/Birds/%s_info.rda", gsub(" ", "_", sp))
    if (file.exists(target)){
      next()
    }
    saveRDS(NULL, target)
    all_item<-list()
    all_item_info<-list()
    for (i_e in c(1:nrow(combs))){
      com_item<-combs[i_e,]
      rda_file<-sprintf(tempelate, gsub(" ", "_", sp), com_item$GCM, com_item$SSP, com_item$exposure)
      if (!file.exists(rda_file)){
        next()
      }
      item<-readRDS(rda_file)
      year=2039
      
      for (year in c(2022:2100)){
        
        y1<-year-1
        y2<-year  
        item1<-item[[as.character(y1)]]
        item2<-item[[as.character(y2)]]
        if (is.null(item2)){
          break()
        }
        me_item<-merge(item1, item2, by=c("x", "y", "mask_100km"))
        me_item$dispersal_dist<-me_item$accumulative_disp.y-me_item$accumulative_disp.x
        me_item[accumulative_disp.x>50000]$dispersal_dist<-me_item[accumulative_disp.x>50000]$accumulative_disp.y
        N_Reset_Cell<-nrow(me_item[accumulative_disp.x>50000])
        N_New_Cell<-nrow(item2[!(mask_100km %in% item1$mask_100km)])
        #me_item<-me_item[dispersal_dist>0]
        cols<-c("x", "y", "mask_100km",  "dispersal_dist", "accumulative_disp.x", "accumulative_disp.y")
        me_item<-me_item[, ..cols]
        me_item$year<-y2
        me_item$iucn_name<-sp
        me_item$SSP<-com_item$SSP
        me_item$GCM<-com_item$GCM
        me_item$exposure<-com_item$exposure
        info_item<-com_item
        info_item$year<-y2
        info_item$iucn_name<-sp
        info_item$N_Reset_Cell<-N_Reset_Cell
        info_item$N_New_Cell<-N_New_Cell
        all_item[[paste(sp, com_item$GCM, com_item$SSP, com_item$exposure, year)]]<-me_item
        all_item_info[[paste(sp, com_item$GCM, com_item$SSP, com_item$exposure, year)]]<-info_item
      }
    }
    all_item<-rbindlist(all_item)
    all_item_info<-rbindlist(all_item_info)
    saveRDS(all_item, target)
    saveRDS(all_item_info, target_info)
  }
}
#all_item<-rbindlist(all_item)
#all_item_se<-all_item[, .(mean_dist=mean(dispersal_dist), sd_dist=sd(dispersal_dist)), 
#                      by=c("iucn_name", "year", "exposure")]
if (T){
  mammal_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  if (F){
    write.csv(mammal_disp, "../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.csv", row.names=F)
  }
  tempelate<-"../../Objects/Dispersal/Mammals/%s/%s_%s_%d_dispersal_1.rda"
  i_e=1
  i=9
  mammal_disp<-mammal_disp[sample(nrow(mammal_disp), nrow(mammal_disp)),]
  sp<-"Caluromysiops irrupta"
  for (i in c(1:nrow(mammal_disp))){
    
    sp<-mammal_disp[i,]$Scientific
    #sp<-"Acerodon_leucotis"
    print(paste(i, nrow(mammal_disp), sp))
    target<-sprintf("../../Objects/Dispersal_distances/Mammals/%s.rda", gsub(" ", "_", sp))
    target_info<-sprintf("../../Objects/Dispersal_distances/Mammals/%s_info.rda", gsub(" ", "_", sp))
    if (file.exists(target)){
      next()
    }
    saveRDS(NULL, target)
    all_item<-list()
    all_item_info<-list()
    for (i_e in c(1:nrow(combs))){
      com_item<-combs[i_e,]
      rda_file<-sprintf(tempelate, gsub(" ", "_", sp), com_item$GCM, com_item$SSP, com_item$exposure)
      if (!file.exists(rda_file)){
        next()
      }
      item<-readRDS(rda_file)
      year=2039
      
      for (year in c(2022:2100)){
        
        y1<-year-1
        y2<-year  
        item1<-item[[as.character(y1)]]
        item2<-item[[as.character(y2)]]
        if (is.null(item2)){
          break()
        }
        me_item<-merge(item1, item2, by=c("x", "y", "mask_100km"))
        me_item$dispersal_dist<-me_item$accumulative_disp.y-me_item$accumulative_disp.x
        me_item[accumulative_disp.x>50000]$dispersal_dist<-me_item[accumulative_disp.x>50000]$accumulative_disp.y
        N_Reset_Cell<-nrow(me_item[accumulative_disp.x>50000])
        N_New_Cell<-nrow(item2[!(mask_100km %in% item1$mask_100km)])
        #me_item<-me_item[dispersal_dist>0]
        cols<-c("x", "y", "mask_100km",  "dispersal_dist", "accumulative_disp.x", "accumulative_disp.y")
        me_item<-me_item[, ..cols]
        me_item$year<-y2
        me_item$iucn_name<-sp
        me_item$SSP<-com_item$SSP
        me_item$GCM<-com_item$GCM
        me_item$exposure<-com_item$exposure
        info_item<-com_item
        info_item$year<-y2
        info_item$iucn_name<-sp
        info_item$N_Reset_Cell<-N_Reset_Cell
        info_item$N_New_Cell<-N_New_Cell
        all_item[[paste(sp, com_item$GCM, com_item$SSP, com_item$exposure, year)]]<-me_item
        all_item_info[[paste(sp, com_item$GCM, com_item$SSP, com_item$exposure, year)]]<-info_item
      }
    }
    all_item<-rbindlist(all_item)
    all_item_info<-rbindlist(all_item_info)
    saveRDS(all_item, target)
    saveRDS(all_item_info, target_info)
  }
}

if (F){
  rm(list=ls())
  bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  tempelate<-"../../Objects/Dispersal/Birds/%s/%s_%s_%d_dispersal_1.rda"
  i_e=1
  i=9
  bird_disp<-bird_disp[sample(nrow(bird_disp), nrow(bird_disp)),]
  all_df<-list()
  for (i in c(1:nrow(bird_disp))){
    
    sp<-bird_disp[i, ]$iucn_name
    print(paste(i, nrow(bird_disp), sp))
    target<-sprintf("../../Objects/Dispersal_distances/Birds/%s.rda", gsub(" ", "_", sp))
    target_info<-sprintf("../../Objects/Dispersal_distances/Birds/%s_info.rda", gsub(" ", "_", sp))
    item<-readRDS(target)
    if (is.null(item)){
      next()
    }
    if (nrow(item)==0){
      next()
    }
    item_info<-readRDS(target_info)
    item_info_se<-item_info[, .(N_Reset_Cell=sum(N_Reset_Cell),
                                N_New_Cell=sum(N_New_Cell),
                                mean_N_New_Cell=mean(N_New_Cell),
                                mean_N_Reset_Cell=mean(N_Reset_Cell)),
                            by=c("iucn_name", "exposure", "SSP", "GCM")]
    all_item_se<-item[, .(mean_dist=mean(dispersal_dist), sd_dist=sd(dispersal_dist)), 
                      by=c("iucn_name", "exposure", "SSP", "GCM")]
    all_item_se<-merge(all_item_se, item_info_se, by=c("iucn_name", "exposure", "SSP", "GCM"))
    all_df[[sp]]<-all_item_se
    
  }
  all_df<-rbindlist(all_df)
  
  mammal_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  
  i_e=1
  i=9
  
  all_df_2<-list()
  for (i in c(1:nrow(mammal_disp))){
    
    sp<-mammal_disp[i,]$Scientific
    print(paste(i, nrow(mammal_disp), sp))
    target<-sprintf("../../Objects/Dispersal_distances/Mammals/%s.rda", gsub(" ", "_", sp))
    target_info<-sprintf("../../Objects/Dispersal_distances/Mammals/%s_info.rda", gsub(" ", "_", sp))
    item<-readRDS(target)
    if (is.null(item)){
      next()
    }
    if (nrow(item)==0){
      next()
    }
    item_info<-readRDS(target_info)
    item_info_se<-item_info[, .(N_Reset_Cell=sum(N_Reset_Cell),
                                N_New_Cell=sum(N_New_Cell),
                                mean_N_New_Cell=mean(N_New_Cell),
                                mean_N_Reset_Cell=mean(N_Reset_Cell)),
                            by=c("iucn_name", "exposure", "SSP", "GCM")]
    all_item_se<-item[, .(mean_dist=mean(dispersal_dist), sd_dist=sd(dispersal_dist)), 
                      by=c("iucn_name", "exposure", "SSP", "GCM")]
    all_item_se<-merge(all_item_se, item_info_se, by=c("iucn_name", "exposure", "SSP", "GCM"))
    all_df_2[[sp]]<-all_item_se
  }
  all_df_2<-rbindlist(all_df_2)
  
  all_df$group<-"Birds"
  all_df_2$group<-"Mammals"
  all_df_2[iucn_name=="Abditomys latidens"]
  
  all_df_all<-rbindlist(list(all_df, all_df_2))
  saveRDS(all_df_all, "../../Objects/Dispersal_distances/all_mean_disp_dist.rda")
}

all_df_all<-readRDS("../../Objects/Dispersal_distances/all_mean_disp_dist.rda")
all_df_all[iucn_name=="Caluromysiops irrupta"]
