library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
exposures<-c(0, 5)

predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, exposure=exposures, group=c("Mammals", "Birds"), stringsAsFactors = F)
i<-1
full_list<-list()
for (i in c(1:nrow(layer_df))){
  itemc<-layer_df[i,]
  f<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_1_10km_2_100km/%s/%s_%s/diversity_df.rda",
             itemc$exposure, itemc$group, itemc$GCM, itemc$SSP)
  print(f)
  df<-readRDS(f)
  all_item<-list()
  for (year in names(df)){
    item<-rbindlist(df[[year]])
    all_item[[year]]<-item
  }
  all_item<-rbindlist(all_item)
  all_item_se<-all_item[, .(max_y=max(y_100km), 
                            min_y=min(y_100km),
                            max_abs_y=max(abs(y_100km)), 
                            min_abs_y=min(abs(y_100km)),
                            quantile_5=quantile(y_100km, 0.05),
                            quantile_95=quantile(y_100km, 0.95),
                            quantile_abs_5=quantile(abs(y_100km), 0.05),
                            quantile_abs_95=quantile(abs(y_100km), 0.95)),
                        by=list(sp, res, YEAR)]
  all_item_last_year<-all_item_se[, .(last_year=max(YEAR)), by=list(sp, res)]
  all_item_se_last_year<-merge(all_item_se, all_item_last_year, by=c("sp", "res"))
  all_item_se_last_last_year<-all_item_se_last_year[YEAR==last_year]
  colnames(all_item_se_last_last_year)[c(4:11)]<-sprintf("last_year_%s", colnames(all_item_se_last_last_year)[c(4:11)])
  all_item_se_last_last_year$YEAR<-NULL
  all_item_se_last_last_year$last_year<-NULL
  all_item_se_last_year_full<-merge(all_item_se_last_year, all_item_se_last_last_year, by=c("sp", "res"))
  year<-2020
  current_list<-list()
  for (year in unique(all_item_se_last_year$YEAR)){
    if (year==2100){
      next()
    }
    current<-all_item_se_last_year_full[YEAR==year]
    next_year<-all_item_se_last_year[YEAR==year+1]
    colnames(next_year)[c(4:11)]<-sprintf("next_year_%s", colnames(next_year)[c(4:11)])
    next_year$YEAR<-NULL
    next_year$last_year<-NULL
    current<-merge(current, next_year, by=c("sp", "res"))
    current_list[[length(current_list)+1]]<-current
  }
  current_list<-rbindlist(current_list)
  current_list$dist_next_year_max_y<-current_list$next_year_max_y - current_list$max_y
  current_list$dist_next_year_min_y<-current_list$next_year_min_y - current_list$min_y
  current_list$dist_next_year_max_abs_y<-current_list$next_year_max_abs_y - current_list$max_abs_y
  current_list$dist_next_year_min_abs_y<-current_list$next_year_min_abs_y - current_list$min_abs_y
  current_list$dist_next_year_quantile_5<-current_list$next_year_quantile_5 - current_list$quantile_5
  current_list$dist_next_year_quantile_95<-current_list$next_year_quantile_95 - current_list$quantile_95
  current_list$dist_next_year_quantile_abs_5<-current_list$next_year_quantile_abs_5 - current_list$quantile_abs_5
  current_list$dist_next_year_quantile_abs_95<-current_list$next_year_quantile_abs_95 - current_list$quantile_abs_95
  
  current_list$dist_last_year_max_y<-current_list$last_year_max_y - current_list$max_y
  current_list$dist_last_year_min_y<-current_list$last_year_min_y - current_list$min_y
  current_list$dist_last_year_max_abs_y<-current_list$last_year_max_abs_y - current_list$max_abs_y
  current_list$dist_last_year_min_abs_y<-current_list$last_year_min_abs_y - current_list$min_abs_y
  current_list$dist_last_year_quantile_5<-current_list$last_year_quantile_5 - current_list$quantile_5
  current_list$dist_last_year_quantile_95<-current_list$last_year_quantile_95 - current_list$quantile_95
  current_list$dist_last_year_quantile_abs_5<-current_list$last_year_quantile_abs_5 - current_list$quantile_abs_5
  current_list$dist_last_year_quantile_abs_95<-current_list$last_year_quantile_abs_95 - current_list$quantile_abs_95
  
  current_list$years<-current_list$last_year-current_list$YEAR
  current_list$GCM<-itemc$GCM
  current_list$SSP<-itemc$SSP
  current_list$exposure<-itemc$exposure
  current_list$group<-itemc$group
  full_list[[length(full_list)+1]]<-current_list
}
full_list_full<-rbindlist(full_list)
saveRDS(full_list_full, "../../Objects/Dispersal_distances/edge_based_dispersal_distance.rda")
