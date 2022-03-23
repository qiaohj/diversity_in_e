library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  bird_property<-data.table(read.csv("../../Objects/Supp.Tables/bird_property_10km.csv", stringsAsFactors = F))
  bird_property$sp<-gsub(" ", "_", bird_property$sp)
  group<-"Birds"
  exposure<-0
  result<-list()
  for (GCM_i in GCMs){
    for (SSP_i in SSPs){
      for (group in c("Birds", "Mammals")){
        for (exposure in c(5)){
          f<-sprintf("../../Objects/cluster_based_pathway_10km/merged/%s_%s_%s_exposure_%d_sub_100.rda",
                     group, GCM_i, SSP_i, exposure)
          print(f)
          df_bird<-readRDS(f)
          
          df_bird_next<-df_bird[c(2:nrow(df_bird))]
          df_bird$next_x<-0
          df_bird$next_y<-0
          df_bird$next_group<-""
          
          df_bird[c(1:(nrow(df_bird))-1)]$next_x<-df_bird_next$x
          df_bird[c(1:(nrow(df_bird))-1)]$next_y<-df_bird_next$y
          df_bird[c(1:(nrow(df_bird))-1)]$next_group<-df_bird_next$line_group
          df_bird$dist<-sqrt((df_bird$x-df_bird$next_x)^2+(df_bird$y-df_bird$next_y)^2)
          dim(df_bird)
          df_bird_filter<-df_bird[line_group==next_group]
          df_bird_filter<-df_bird_filter[next_x!=0]
          df_bird_filter$YEAR_group<-floor(df_bird_filter$YEAR)
          quantiles<-quantile(df_bird_filter$dist, c(0.25, 0.75))
          iqr<-quantiles[2]-quantiles[1]
          upperlimit<-quantiles[2]+1.5*iqr
          df_bird_filter<-df_bird_filter[dist<=upperlimit]
          df_bird_filter_se<-df_bird_filter[, .(sum_dist=sum(dist)), by=list(sp, YEAR_group)]
          
          df_bird_filter_se<-df_bird_filter_se[,.(mean_dist=mean(sum_dist)), by=list(sp)]
          item<-data.frame(is_migratory_bird="ALL", 
                           mean_dist=mean(df_bird_filter_se$mean_dist)/1000,
                           sd_dist=sd(df_bird_filter_se$mean_dist/1000))
          if (group=="Birds"){
            
            df_bird_filter_se<-merge(df_bird_filter_se, bird_property, by="sp")
            
            item_mig<-df_bird_filter_se[, .(mean_dist=mean(mean_dist/1000), sd_dist=sd(mean_dist/1000)),
                                        by=(is_migratory_bird)]
            item_mig$is_migratory_bird<-as.character(item_mig$is_migratory_bird)
            item<-rbindlist(list(item,item_mig))
          }
          item$group<-group
          item$GCM<-GCM_i
          item$SSP<-SSP_i
          item$exposure<-exposure
          result[[length(result)+1]]<-item
        }
      }
    }
  }
  result<-rbindlist(result)
  saveRDS(result, "../../Objects/Dispersal_distances/cluster_based_disp_distance.rda")
  
}
result<-readRDS("../../Objects/Dispersal_distances/cluster_based_disp_distance.rda")
result<-result[, .(mean_dist=mean(mean_dist), sd_dist=mean(sd_dist)),
               by=list(is_migratory_bird, group)]
