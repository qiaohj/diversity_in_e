library(raster)
library(data.table)
library(gdalUtilities)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
setDTthreads(1)
if (F){
  r<-raster("../../Raster/2020_LC_Type5.tif")
  r_sinu<-projectRaster(r, res=c(500, 500), crs=crs(mask_10km), method="ngb")
  writeRaster(r_sinu, filename="../../Raster/2020_LC_Type5_sinu.tif")
}
r_sinu<-raster("../../Raster/2020_LC_Type5_sinu_1km.tif")
print(sprintf("Current core number is %d", getDTthreads()))
esm_ssp<-c("EC-Earth3-Veg_SSP119", "MRI-ESM2-0_SSP119", "UKESM1_SSP119", 
           "EC-Earth3-Veg_SSP245", "MRI-ESM2-0_SSP245", "UKESM1_SSP245",
           "EC-Earth3-Veg_SSP585", "MRI-ESM2-0_SSP585", "UKESM1_SSP585")
group_df_birds<-readRDS("../../Data/Birds/bird_df.rda")
group_df_mammals<-readRDS("../../Data/Mammals/mammal_df.rda")
group_df_birds<-data.table(group="Birds", sp=unique(group_df_birds$SCINAME))
group_df_mammals<-data.table(group="Mammals", sp=unique(group_df_mammals$binomial))
group_df<-rbindlist(list(group_df_birds, group_df_mammals))

bi="Colius striatus"
coms<-expand.grid(exposure_threshold=c(0, 5), dispersal=c(0, 1))

i=1
j=1
group_df<-group_df[sample(nrow(group_df), nrow(group_df))]
threshold<-0.1
if (T){
  #group_df<-group_df[group=="Birds"]
  for (i in 1:length(group_df$sp)) {
    start_time<-Sys.time()
    bi<-group_df$sp[i]
    group<-group_df$group[i]
    print(paste(bi, group, i, length(group_df$sp), "Start time:", start_time))
    for (j in 1:nrow(coms)){
      dispersal<-coms[j, "dispersal"]
      exposure_threshold<-coms[j, "exposure_threshold"]
      
      
      target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
      fit_str<-sprintf("%s/fit.rda", target_folder)
      if (!file.exists(fit_str)){
        next()
      }
      item_str<-esm_ssp[1]
      for (item_str in esm_ssp){
        print(paste(i, nrow(group_df), "disp", dispersal, "expo", exposure_threshold, item_str))
        source<-sprintf("%s/%s_%d_dispersal_%d_lc.rda", target_folder, item_str,
                        exposure_threshold, dispersal)
        if (file.exists(source)){
          target<-sprintf("%s/%s_%d_dispersal_%d_lc_log.rda", target_folder, item_str,
                          exposure_threshold, dispersal)
          if (file.exists(target)){
            print("skip")
            #if (file.size(target)>100){
              next()
            #}
          }
          saveRDS(NULL, target)
          rr<-readRDS(source)
          if (length(rr)==2){
            if (raster::inMemory(rr[[1]])){
              r_2020<-rr[[1]]
            }else{
              r_2020<-raster(sprintf("%s/exposure_%s_dispersal_%s_lc_2020.tif", target_folder, exposure_threshold, dispersal))
            }
            p_2020<-data.table(rasterToPoints(r_2020))
            colnames(p_2020)[3]<-"X2020_LC_Type5_sinu_1km"
            p_2020_se<-p_2020[, .(N=.N), by=list(X2020_LC_Type5_sinu_1km)]
            p_2020_se$ALL<-nrow(p_2020)
            p_2020_se$PER<-p_2020_se$N/p_2020_se$ALL
            target_lc<-p_2020_se[PER>=threshold]$X2020_LC_Type5_sinu_1km
            p_2020_se$is_meaningfull<-p_2020_se$X2020_LC_Type5_sinu_1km %in% target_lc
            meaningful_area_2020<-sum(p_2020_se[X2020_LC_Type5_sinu_1km %in% target_lc]$N)
            if (raster::inMemory(rr[[2]])){
              r_2100<-rr[[2]]
            }else{
              if (file.exists(sprintf("%s/%s_exposure_%d_dispersal_%d_lc_2100.tif", target_folder, item_str, exposure_threshold, dispersal))){
                r_2100<-raster(sprintf("%s/%s_exposure_%d_dispersal_%d_lc_2100.tif", target_folder, item_str, exposure_threshold, dispersal))
              }else{
                r_2100<-rr[[2]]
              }
              
            }
            p_2100<-data.table(rasterToPoints(r_2100))
            colnames(p_2100)[3]<-"X2020_LC_Type5_sinu_1km"
            p_2100_se<-p_2100[, .(N=.N), by=list(X2020_LC_Type5_sinu_1km)]
            p_2100_se$ALL<-nrow(p_2100)
            p_2100_se$PER<-p_2100_se$N/p_2100_se$ALL
            p_2100_se$is_meaningfull<-p_2100_se$X2020_LC_Type5_sinu_1km %in% target_lc
            meaningful_area_2100<-sum(p_2100_se[X2020_LC_Type5_sinu_1km %in% target_lc]$N)
            item<-data.table(meaningful_area_2020=meaningful_area_2020, meaningful_area_2100=meaningful_area_2100,
                             change=meaningful_area_2100-meaningful_area_2020)
            item$sp<-bi
            item$group<-group
            item$dispersal<-dispersal
            item$exposure<-exposure_threshold
            esm_ssp_s<-strsplit(item_str, "_")[[1]]
            item$esm<-esm_ssp_s[1]
            item$ssp<-esm_ssp_s[2]
            item$change_per<-item$change/item$meaningful_area_2020
            log<-list(log=item, p_2020=p_2020_se, p_2100=p_2100_se)
            saveRDS(log, target)
          }
        }
      }
    }
  }
}

asdf

if (F){
  
  all_df<-list()
  for (i in 1:length(group_df$sp)) {
    start_time<-Sys.time()
    bi<-group_df$sp[i]
    group<-group_df$group[i]
    print(paste(bi, group, i, length(group_df$sp), "Start time:", start_time))
    for (j in 1:nrow(coms)){
      dispersal<-coms[j, "dispersal"]
      exposure_threshold<-coms[j, "exposure_threshold"]
      
      
      target_folder<-sprintf("../../Objects/Dispersal/%s/%s", group, gsub(" ", "_", bi))
      fit_str<-sprintf("%s/fit.rda", target_folder)
      if (!file.exists(fit_str)){
        next()
      }
      item_str<-esm_ssp[1]
      for (item_str in esm_ssp){
        print(paste(i, nrow(group_df), "disp", dispersal, "expo", exposure_threshold, item_str))
        source<-sprintf("%s/%s_%d_dispersal_%d_lc.rda", target_folder, item_str,
                        exposure_threshold, dispersal)
        if (file.exists(source)){
          target<-sprintf("%s/%s_%d_dispersal_%d_lc_log.rda", target_folder, item_str,
                          exposure_threshold, dispersal)
          if (file.exists(target)){
            if (file.size(target)>100){
              item<-readRDS(target)
              all_df[[length(all_df)+1]]<-item$log
            }
          }
        }
      }
    }
  }
  
  all_df<-rbindlist(all_df)
  saveRDS(all_df, "../../Objects/landcover_change.rda")
}
all_df<-readRDS("../../Objects/landcover_change.rda")
all_df<-all_df[!is.na(change_per)]
quantiles_iqr<-quantile(all_df$change_per, c(0.25, 0.75))
iqr<-quantiles_iqr[2]-quantiles_iqr[1]
iqr_range<-c(quantiles_iqr[1]-1.5*iqr, quantiles_iqr[2]+1.5*iqr)
all_df_iqr<-all_df[between(change_per, iqr_range[1], iqr_range[2])]
hist(all_df_iqr$change_per)

all_df_iqr_se<-all_df_iqr[, .(change_per=mean(change_per), sd_change_per=sd(change_per)),
                          by=list(group, dispersal, exposure, ssp)]

ggplot(all_df_iqr)+geom_density(aes(x=change_per, fill=ssp, color=ssp), alpha=0.2)+facet_grid(dispersal~exposure)+theme_bw()

all_df_iqr$loss_99<-all_df_iqr$change_per<=-0.99
all_df_iqr$loss_95<-all_df_iqr$change_per<=-0.95
all_df_iqr$loss_90<-all_df_iqr$change_per<=-0.90
all_df_iqr_sp<-all_df_iqr[, .(N_sp=length(unique(sp))), by=list(dispersal, exposure, esm, ssp)]
all_df_iqr_99<-all_df_iqr[, .(N=.N/3), by=list(dispersal, exposure, esm, ssp, loss_99)]
all_df_iqr_99<-merge(all_df_iqr_99, all_df_iqr_sp, by=c("dispersal", "exposure", "esm", "ssp"))
all_df_iqr_99$per<-all_df_iqr_99$N/all_df_iqr_99$N_sp
all_df_iqr_99<-all_df_iqr_99[loss_99==T]
all_df_iqr_99_se<-all_df_iqr_99[, .(per=mean(per)), by=list(dispersal, exposure, ssp)]
all_df_iqr_99_se$threshold<-"99%"

all_df_iqr_95<-all_df_iqr[, .(N=.N/3), by=list(dispersal, exposure, esm, ssp, loss_95)]
all_df_iqr_95<-merge(all_df_iqr_95, all_df_iqr_sp, by=c("dispersal", "exposure", "esm", "ssp"))
all_df_iqr_95$per<-all_df_iqr_95$N/all_df_iqr_95$N_sp
all_df_iqr_95<-all_df_iqr_95[loss_95==T]
all_df_iqr_95_se<-all_df_iqr_95[, .(per=mean(per)), by=list(dispersal, exposure, ssp)]
all_df_iqr_95_se$threshold<-"95%"

all_df_iqr_90<-all_df_iqr[, .(N=.N/3), by=list(dispersal, exposure, esm, ssp, loss_90)]
all_df_iqr_90<-merge(all_df_iqr_90, all_df_iqr_sp, by=c("dispersal", "exposure", "esm", "ssp"))
all_df_iqr_90$per<-all_df_iqr_90$N/all_df_iqr_90$N_sp
all_df_iqr_90<-all_df_iqr_90[loss_90==T]
all_df_iqr_90_se<-all_df_iqr_90[, .(per=mean(per)), by=list(dispersal, exposure, ssp)]
all_df_iqr_90_se$threshold<-"90%"

all_df_iqr_all<-rbindlist(list(all_df_iqr_90_se, all_df_iqr_95_se, all_df_iqr_99_se))

p<-ggplot()+
  geom_bar(data=all_df_iqr_all, 
           stat="identity", position="dodge2", 
           aes(y=per*100, x=ssp, fill=threshold, group=group), width=2, color="grey")+
  xlab("SSP scenario")+
  #ggtitle(sprintf("Distribution>%d", ttt))+
  #ylim(c(0, 1))+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  facet_grid(exposure~dispersal)+
  labs(fill = "Threshold")+
  ylab("Endangrous because of loss of landcover")
p
