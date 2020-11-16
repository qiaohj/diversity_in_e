library(dplyr)
library(ggplot2)
library(raster)
library(Rmisc)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
group<-"Mammals"
if (F){
  
  setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
  group<-"Birds"
  df<-readRDS(sprintf("../../Objects/Species_exposure/%s_accumulative.rda", group))
  #1: in range, but out mve
  #2: out of range because of the max temp
  #3: out of range because of the prec
  #4: out of range because of the prec and temp
  
  df$exposure_ratio<-df$N_exposure/df$N_SP
  df$range_type_0_ratio<-df$range_type_0/df$N_SP
  df$range_type_1_ratio<-df$range_type_1/df$N_SP
  df$range_type_2_ratio<-df$range_type_2/df$N_SP
  df$range_type_3_ratio<-df$range_type_3/df$N_SP
  df$range_type_4_ratio<-df$range_type_4/df$N_SP
  
  
  exposure_se<-df%>%dplyr::group_by(mask_index, SSP, year)%>%
    dplyr::summarise(exposure_ratio_mean=mean(exposure_ratio),
                     exposure_ratio_sd=sd(exposure_ratio),
                     range_type_0_ratio_mean=mean(range_type_0_ratio),
                     range_type_0_ratio_sd=sd(range_type_0_ratio),
                     range_type_1_ratio_mean=mean(range_type_1_ratio),
                     range_type_1_ratio_sd=sd(range_type_1_ratio),
                     range_type_2_ratio_mean=mean(range_type_2_ratio),
                     range_type_2_ratio_sd=sd(range_type_2_ratio),
                     range_type_3_ratio_mean=mean(range_type_3_ratio),
                     range_type_3_ratio_sd=sd(range_type_3_ratio),
                     range_type_4_ratio_mean=mean(range_type_4_ratio),
                     range_type_4_ratio_sd=sd(range_type_4_ratio),
                     N_SP_mean=mean(N_SP),
                     N_SP_sd=sd(N_SP),
                     N_exposure_mean=mean(N_exposure),
                     N_exposure_sd=sd(N_exposure),
                     range_type_0_mean=mean(range_type_0),
                     range_type_0_ratio_sd=sd(range_type_0),
                     range_type_1_mean=mean(range_type_1),
                     range_type_1_ratio_sd=sd(range_type_1),
                     range_type_2_mean=mean(range_type_2),
                     range_type_2_ratio_sd=sd(range_type_2),
                     range_type_3_mean=mean(range_type_3),
                     range_type_3_ratio_sd=sd(range_type_3),
                     range_type_4_mean=mean(range_type_4),
                     range_type_4_ratio_sd=sd(range_type_4)
                     )
  
  mask<-raster("../../Raster/mask_index.tif")
  mask_p<-data.frame(rasterToPoints(mask))
  exposure_se<-inner_join(exposure_se, mask_p, by=c("mask_index"))
  
  saveRDS(exposure_se, sprintf("../../Objects/Species_exposure/%s_exposure_se_acc.rda", group))
  
  
  head(exposure_se[which(exposure_se$exposure_ratio_mean>0.334),])
  hist(exposure_se$exposure_ratio_mean)
  
}
df<-readRDS(sprintf("../../Objects/Species_exposure/%s_exposure_se_acc.rda", group))

sp_richness<-readRDS(sprintf("../../Objects/Diversity/%s/EC-Earth3-Veg_SSP119_0_1/indices_df.rda", group))
sp_richness_st<-sp_richness[["2014"]]
head(sp_richness_st$species.richness)
species.richness<-sp_richness_st$species.richness

df_with_richness<-left_join(df, species.richness, by=c("mask_index"="index"))

mask<-raster("../../Raster/mask_index.tif")
p_mask<-data.frame(rasterToPoints(mask))
no_na<-!is.na(values(mask))
yyy<-2015
SSP_item<-"SSP585"
dir.create(sprintf("../../Objects/Species_exposure/Exposure_by_year/%s", group))
dir.create(sprintf("../../Objects/Species_exposure/Exposure_by_year/%s/TIF", group))
dir.create(sprintf("../../Objects/Species_exposure/Exposure_by_year/%s/PNG", group))
pal <- colorRampPalette(c("#0C7BDC","#E1BE6A", "#DC3220"))
for (yyy in c(2015:2100)){
  print(yyy)
  png(sprintf("../../Objects/Species_exposure/Exposure_by_year/%s/PNG/%d.png", group, yyy), width=2000, height=700)
  par(mfrow = c(1, 3))
  for (SSP_item in c("SSP119", "SSP245", "SSP585")){
    df_with_richness_y<-df_with_richness%>%dplyr::filter((year==yyy)&(SSP==SSP_item))
    df_with_richness_y<-left_join(p_mask, df_with_richness_y, by="mask_index")
    #head(df_with_richness_y[which((!is.na(df_with_richness_y$metric))&(is.na(df_with_richness_y$exposure_ratio_mean))),])
    if (yyy==2015){
      r_richness<-mask
      values(r_richness)[no_na]<-df_with_richness_y$metric
      writeRaster(r_richness, sprintf("../../Objects/Species_Richness/%s.tif", group), overwrite=T)
    }
    r_richness<-mask
    values(r_richness)[no_na]<-df_with_richness_y$exposure_ratio_mean
    writeRaster(r_richness, 
                sprintf("../../Objects/Species_exposure/Exposure_by_year/%s/TIF/%s_%d.tif", group, SSP_item, yyy),
                overwrite=T)
    
    
    plot(r_richness, breaks=seq(0, 1, 0.01), col=pal(100), main=sprintf("%s %s %d", group, SSP_item, yyy))
    
  }
  dev.off()
}

ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../Amphibians_Species_Exposure.mp4
ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../Birds_Species_Exposure.mp4
ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../Reptiles_Species_Exposure.mp4
ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../../Mammals_Species_Exposure.mp4

