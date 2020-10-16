library(ggplot2)
library(raster)
library(dplyr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
group<-"Mammals"
df<-readRDS(sprintf("../../Objects/Species_exposure/%s.rda", group))
df$exposure_ratio<-df$N_exposure/df$N_SP

exposure_se<-df%>%dplyr::group_by(mask_index, SSP, year)%>%
  dplyr::summarise(exposure_ratio_mean=mean(exposure_ratio),
                   exposure_ratio_sd=sd(exposure_ratio))

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
exposure_se<-inner_join(exposure_se, mask_p, by=c("mask_index"))

saveRDS(exposure_se, sprintf("../../Objects/Species_exposure/%s_exposure_se.rda", group))


head(exposure_se[which(exposure_se$exposure_ratio>0.334),])
hist(df$exposure_ratio)
