library(raster)
library(ggplot2)
library(dplyr)
library(dplyr)
library(Rmisc)
g<-"Amphibians"
final_df<-readRDS(sprintf("../../Figures/dispersal_usage/%s.rda", g))
mask<-raster("../../Raster/mask_index.tif")
final_df_sum<-final_df%>%dplyr::group_by(mask_index, GCM, SSP)%>%
  dplyr::summarise(N_SP_SUM=sum(N_SP))

final_df_sum_se<-final_df_sum%>%dplyr::group_by(mask_index, SSP)%>%
  dplyr::summarise(MEAN_N_SP_SUM=mean(N_SP_SUM),
                   SD_N_SP_SUM=sd(N_SP_SUM),
                   CI_N_SP_SUM=CI(N_SP_SUM)[3]-CI(N_SP_SUM)[2])
mask_p<-data.frame(rasterToPoints(mask))
p<-left_join(mask_p, final_df_sum_se, by="mask_index")
r<-mask
values(r)[!is.na(values(r))]<-p$MEAN_N_SP_SUM
plot(r)
