library(raster)
library(ggplot2)
library(dplyr)

library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
g<-"Reptiles"

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")

p_list<-list()
for (g in groups){
  print(g)
  final_df<-readRDS(sprintf("../../Figures/dispersal_usage_5/%s.rda", g))
  
  final_df_sum<-final_df%>%dplyr::group_by(mask_index, GCM, SSP)%>%
    dplyr::summarise(N_SP_SUM=sum(N_SP))
  
  final_df_sum_se<-final_df_sum%>%dplyr::group_by(mask_index, SSP)%>%
    dplyr::summarise(MEAN_N_SP_SUM=mean(N_SP_SUM),
                     SD_N_SP_SUM=sd(N_SP_SUM, na.rm=T),
                     CI_N_SP_SUM=CI(N_SP_SUM)[3]-CI(N_SP_SUM)[2])
  final_df_sum_se<-inner_join(final_df_sum_se, mask_p, by="mask_index")
  #final_df_sum_se<-final_df_sum_se[which(final_df_sum_se$MEAN_N_SP_SUM>=quantile(final_df_sum_se$MEAN_N_SP_SUM, 0.5)),]
  final_df_sum_se$MEAN_N_SP_SUM_log<-log2(final_df_sum_se$MEAN_N_SP_SUM)
  
  
  p<-ggplot()+geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
    geom_tile(data=final_df_sum_se, aes(x=x, y=y, fill=MEAN_N_SP_SUM_log))+
    scale_fill_gradient(low=color_two_map[2], high=color_two_map[1])+
    map_theme+facet_wrap(~SSP)+
    ggtitle(g)
  p_list[[g]]<-p
}
p<-ggarrange(p_list[[groups[1]]],
          p_list[[groups[2]]],
          p_list[[groups[3]]],
          p_list[[groups[4]]],
          ncol=1)
p
ggsave(p, filename="../../Figures/dispersal_usage_5/dispersal_usage.png", width=10, height=12)
ggsave(p, filename="../../Figures/dispersal_usage_5/dispersal_usage.pdf", width=10, height=12)
