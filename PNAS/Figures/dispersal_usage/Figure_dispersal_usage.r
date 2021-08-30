library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)

library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/colors.r")
g<-"Reptiles"

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")

if (F){
  threshold<-5
  p_list<-list()
  final_df_sum_all<-NULL
  for (g in groups){
    print(g)
    final_df<-readRDS(sprintf("../../Figures_Full_species/dispersal_usage_%d/%s.rda", threshold, g))
    
    final_df_sum<-final_df%>%dplyr::group_by(mask_index, GCM, SSP, dispersal)%>%
      dplyr::summarise(N_SP_SUM=sum(N_SP))
    final_df_sum$group<-g
    
    final_df_sum_se<-final_df_sum%>%dplyr::group_by(mask_index, SSP, dispersal)%>%
      dplyr::summarise(MEAN_N_SP_SUM=mean(N_SP_SUM),
                       SD_N_SP_SUM=sd(N_SP_SUM, na.rm=T),
                       CI_N_SP_SUM=CI(N_SP_SUM)[3]-CI(N_SP_SUM)[2])
    final_df_sum_se<-inner_join(final_df_sum_se, mask_p, by="mask_index")
    #final_df_sum_se<-final_df_sum_se[which(final_df_sum_se$MEAN_N_SP_SUM>=quantile(final_df_sum_se$MEAN_N_SP_SUM, 0.5)),]
    final_df_sum_se$MEAN_N_SP_SUM_log<-log2(final_df_sum_se$MEAN_N_SP_SUM)
    
    
    if (is.null(final_df_sum_all)){
      final_df_sum_all<-final_df_sum
    }else{
      final_df_sum_all<-bind_rows(final_df_sum_all, final_df_sum)
    }
    p<-ggplot()+geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(data=final_df_sum_se%>%dplyr::filter(dispersal==1), aes(x=x, y=y, fill=MEAN_N_SP_SUM_log))+
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
  ggsave(p, filename=sprintf("../../Figures/dispersal_usage_%d/dispersal_usage_1.png", threshold), width=10, height=10)
  ggsave(p, filename=sprintf("../../Figures/dispersal_usage_%d/dispersal_usage_1.pdf", threshold), width=10, height=10)
  
  
  
  final_df_sum_all_se<-final_df_sum_all%>%dplyr::group_by(mask_index, SSP, dispersal)%>%
    dplyr::summarise(MEAN_N_SP_SUM=mean(N_SP_SUM),
                     SD_N_SP_SUM=sd(N_SP_SUM, na.rm=T),
                     CI_N_SP_SUM=CI(N_SP_SUM)[3]-CI(N_SP_SUM)[2])
  final_df_sum_all_se<-inner_join(final_df_sum_all_se, mask_p, by="mask_index")
  #final_df_sum_se<-final_df_sum_se[which(final_df_sum_se$MEAN_N_SP_SUM>=quantile(final_df_sum_se$MEAN_N_SP_SUM, 0.5)),]
  final_df_sum_all_se$MEAN_N_SP_SUM_log<-log2(final_df_sum_all_se$MEAN_N_SP_SUM)
  saveRDS(final_df_sum_all_se, file=sprintf("../../Figures/dispersal_usage_all/data_%d.rda", threshold))
  
}

df1<-readRDS(sprintf("../../Figures/dispersal_usage_all/data_%d.rda", 1))
df1$exposure<-" no climate reslience"

df5<-readRDS(sprintf("../../Figures/dispersal_usage_all/data_%d.rda", 5))
df5$exposure<-"climate resilience"

df<-bind_rows(df1, df5)

p<-ggplot()+geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  geom_tile(data=df%>%dplyr::filter(dispersal==1), aes(x=x, y=y, fill=MEAN_N_SP_SUM_log))+
  scale_fill_gradient(low=color_two_map[2], high=color_two_map[1])+
  map_theme+
  facet_grid(exposure~SSP)
p
ggsave(p, filename="../../Figures/dispersal_usage_all/dispersal_usage.png", width=12, height=6)
ggsave(p, filename="../../Figures/dispersal_usage_all/dispersal_usage.pdf", width=12, height=6)
