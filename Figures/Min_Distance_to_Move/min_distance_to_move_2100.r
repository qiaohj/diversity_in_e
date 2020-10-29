library(ggplot2)
library(dplyr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
g<-"Amphibians"
source("functions.r")
if (F){
  df_all<-NULL
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    df<-readRDS(sprintf("../../Objects/Min_distance_to_Dispersal/%s.rda", g))
    #hist(df$dist_min)
    df_filter<-df%>%dplyr::filter((dist_min!=0)&!is.infinite(dist_min))
    df_filter<-df%>%dplyr::filter(!is.infinite(dist_min))
    df_se<-df_filter%>%dplyr::group_by(SSP, sp, group)%>%
      dplyr::summarise(dist_min_mean=mean(dist_min,),
                       dist_min_sd=sd(dist_min),
                       dist_min_CI=CI(dist_min)[2]-CI(dist_min)[3])
    df_all<-bind(df_all, df_se)
  }
  saveRDS(df_all, "../../Figures/Min_distance_to_Dispersal/2100.rda")
}
df_all<-readRDS("../../Figures/Min_distance_to_Dispersal/2100.rda")
p<-ggplot(df_all, 
       aes(x=dist_min_mean, y = ..density.., color=factor(group))) +
  geom_density()+scale_x_log10()+theme_bw()+
  facet_wrap(~SSP, ncol=1)
  

ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/2100.png")
