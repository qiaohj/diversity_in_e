library(ggplot2)
library(dplyr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
g<-"Amphibians"
source("functions.r")
if (F){
  df_all<-NULL
  for (y in c(2015:2100)){
    print(y)
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      df<-readRDS(sprintf("../../Objects/Min_distance_to_Dispersal/%s/%d.rda", g, y))
      #hist(df$dist_min)
      #df_filter<-df%>%dplyr::filter((dist_min!=0)&!is.infinite(dist_min))
      df_filter<-df%>%dplyr::filter(!is.infinite(dist_min))
      df_se<-df_filter%>%dplyr::group_by(SSP, sp, group)%>%
        dplyr::summarise(dist_min_mean=mean(dist_min,),
                         dist_min_sd=sd(dist_min),
                         dist_min_CI=CI(dist_min)[2]-CI(dist_min)[3])
      df_se$year<-y
      df_all<-bind(df_all, df_se)
    }
  }
  saveRDS(df_all, "../../Figures/Min_distance_to_Dispersal/full.rda")
}
df_all<-readRDS("../../Figures/Min_distance_to_Dispersal/full.rda")
p<-ggplot(df_all%>%dplyr::filter(year %in% c(2015, 2100)), 
       aes(x=dist_min_mean, y = ..density.., color=factor(group))) +
  geom_density(aes(linetype=factor(year)))+scale_x_log10()+theme_bw()+
  facet_wrap(~SSP, ncol=1)
  

ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/start_end.png")

df_all_se<-df_all%>%dplyr::group_by(SSP, group, year)%>%
  dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                   sd_dist_mean=sd(dist_min_mean))

p<-ggplot(df_all_se, aes(x=year, y = mean_dist_mean, color=factor(group)))+
  geom_line()+theme_bw()+
  facet_wrap(~SSP, ncol=1)
ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/mean_by_year.png")

#p<-ggplot(df_all, aes(x=year, y = dist_min_mean, color=factor(group)))+
#  geom_point()+theme_bw()+
#  facet_wrap(~SSP, ncol=1)

df_all_se<-df_all%>%dplyr::filter(dist_min_mean>=1)%>%dplyr::group_by(SSP, group, year)%>%
  dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                   sd_dist_mean=sd(dist_min_mean))

p<-ggplot(df_all_se, aes(x=year, y = mean_dist_mean, color=factor(group)))+
  geom_line()+theme_bw()+
  facet_wrap(~SSP, ncol=1)
ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/mean_by_year_without_zero.png")

