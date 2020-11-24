library(ggplot2)
library(dplyr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
g<-"Amphibians"
source("commonFuns/functions.r")
if (F){
  df_all<-NULL
  for (y in c(2021:2100)){
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
source("commonFuns/colors.r")
df_all<-readRDS("../../Figures/Min_distance_to_Dispersal/full.rda")
p<-ggplot(df_all%>%dplyr::filter((year %in% c(2100))&(dist_min_mean>-1)), 
       aes(x=dist_min_mean, fill=group)) +
  geom_histogram(bins=20, position="dodge")+theme_bw()+
  scale_fill_manual(values=color_groups)+
  xlab("Minimum distance need to dispersal")+
  ylab("Number of species")+
  facet_wrap(~SSP, ncol=1)
  
p
ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/start_end.png")

df_all_se<-df_all%>%dplyr::group_by(SSP, group, year)%>%
  dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                   sd_dist_mean=sd(dist_min_mean))

p<-ggplot(df_all_se, aes(x=year, y = mean_dist_mean, color=group))+
  geom_line()+theme_bw()+
  xlab("Year")+
  ylab("Average distance need to dispersal")+
  scale_color_manual(values=color_groups)+
  facet_wrap(~SSP, ncol=1)
p
ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/mean_by_year.png")

#p<-ggplot(df_all, aes(x=year, y = dist_min_mean, color=factor(group)))+
#  geom_point()+theme_bw()+
#  facet_wrap(~SSP, ncol=1)

df_all_se<-df_all%>%dplyr::filter(dist_min_mean>=1)%>%dplyr::group_by(SSP, group, year)%>%
  dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                   sd_dist_mean=sd(dist_min_mean))

p<-ggplot(df_all_se, aes(x=year, y = mean_dist_mean, color=group))+
  geom_line()+theme_bw()+
  xlab("Year")+
  ylab("Average distance need to dispersal")+
  scale_color_manual(values=color_groups)+
  facet_wrap(~SSP, ncol=1)
ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/mean_by_year_without_zero.png")

