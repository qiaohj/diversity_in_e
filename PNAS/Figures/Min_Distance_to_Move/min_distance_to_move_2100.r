library(ggplot2)
library(dplyr)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
g<-"Birds"
source("commonFuns/functions.r")

if (F){
  exposure=0
  df_extincts<-list()
  for (g in c("Birds", "Mammals")){
    df_extinct<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
    df_extinct<-df_extinct%>%dplyr::filter(dispersal==1)
    df_extinct<-df_extinct%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
    df_extinct[is.infinite(df_extinct$extinct_year), "extinct_year"]<-2020
    df_extincts[[g]]<-df_extinct
  }
  df_all<-NULL
  for (y in c(2021:2100)){
    print(y)
    for (g in c("Birds", "Mammals")){
      
      df<-readRDS(sprintf("../../Objects/Min_distance_to_Dispersal/%s/%d.rda", g, y))
      df_extinct<-df_extincts[[g]]
      df$sp<-df$SP
      df$group<-g
      df<-left_join(df, df_extinct, by=c("group", "sp", "GCM", "SSP"))
      df[is.na(df$extinct_year), "extinct_year"]<-9999
      #hist(df$dist_min)
      #df_filter<-df%>%dplyr::filter((dist_min!=0)&!is.infinite(dist_min))
      df_filter<-df%>%dplyr::filter(!is.infinite(dist_min))
      #df_se<-df_filter%>%dplyr::group_by(SSP, sp, group, extinct_year)%>%
      #  dplyr::summarise(dist_min_mean=mean(dist_min),
      #                   dist_min_sd=sd(dist_min),
      #                   dist_min_CI=CI(dist_min)[2]-CI(dist_min)[3])
      #df_se$year<-y
      df_all<-bind(df_all, df_filter)
    }
  }
  saveRDS(df_all, "../../Figures/Min_distance_to_Dispersal/full.rda")
}
source("commonFuns/colors.r")
df_all<-readRDS("../../Figures/Min_distance_to_Dispersal/full.rda")

df_sp_list<-list()
exposure<-0
for (group in c("Birds", "Mammals")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list, fill=T)

df_sp_list$sp<-gsub(" ", "_", df_sp_list$SP)
df_all<-df_all[sp %in% df_sp_list$sp]

df_all_SSP<-df_all%>%dplyr::group_by(SSP, group, year, extinct_year)%>%
  dplyr::summarise(dist_min_mean=mean(dist_min),
                   dist_min_sd=sd(dist_min),
                   dist_min_CI=CI(dist_min)[2]-CI(dist_min)[3])
                   
p<-ggplot(df_all_SSP%>%dplyr::filter((year %in% c(2100))&(dist_min_mean>-1)), 
       aes(x=dist_min_mean, fill=group)) +
  geom_histogram(bins=20, position="dodge")+theme_bw()+
  scale_fill_manual(values=color_groups)+
  xlab("Minimum distance need to dispersal")+
  ylab("Number of species")+
  facet_wrap(~SSP, ncol=1)
  
p
ggsave(p, filename=sprintf("../../Figures/Min_distance_to_Dispersal/start_end.png"))
yyy=2040

df_all_se<-df_all%>%dplyr::group_by(SSP, group, year)%>%
  dplyr::summarise(mean_dist_mean=mean(dist_min),
                   sd_dist_mean=sd(dist_min),
                   ci_dist_mean=CI(dist_min)[1]-CI(dist_min)[2])
write.csv(df_all_se, sprintf("../../Figures/Min_distance_to_Dispersal/mean.csv"))


for (yyy in c(2040, 2100)){
  df_all_se<-df_all%>%dplyr::filter(extinct_year>yyy)%>%dplyr::group_by(SSP, group, year)%>%
    dplyr::summarise(mean_dist_mean=mean(dist_min),
                     sd_dist_mean=sd(dist_min),
                     ci_dist_mean=CI(dist_min)[1]-CI(dist_min)[2])
  df_all_se%>%filter(year==yyy)
  write.csv(df_all_se, sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_%d.csv", yyy))
  
  df_all_se_grouped<-df_all%>%dplyr::filter(extinct_year>yyy)%>%group_by(group, year)%>%
    dplyr::summarise(mean_dist_mean=mean(dist_min),
                     sd_dist_mean=sd(dist_min),
                     ci_dist_mean=CI(dist_min)[1]-CI(dist_min)[2])
  df_all_se_grouped%>%filter(year==yyy)
  
  
  write.csv(df_all_se_grouped, 
            sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_all_SSP_%d.csv", yyy))
  
  p<-ggplot(df_all_se)+
    geom_ribbon(data=df_all_se, 
                aes(x=year,
                    ymin=mean_dist_mean-ci_dist_mean, 
                    ymax=mean_dist_mean+ci_dist_mean, 
                    fill=group), alpha=0.2)+
    geom_line(data=df_all_se, aes(x=year, y = mean_dist_mean, color=group))+theme_bw()+
    xlab("Year")+
    ylab("Average distance need to dispersal (KM)")+
    #ggtitle(yyy)+
    labs(color="Group", fill="Group")+
    scale_color_manual(values=color_groups)+
    scale_fill_manual(values=color_groups)+
    facet_wrap(~SSP, ncol=1)
  p
  ggsave(p, filename=sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_%d.png", yyy))
  ggsave(p, filename=sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_%d.png", yyy))
  
  p<-ggplot(df_all_se)+
    geom_ribbon(data=df_all_se, 
                aes(x=year,
                    ymin=mean_dist_mean-ci_dist_mean, 
                    ymax=mean_dist_mean+ci_dist_mean, 
                    fill=group), alpha=0.2)+
    geom_line(data=df_all_se, aes(x=year, y = mean_dist_mean, color=group))+
    theme_bw()+
    xlab("Year")+
    #ggtitle(yyy)+
    ylab("Average distance need to dispersal (KM)")+
    labs(color="Group", fill="Group")+
    scale_color_manual(values=color_groups)+
    scale_fill_manual(values=color_groups)+
    xlim(2020, 2100)+
    facet_wrap(~SSP, ncol=1, scale="free")
  p
  ggsave(p, filename=sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_scale_free_%d.png", yyy))
  ggsave(p, filename=sprintf("../../Figures/Min_distance_to_Dispersal/mean_by_year_scale_free_%d.pdf", yyy))
}

#p<-ggplot(df_all, aes(x=year, y = dist_min_mean, color=factor(group)))+
#  geom_point()+theme_bw()+
#  facet_wrap(~SSP, ncol=1)

if (F){
  df_all_se<-df_all%>%dplyr::filter(dist_min_mean>=1)%>%dplyr::group_by(SSP, group, year)%>%
    dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                     CI_dist_mean=sd(dist_min_mean))
  
  p<-ggplot(df_all_se, aes(x=year, y = mean_dist_mean, color=group))+
    geom_line()+theme_bw()+
    xlab("Year")+
    ylab("Average distance need to dispersal")+
    scale_color_manual(values=color_groups)+
    facet_wrap(~SSP, ncol=1)
  ggsave(p, filename="../../Figures/Min_distance_to_Dispersal/mean_by_year_without_zero.png")
  
  df_all_se_all<-df_all%>%dplyr::filter(dist_min_mean>=1)%>%dplyr::group_by(group, year)%>%
    dplyr::summarise(mean_dist_mean=mean(dist_min_mean),
                     CI_dist_mean=sd(dist_min_mean))
  
  
  df_all_se%>%filter(year==2100)
  
  
}