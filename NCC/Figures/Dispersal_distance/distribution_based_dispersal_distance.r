library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
full_set<-readRDS("../../Objects/Dispersal_distances/edge_based_dispersal_distance.rda")
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_quantile_abs_95)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_quantile_abs_5)/years)), by=list(group, SSP)]

full_set[YEAR==2020, .(mean_dist=smean(abs(dist_last_year_max_abs_y)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_max_y)/years)), by=list(group, SSP)]

full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_min_abs_y)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_min_y)/years)), by=list(group, SSP)]

full_set[, .(mean_dist=mean(abs(dist_next_year_max_abs_y)/years)), by=list(group, SSP)]

full_set[, .(mean_dist=mean(abs(dist_next_year_max_abs_y))), by=list(group, SSP)]

p<-ggplot(full_set[YEAR==2020])+geom_histogram(aes(x=abs(dist_last_year_max_abs_y/1000)/years), bins=20)+
  xlim(0, 50)+ylim(0, 10000)+
  facet_grid(group~SSP, scale="free")+theme_bw()+
  xlab("dispersal distance per year (km/yr)")
p
ggsave(p, filename="../../Figures/dispersal_distance/edge_based_dispersal_distance.png", width=6, height=4)

full_set_yearly_max<-full_set[, .(mean_dist=mean(abs(dist_next_year_max_abs_y/1000))),
                          by=list(group, SSP, YEAR)]
full_set_yearly_max$type<-"Max abs latitude"
full_set_yearly_min<-full_set[, .(mean_dist=mean(abs(dist_next_year_min_abs_y/1000))),
                              by=list(group, SSP, YEAR)]
full_set_yearly_min$type<-"Min abs latitude"
full_set_yearly<-rbindlist(list(full_set_yearly_max, full_set_yearly_min))

full_set[, .(mean_dist=mean(abs(dist_next_year_max_abs_y))),
         by=list(group, SSP, YEAR)]
full_set[, .(mean_dist=mean(abs(dist_next_year_quantile_abs_95))),
         by=list(group, SSP, YEAR)]
source("commonFuns/colors.r")
p<-ggplot(full_set_yearly)+geom_line(aes(x=YEAR, y=mean_dist, color=group))+
  facet_grid(SSP~type, scale="free")+theme_bw()+
  ylab("dispersal distance per year (km/yr)")+
  scale_color_manual(values=color_groups)
p
ggsave(p, filename="../../Figures/dispersal_distance/edge_based_dispersal_distance_details_abs.png", width=8, height=4)


full_set_yearly_max<-full_set[, .(mean_dist=mean(dist_next_year_max_y/1000)),
                              by=list(group, SSP, YEAR)]
full_set_yearly_max$type<-"Max abs latitude"
full_set_yearly_min<-full_set[, .(mean_dist=mean(dist_next_year_min_abs_y/1000)),
                              by=list(group, SSP, YEAR)]
full_set_yearly_min$type<-"Min abs latitude"
full_set_yearly<-rbindlist(list(full_set_yearly_max, full_set_yearly_min))

full_set[, .(mean_dist=mean(abs(dist_next_year_max_abs_y))),
         by=list(group, SSP, YEAR)]
full_set[, .(mean_dist=mean(abs(dist_next_year_quantile_abs_95))),
         by=list(group, SSP, YEAR)]
source("commonFuns/colors.r")
p<-ggplot(full_set_yearly)+geom_line(aes(x=YEAR, y=mean_dist, color=group))+
  facet_grid(SSP~type, scale="free")+theme_bw()+
  ylab("dispersal distance per year (km/yr)")+
  scale_color_manual(values=color_groups)
p
ggsave(p, filename="../../Figures/dispersal_distance/edge_based_dispersal_distance_details_raw.png", width=8, height=4)
