library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
full_set<-readRDS("../../Objects/Dispersal_distances/edge_based_dispersal_distance.rda")
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_quantile_abs_95)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_quantile_abs_5)/years)), by=list(group, SSP)]

full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_max_abs_y)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_max_y)/years)), by=list(group, SSP)]

full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_min_abs_y)/years)), by=list(group, SSP)]
full_set[YEAR==2020, .(mean_dist=mean(abs(dist_last_year_min_y)/years)), by=list(group, SSP)]


full_set[, .(mean_dist=mean(abs(dist_next_year_max_y)/years)), by=list(group, SSP)]

p<-ggplot(full_set[YEAR==2020])+geom_histogram(aes(x=abs(dist_last_year_max_abs_y/1000)/years), bins=20)+
  xlim(0, 50)+ylim(0, 10000)+
  facet_grid(group~SSP, scale="free")+theme_bw()+
  xlab("dispersal distance per year, km/year")
ggsave(p, filename="../../Figures/dispersal_distance/edge_based_dispersal_distance.png", width=6, height=4)
