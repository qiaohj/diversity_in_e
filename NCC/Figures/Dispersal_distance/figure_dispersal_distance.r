library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
all_disp<-readRDS("../../Objects/Dispersal_distances/all_mean_disp_dist.rda")
group_disp_bird<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
group_disp_mammal<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
cols_bird<-c("iucn_name", "Diet", "estimated_disp", "Migration_3")
group_disp_bird<-group_disp_bird[, ..cols_bird]
colnames(group_disp_bird)<-c("sp", "Diet", "estimated_disp", "Migration")
group_disp_bird$group<-"Birds"
group_disp_bird[is.na(group_disp_bird$Migration)]$Migration<-"unknown"
cols_mammals<-c("Scientific", "Diet", "estimated_disp")
group_disp_mammal<-group_disp_mammal[, ..cols_mammals]
colnames(group_disp_mammal)<-c("sp", "Diet", "estimated_disp")
group_disp_mammal$Migration<-"unknown"
group_disp_mammal$group<-"Mammals"

group_disp<-rbindlist(list(group_disp_mammal, group_disp_bird))

all_disp_all<-merge(all_disp, group_disp, by=c("sp", "group"))
all_disp_all_mammals<-all_disp_all[group=="Mammals"]
all_disp_all_birds<-all_disp_all[group=="Birds"]

quantiles<-quantile(all_disp_all_mammals$mean_disp, c(0.25, 0.75))
iqr<-quantiles[2]-quantiles[1]
iqr_range<-c(quantiles[1]-1.5*iqr, quantiles[2]+1.5*iqr)
quantiles<-quantile(all_disp_all_birds$mean_disp, c(0.25, 0.75))
iqr<-quantiles[2]-quantiles[1]
iqr_range<-c(quantiles[1]-1.5*iqr, quantiles[2]+1.5*iqr)
all_disp_all_mammals<-all_disp_all_mammals[between(mean_disp, iqr_range[1], iqr_range[2])]


all_disp_all_birds<-all_disp_all_birds[between(mean_disp, iqr_range[1], iqr_range[2])]

all_disp_all2<-rbindlist(list(all_disp_all_birds, all_disp_all_mammals))

colorBlindBlack8  <- c("plants"="#E69F00", 
                       "invertebrates"="#56B4E9",
                       "omnivore"="#009E73", 
                       "vertebrates"="#F0E442", 
                       "herbivore"="#0072B2", 
                       "carnivore"="#D55E00")

p<-ggplot()+geom_histogram(data=all_disp_all2, 
                      aes(x=mean_disp/1000, fill=Diet), 
                      alpha=0.5, bins = 50)+
  facet_wrap(~group, scale="free")+ylab("Number of species")+xlab("mean dispersal distance")+
  scale_fill_manual(values=colorBlindBlack8)+
  theme_bw()
p
ggsave(p, filename="../../Figures/dispersal_distance/dispersal_distance_diet.png", width=12, height=6)

unique(all_disp_all2$Migration)
all_disp_all2$Migration_str<-all_disp_all2$Migration
all_disp_all2[Migration %in% c("directional migratory","dispersive migratory")]$Migration_str<-"migratory"
colorBlindBlack8  <- c("resident"="#E69F00", 
                       "nomadic"="#56B4E9",
                       "migratory"="#009E73", 
                       "unknown"="#D55E00")


p<-ggplot()+geom_histogram(data=all_disp_all2, 
                           aes(x=mean_disp/1000, fill=Migration_str), 
                           alpha=0.5, bins = 50)+
  facet_wrap(~group, scale="free")+ylab("Number of species")+xlab("mean dispersal distance")+
  scale_fill_manual(values=colorBlindBlack8)+
  theme_bw()
p
ggsave(p, filename="../../Figures/dispersal_distance/dispersal_distance_migration.png", width=12, height=6)
all_disp_all2[,.(mean_disp=mean(mean_disp),
                 sd_disp=sd(mean_disp),
                 max=max(mean_disp),
                 min=min(mean_disp)), by=list(group,Migration_str)]
range()
ggsave(p, filename="../../Figures/dispersal_distance/dispersal_distance_diet.png", width=12, height=6)
