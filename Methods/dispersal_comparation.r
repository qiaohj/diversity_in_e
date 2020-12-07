library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(ggpubr)
library(data.table)
library(batchtools)
group<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
example_sp<-"Dendropsophus_walfordi"
target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
target<-sprintf("%s/dispersal_1", target_folder)
dispersal_1<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1"))
target<-sprintf("%s/dispersal_5", target_folder)
dispersal_5<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1"))
start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
xrange<-c(-8000000, -3000000)
yrange<-c(-4500000, 1500000)
predict_range<-c(2021:2100)
year<-2100
p<-ggplot()+
  xlim(xrange)+
  ylim(yrange)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")

for (year in predict_range){
  print(year)
  dis1<-dispersal_1%>%dplyr::filter(YEAR==year)
  p1<-p+geom_tile(data=dis1, aes(x=x, y=y, fill=factor(exposure)))+
    geom_text(data=dis1, aes(x=x, y=y, label=exposure), size=2)+
    ggtitle(paste(year, ", Exposure year: 1"))
  dis5<-dispersal_5%>%dplyr::filter(YEAR==year)
  p5<-p+geom_tile(data=dis5, aes(x=x, y=y, fill=factor(exposure)))+
    geom_text(data=dis5, aes(x=x, y=y, label=exposure), size=2)+
    ggtitle(paste(year, ", Exposure year: 5"))
  gg<-ggarrange(p1, p5, nrow = 1, ncol=2)
  ggsave(gg, filename=sprintf("../../Figures/Methods/dispersal_comparison/%d.png", year), width=8, height=5)
}


#ffmpeg -r 2 -start_number 2021 -i %04d.png -y ../dispersal_comparison.mp4
#ffmpeg -r 2 -start_number 2021 -i %04d.png -y ../dispersal_comparison.gif
