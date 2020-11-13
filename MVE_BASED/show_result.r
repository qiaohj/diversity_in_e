library(raster)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
group<-"Mammals"
sp<-"Akodon_albiventer"
GCM<-"UKESM1"
SSP<-"SSP585"
M<-5
N<-1
target<-sprintf("../../Figures/Dispersal_Example/%s/%s/%s/%s/%d_%d", group, sp, GCM, SSP, M, N)

dir.create(target, showWarnings = T, recursive = T)

print("Init all potential distributions")
distributoins<-list()
enm_folter<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s/predict", group, sp)
for (year in predict_range){
  #print(year)
  env_item<-readRDS(sprintf("%s/%s_%s_%d.rda", enm_folder, GCM, SSP, year))
  distributoins[[as.character(year)]]<-env_item
}

dispersal_log<-readRDS(sprintf("%s/../dispersal/%s_%s_%d_%d.rda", enm_folter, GCM, SSP, M, N))
if (F){
  for (year in predict_range){
    env_item<-distributoins[[as.character(year)]]
    p<-ggplot(env_item, aes(x=x, y=y), color="grey") + geom_point() +
      geom_point(data=dispersal_log%>%dplyr::filter(YEAR==year), aes(x=x, y=y), color="red")+
      theme_bw()+ggtitle(year)
    #print(p)
    ggsave(p, file=sprintf("%s/%d.png", target, year))
    
  }
}