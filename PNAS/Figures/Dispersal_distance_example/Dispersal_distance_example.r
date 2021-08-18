library(data.table)
library(ggplot2)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
species_list<-c("Clibanornis rubiginosus", "Dysithamnus mentalis",
                "Myrmotherula schisticolor", "Myiobius sulphureipygius",
                "Henicorhina leucosticta", "Microcerculus marginatus",
                "Habia rubica", "Eucometis penicillata")

species_list<-c("Clibanornis rubiginosus", "Dysithamnus mentalis",
                "Myrmotherula schisticolor", "Myiobius sulphureipygius",
                 "Microcerculus marginatus"
                )
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
sp<-species_list[5]

sub_bird_disp<-bird_disp[iucn_name %in% species_list]
sub_bird_disp$label<-sprintf("%s, %.2fkm", sub_bird_disp$iucn_name, sub_bird_disp$estimated_disp)

p<-ggplot(bird_disp)+geom_point(aes(x=Diet, y=estimated_disp))+
  geom_point(data=sub_bird_disp, aes(x=Diet, y=estimated_disp), color="red")+
  theme_bw()

ggsave(p, filename="../../Figures/Dispersal_distance_example/estimated_dispersal_distance.png")
bird_disp_item<-bird_disp[iucn_name==sp]

exposure<-c(0, 5)
tempelate<-"../../Objects/Dispersal/Birds/%s/UKESM1_SSP245_%d_dispersal_1.rda"

e<-5


all_item<-list()
i=1
for (i in c(1:nrow(sub_bird_disp))){
  item_sp<-sub_bird_disp[i,]
  sp<-item_sp$iucn_name
  if (sp %in% c("Henicorhina leucosticta", "Habia rubica", "Eucometis penicillata")){
    next()
  }
  for (e in exposure){
    item<-readRDS(sprintf(tempelate, gsub(" ", "_", sp), e))
    bird_disp_item<-bird_disp[iucn_name==sp]
    year=2022
    
    for (year in c(2022:2100)){
      print(paste(sp, e, year))
      y1<-year-1
      y2<-year  
      item1<-item[[as.character(y1)]]
      item2<-item[[as.character(y2)]]
      me_item<-merge(item1, item2, by=c("x", "y", "mask_100km"))
      
      me_item$dispersal_dist<-me_item$accumulative_disp.y-me_item$accumulative_disp.x
      me_item<-me_item[dispersal_dist>0]
      cols<-c("x", "y", "mask_100km",  "dispersal_dist")
      me_item<-me_item[, ..cols]
      me_item$year<-y2
      me_item$iucn_name<-sp
      me_item$exposure<-e
      me_item$label<-item_sp$label
      all_item[[paste(sp, e, year)]]<-me_item
    }
  }
}
all_item<-rbindlist(all_item)
all_item_se<-all_item[, .(mean_dist=mean(dispersal_dist), sd_dist=sd(dispersal_dist)), 
                      by=c("iucn_name", "year", "exposure", "label")]

p1<-ggplot(all_item)+
  geom_histogram(aes(x=dispersal_dist/1000, fill=factor(exposure)), bins=50)+
  theme_bw()+
  facet_wrap(~label, nrow=1)+
  labs(x="Mean dispersal distance (km)", y="count", fill="Exposure")

p1
p2<-ggplot(all_item_se)+
  geom_line(aes(x=year,  y=mean_dist/1000, color=factor(exposure)))+
  #geom_errorbar(aes(x=year, ymin=mean_dist-sd_dist, ymax=mean_dist+sd_dist))+
  theme_bw()+
  facet_wrap(~label, nrow=1)+
  labs(x="YEAR", y="Mean dispersal distance (km)", color="Exposure")

p<-ggarrange(p1, p2, nrow=2)
p
ggsave(p, filename="../../Figures/Dispersal_distance_example/Dispersal_distance_example.png", width=15, height=6)
ggsave(p, filename="../../Figures/Dispersal_distance_example/Dispersal_distance_example.pdf", width=15, height=6)



examples<-quantile(bird_disp$estimated_disp, c(0.1, 0.5, 0.9))


examples<-bird_disp[estimated_disp %in% examples][c(1:3),]


exposure<-c(0, 5)
tempelate<-"../../Objects/Dispersal/Birds/%s/UKESM1_SSP245_%d_dispersal_1.rda"
examples$label<-sprintf("%s, %.2fkm", examples$iucn_name, examples$estimated_disp)
all_item<-list()
i=1
for (i in c(1:3)){
  iii<-examples[i,]
  sp<-iii$iucn_name
  if (sp %in% c("Henicorhina leucosticta", "Habia rubica", "Eucometis penicillata")){
    next()
  }
  for (e in exposure){
    item<-readRDS(sprintf(tempelate, gsub(" ", "_", sp), e))
    bird_disp_item<-bird_disp[iucn_name==sp]
    year=2022
    
    for (year in c(2022:2100)){
      print(paste(sp, e, year))
      y1<-year-1
      y2<-year  
      item1<-item[[as.character(y1)]]
      item2<-item[[as.character(y2)]]
      me_item<-merge(item1, item2, by=c("x", "y", "mask_100km"))
      
      me_item$dispersal_dist<-me_item$accumulative_disp.y-me_item$accumulative_disp.x
      me_item<-me_item[dispersal_dist>0]
      cols<-c("x", "y", "mask_100km",  "dispersal_dist")
      me_item<-me_item[, ..cols]
      me_item$year<-y2
      me_item$iucn_name<-sp
      me_item$exposure<-e
      me_item$label<-iii$label
      all_item[[paste(sp, e, year)]]<-me_item
    }
  }
}
all_item<-rbindlist(all_item)
all_item_se<-all_item[, .(mean_dist=mean(dispersal_dist), sd_dist=sd(dispersal_dist)), 
                      by=c("iucn_name", "year", "exposure", "label")]

p1<-ggplot(all_item)+
  geom_histogram(aes(x=dispersal_dist/1000, fill=factor(exposure)), bins=50)+
  theme_bw()+
  facet_wrap(~label, nrow=1)+
  labs(x="Mean dispersal distance (km)", y="count", fill="Exposure")

p1
p2<-ggplot(all_item_se)+
  geom_line(aes(x=year,  y=mean_dist/1000, color=factor(exposure)))+
  #geom_errorbar(aes(x=year, ymin=mean_dist-sd_dist, ymax=mean_dist+sd_dist))+
  theme_bw()+
  facet_wrap(~label, nrow=1)+
  labs(x="YEAR", y="Mean dispersal distance (km)", color="Exposure")

p<-ggarrange(p1, p2, nrow=2)
p
ggsave(p, filename="../../Figures/Dispersal_distance_example/Dispersal_distance_example_2.png", width=8, height=4)
ggsave(p, filename="../../Figures/Dispersal_distance_example/Dispersal_distance_example_2.pdf", width=8, height=4)


          