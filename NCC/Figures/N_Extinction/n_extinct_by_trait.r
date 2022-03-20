library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/colors.r")
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
colnames(bird_disp)[5]<-"sp"
bird_disp$sp<-gsub(" ", "_", bird_disp$sp)
mammal_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
colnames(mammal_disp)[c(1)]<-"sp"

sp_dis_all_l<-list()
for (exposure in c(0, 5)){
  rda<-sprintf("../../Figures/N_Extinction/sp_dis_all_%d_10km_2_100km.rda", exposure)
  print(paste("Reading", rda))
  sp_dis_all<-readRDS(rda)
  sp_dis_all$exposure<-exposure
  sp_dis_all_l[[length(sp_dis_all_l)+1]]<-sp_dis_all
}
sp_dis_all_l<-rbindlist(sp_dis_all_l)
df_with_family<-data.table(read.csv("../../Data/Dispersal_distance/bird.csv", stringsAsFactors = F))

cols<-c("max_dis", "HWI", "log_body_mass", "Diet", "sp", "Order.x", "Family.name")
model_df_birds<-df_with_family[, ..cols]
min_HWI<-min(model_df_birds$HWI)
max_HWI<-max(model_df_birds$HWI)
min_body_mass<-min(model_df_birds$log_body_mass)
max_body_mass<-max(model_df_birds$log_body_mass)

sp_list<-sp_dis_all_l[, .(S=.N), by=list(sp, group)]
bird_disp_filter<-bird_disp[sp %in% sp_list$sp]
dim(bird_disp_filter[between(HWI, min_HWI, max_HWI)])
dim(bird_disp_filter[between(log_body_mass, min_body_mass, max_body_mass)|between(HWI, min_HWI, max_HWI)])
sp1<-bird_disp_filter[HWI<min_HWI]$sp
length(sp1)
sp2<-bird_disp_filter[log_body_mass<min_body_mass]$sp
length(sp2)
length(sp1[sp1 %in% sp2])
sp1<-bird_disp_filter[HWI>max_HWI]$sp
sp2<-bird_disp_filter[log_body_mass>max_body_mass]$sp
length(sp1[sp1 %in% sp2])

dim(bird_disp_filter[HWI>=max_HWI&log_body_mass>=min_body_mass])

df_with_family<-data.table(read.csv("../../Data/Dispersal_distance/mammal.csv", stringsAsFactors = F))

model_df_mammal<-df_with_family
min_body_mass<-min(model_df_mammal$log_body_mass, na.rm=T)
max_body_mass<-max(model_df_mammal$log_body_mass, na.rm=T)
mammal_disp$sp<-gsub(" ", "_", mammal_disp$sp)
mammals_disp_filter<-mammal_disp[sp %in% sp_list$sp]
dim(mammals_disp_filter[between(log_body_mass, min_body_mass, max_body_mass)])
sp2<-mammals_disp_filter[log_body_mass<min_body_mass]$sp
length(sp2)
sp2<-mammals_disp_filter[log_body_mass>max_body_mass]$sp
length(sp2)


dim(bird_disp_filter[HWI>=max_HWI&log_body_mass>=min_body_mass])


extinct<-sp_dis_all_l[N_type=="EXTINCT"]
extinct<-extinct[, .(N=.N), by=c("sp", "GCM", "SSP", "M", "exposure", "N_type", "group")]
table(extinct$group)
#Birds
extinct_bird<-merge(extinct, bird_disp, by="sp")
bird_disp_target<-bird_disp[sp %in% sp_dis_all_l$sp]
#For diet
bird_disp_diet<-bird_disp_target[, .(N_all_sp=.N), by=list(Diet)]
extinct_bird_diet<-extinct_bird[, .(N=.N), by=list(GCM, SSP, M, exposure, Diet)]
full_Com<-data.table(expand.grid(SSP=c("SSP245", "SSP585", "SSP119"),
                      GCM=c("UKESM1", "MRI-ESM2-0", "EC-Earth3-Veg"),
                      M=c(0, 1),
                      exposure=c(0, 5),
                      Diet=unique(extinct_bird_diet$Diet),
                      stringsAsFactors = F))
extinct_bird_diet<-merge(full_Com, extinct_bird_diet, by=c("SSP", "GCM", "M", "exposure", "Diet"), all=T)
extinct_bird_diet[is.na(N)]$N<-0
extinct_bird_diet<-merge(extinct_bird_diet, bird_disp_diet, by="Diet")
extinct_bird_diet$extinct_proportion<-extinct_bird_diet$N/extinct_bird_diet$N_all_sp * 100
extinct_bird_diet<-extinct_bird_diet[, .(N=mean(N), sd_N=sd(N),
                                         extinct_proportion=mean(extinct_proportion), 
                                         sd_extinct_proportion=sd(extinct_proportion)), 
                                     by=list(SSP, M, exposure, Diet)]
extinct_bird_diet$dispersal<-ifelse(extinct_bird_diet$M==0, "no dispersal", "with dispersal")
extinct_bird_diet$exposure<-ifelse(extinct_bird_diet$exposure==0, " no climate resilience", "climate resilience")


p<-ggplot(extinct_bird_diet)+
  geom_bar(aes(x=Diet, y=extinct_proportion,fill=SSP), stat="identity", position="dodge")+
  geom_errorbar(aes(x=Diet, 
                    ymin=extinct_proportion-sd_extinct_proportion,
                    ymax=extinct_proportion+sd_extinct_proportion,
                    group=SSP), width=0.5, position=position_dodge(.9))+
  scale_fill_manual(values=color_ssp)+
  facet_grid(exposure~dispersal)+
  theme_bw()+
  xlab("Diet")+ylab("Mean extinction proportion")
ggsave(p, filename="../../Figures/N_Extinction/Extinction_bird_diet_10km_2_100km_extinct_only.png", width=10, height=6)


#For $Migration

bird_disp_target$Migration<-bird_disp_target$Migration_3
bird_disp_target[Migration %in% c("dispersive migratory", "directional migratory")]$Migration<-"migratory"
bird_disp_target[is.na(Migration)]$Migration<-"unknown"
bird_disp_Migration<-bird_disp_target[, .(N_all_sp=.N), by=list(Migration)]


extinct_bird$Migration<-extinct_bird$Migration_3
extinct_bird[Migration %in% c("dispersive migratory", "directional migratory")]$Migration<-"migratory"
extinct_bird[is.na(Migration)]$Migration<-"unknown"

extinct_bird_Migration<-extinct_bird[, .(N=.N), by=list(GCM, SSP, M, exposure, Migration)]
full_Com<-data.table(expand.grid(SSP=c("SSP245", "SSP585", "SSP119"),
                                 GCM=c("UKESM1", "MRI-ESM2-0", "EC-Earth3-Veg"),
                                 M=c(0, 1),
                                 exposure=c(0, 5),
                                 Migration=unique(extinct_bird_Migration$Migration),
                                 stringsAsFactors = F))
extinct_bird_Migration<-merge(full_Com, extinct_bird_Migration, by=c("SSP", "GCM", "M", "exposure", "Migration"), all=T)
extinct_bird_Migration[is.na(N)]$N<-0
extinct_bird_Migration<-merge(extinct_bird_Migration, bird_disp_Migration, by="Migration")
extinct_bird_Migration$extinct_proportion<-extinct_bird_Migration$N/extinct_bird_Migration$N_all_sp * 100
extinct_bird_Migration<-extinct_bird_Migration[, .(N=mean(N), sd_N=sd(N),
                                         extinct_proportion=mean(extinct_proportion), 
                                         sd_extinct_proportion=sd(extinct_proportion)), 
                                     by=list(SSP, M, exposure, Migration)]
extinct_bird_Migration$dispersal<-ifelse(extinct_bird_Migration$M==0, "no dispersal", "with dispersal")
extinct_bird_Migration$exposure<-ifelse(extinct_bird_Migration$exposure==0, " no climate resilience", "climate resilience")


p<-ggplot(extinct_bird_Migration)+
  geom_bar(aes(x=Migration, y=extinct_proportion,fill=SSP), stat="identity", position="dodge")+
  geom_errorbar(aes(x=Migration, 
                    ymin=extinct_proportion-sd_extinct_proportion,
                    ymax=extinct_proportion+sd_extinct_proportion,
                    group=SSP), width=0.5, position=position_dodge(.9))+
  scale_fill_manual(values=color_ssp)+
  facet_grid(exposure~dispersal)+
  theme_bw()+
  xlab("Migration type")+ylab("Mean extinction proportion")
p
ggsave(p, filename="../../Figures/N_Extinction/Extinction_bird_Migration_10km_2_100km_extinct_only.png", width=10, height=6)



#For HWI
range(bird_disp_target$HWI)
HWI_cuts<-seq(from=0, to=75, by=75/100)
bird_disp_target$HWI_bins<-Hmisc::cut2(bird_disp_target$HWI, cuts=HWI_cuts)
bird_disp_HWI<-bird_disp_target[, .(N_all_sp=.N), by=list(HWI_bins)]
extinct_bird$HWI_bins<-Hmisc::cut2(extinct_bird$HWI, cuts=HWI_cuts)
extinct_bird_HWI<-extinct_bird[, .(N=.N), by=list(GCM, SSP, M, exposure, HWI_bins)]
full_Com<-data.table(expand.grid(SSP=c("SSP245", "SSP585", "SSP119"),
                                 GCM=c("UKESM1", "MRI-ESM2-0", "EC-Earth3-Veg"),
                                 M=c(0, 1),
                                 exposure=c(0, 5),
                                 HWI_bins=unique(extinct_bird_HWI$HWI_bins),
                                 stringsAsFactors = F))
extinct_bird_HWI<-merge(full_Com, extinct_bird_HWI, by=c("SSP", "GCM", "M", "exposure", "HWI_bins"), all=T)
extinct_bird_HWI[is.na(N)]$N<-0
extinct_bird_HWI<-merge(extinct_bird_HWI, bird_disp_HWI, by="HWI_bins")
extinct_bird_HWI$extinct_proportion<-extinct_bird_HWI$N/extinct_bird_HWI$N_all_sp * 100
#extinct_bird_HWI<-extinct_bird_HWI[, .(N=mean(N), sd_N=sd(N),
#                                         extinct_proportion=mean(extinct_proportion), 
#                                         sd_extinct_proportion=sd(extinct_proportion)), 
#                                     by=list(SSP, M, exposure, HWI_bins)]
extinct_bird_HWI$dispersal<-ifelse(extinct_bird_HWI$M==0, "no dispersal", "with dispersal")
extinct_bird_HWI$exposure<-ifelse(extinct_bird_HWI$exposure==0, " no climate resilience", "climate resilience")

extinct_bird_HWI$HWI_v<-HWI_cuts[as.numeric(extinct_bird_HWI$HWI_bins)]
p<-ggplot(extinct_bird_HWI)+
  geom_point(aes(x=HWI_v, y=extinct_proportion,color=SSP))+
  geom_smooth(aes(x=HWI_v, y=extinct_proportion,color=SSP), method="loess")+
  #geom_errorbar(aes(x=HWI_v, 
  #                  ymin=extinct_proportion-sd_extinct_proportion,
  #                  ymax=extinct_proportion+sd_extinct_proportion,
  #                  group=SSP), width=0.5, position=position_dodge(.9))+
  scale_color_manual(values=color_ssp)+
  facet_grid(exposure~dispersal)+
  theme_bw()+
  xlab("HWI")+ylab("Mean extinction proportion")
p
ggsave(p, filename="../../Figures/N_Extinction/Extinction_bird_hwi_10km_2_100km_extinct_only.png", width=10, height=6)


#For body mass
range(bird_disp_target$body_mass)
body_mass_cuts<-seq(from=0, to=50, by=160/100)
bird_disp_target$body_mass_bins<-Hmisc::cut2(bird_disp_target$body_mass, cuts=body_mass_cuts)
bird_disp_body_mass<-bird_disp_target[, .(N_all_sp=.N), by=list(body_mass_bins)]
extinct_bird$body_mass_bins<-Hmisc::cut2(extinct_bird$body_mass, cuts=body_mass_cuts)
extinct_bird_body_mass<-extinct_bird[, .(N=.N), by=list(GCM, SSP, M, exposure, body_mass_bins)]
full_Com<-data.table(expand.grid(SSP=c("SSP245", "SSP585", "SSP119"),
                                 GCM=c("UKESM1", "MRI-ESM2-0", "EC-Earth3-Veg"),
                                 M=c(0, 1),
                                 exposure=c(0, 5),
                                 body_mass_bins=unique(extinct_bird_body_mass$body_mass_bins),
                                 stringsAsFactors = F))
extinct_bird_body_mass<-merge(full_Com, extinct_bird_body_mass, by=c("SSP", "GCM", "M", "exposure", "body_mass_bins"), all=T)
extinct_bird_body_mass[is.na(N)]$N<-0
extinct_bird_body_mass<-merge(extinct_bird_body_mass, bird_disp_body_mass, by="body_mass_bins")
extinct_bird_body_mass$extinct_proportion<-extinct_bird_body_mass$N/extinct_bird_body_mass$N_all_sp * 100
#extinct_bird_body_mass<-extinct_bird_body_mass[, .(N=mean(N), sd_N=sd(N),
#                                       extinct_proportion=mean(extinct_proportion), 
#                                       sd_extinct_proportion=sd(extinct_proportion)), 
#                                   by=list(SSP, M, exposure, body_mass_bins)]
extinct_bird_body_mass$dispersal<-ifelse(extinct_bird_body_mass$M==0, "no dispersal", "with dispersal")
extinct_bird_body_mass$exposure<-ifelse(extinct_bird_body_mass$exposure==0, " no climate resilience", "climate resilience")

extinct_bird_body_mass$body_mass_v<-body_mass_cuts[as.numeric(extinct_bird_body_mass$body_mass_bins)]
p<-ggplot(extinct_bird_body_mass)+
  geom_point(aes(x=body_mass_v, y=extinct_proportion,color=SSP))+
  geom_smooth(aes(x=body_mass_v, y=extinct_proportion,color=SSP), method="loess")+
  #geom_errorbar(aes(x=body_mass_v, 
  #                  ymin=extinct_proportion-sd_extinct_proportion,
  #                  ymax=extinct_proportion+sd_extinct_proportion,
  #                  group=SSP), width=0.5, position=position_dodge(.9))+
  scale_color_manual(values=color_ssp)+
  facet_grid(exposure~dispersal)+
  theme_bw()+
  xlab("Body mass")+ylab("Mean extinction proportion")
p
ggsave(p, filename="../../Figures/N_Extinction/Extinction_bird_body_mass_10km_2_100km_extinct_only.png", width=10, height=6)

#For estimated max natal dispersal distance
hist(bird_disp_target$estimated_disp)
estimated_disp_cuts<-seq(from=0, to=500, by=500/100)
bird_disp_target$estimated_disp_bins<-Hmisc::cut2(bird_disp_target$estimated_disp, cuts=estimated_disp_cuts)
bird_disp_estimated_disp<-bird_disp_target[, .(N_all_sp=.N), by=list(estimated_disp_bins)]
extinct_bird$estimated_disp_bins<-Hmisc::cut2(extinct_bird$estimated_disp, cuts=estimated_disp_cuts)
extinct_bird_estimated_disp<-extinct_bird[, .(N=.N), by=list(GCM, SSP, M, exposure, estimated_disp_bins)]
full_Com<-data.table(expand.grid(SSP=c("SSP245", "SSP585", "SSP119"),
                                 GCM=c("UKESM1", "MRI-ESM2-0", "EC-Earth3-Veg"),
                                 M=c(0, 1),
                                 exposure=c(0, 5),
                                 estimated_disp_bins=unique(extinct_bird_estimated_disp$estimated_disp_bins),
                                 stringsAsFactors = F))
extinct_bird_estimated_disp<-merge(full_Com, extinct_bird_estimated_disp, 
                                   by=c("SSP", "GCM", "M", "exposure", "estimated_disp_bins"), all=T)
extinct_bird_estimated_disp[is.na(N)]$N<-0
extinct_bird_estimated_disp<-merge(extinct_bird_estimated_disp, bird_disp_estimated_disp, by="estimated_disp_bins")
extinct_bird_estimated_disp$extinct_proportion<-extinct_bird_estimated_disp$N/extinct_bird_estimated_disp$N_all_sp * 100

#extinct_bird_estimated_disp<-extinct_bird_estimated_disp[, .(N=mean(N), sd_N=sd(N),
#                                                   extinct_proportion=mean(extinct_proportion), 
#                                                   sd_extinct_proportion=sd(extinct_proportion)), 
#                                               by=list(SSP, M, exposure, estimated_disp_bins)]
extinct_bird_estimated_disp$dispersal<-ifelse(extinct_bird_estimated_disp$M==0, "no dispersal", "with dispersal")
extinct_bird_estimated_disp$exposure<-ifelse(extinct_bird_estimated_disp$exposure==0, " no climate resilience", "climate resilience")

extinct_bird_estimated_disp$estimated_disp_v<-estimated_disp_cuts[as.numeric(extinct_bird_estimated_disp$estimated_disp_bins)]
p<-ggplot(extinct_bird_estimated_disp)+
  geom_point(aes(x=estimated_disp_v, y=extinct_proportion,color=SSP))+
  geom_smooth(aes(x=estimated_disp_v, y=extinct_proportion,color=SSP), method="loess")+
  #geom_errorbar(aes(x=estimated_disp_v, 
  #                  ymin=extinct_proportion-sd_extinct_proportion,
  #                  ymax=extinct_proportion+sd_extinct_proportion,
  #                  group=SSP), width=0.5, position=position_dodge(.9))+
  scale_color_manual(values=color_ssp)+
  facet_grid(exposure~dispersal)+
  theme_bw()+
  xlab("Estimated maximum natal dispersal distance")+ylab("Mean extinction proportion")
p
ggsave(p, filename="../../Figures/N_Extinction/Extinction_bird_estimated_disp_10km_2_100km_extinct_only.png", width=10, height=6)
