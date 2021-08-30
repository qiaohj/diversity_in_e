library(data.table)
library(ggplot2)
library(Hmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

islands<-readRDS("../../Objects/Island/species_island.rda")
mountains<-readRDS("../../Objects/Island/species_mountain.rda")
birds<-readRDS("../../Objects/IUCN_List/Birds_df_with_family.rda")
mammals<-readRDS("../../Objects/IUCN_List/Mammals_df_with_family.rda")

colnames(mammals)[8]<-"Order"
colnames<-c("Order", "family", "SP", "estimated_disp", "N_CELL")
all_sp<-rbindlist(list(birds[, ..colnames], mammals[, ..colnames]))
colnames(all_sp)[3]<-"sp"
all_sp$sp<-gsub(" ", "_", all_sp$sp)
all_sp<-merge(all_sp, islands, by="sp")
all_sp[group=="Bird"]$group<-"Birds"
all_sp[(type=="continent")&(sp %in% mountains[type=="mountain"]$sp)]$type<-"mountain"
all_sp[(type=="continent")&(sp %in% mountains[type=="plain"]$sp)]$type<-"plain"

table(all_sp$type)

mean_disp<-readRDS("../../Objects/Dispersal_distances/all_mean_disp_dist.rda")
mean_disp$sp<-gsub(" ", "_", mean_disp$iucn_name)

mean_disp_all<-merge(mean_disp, all_sp, by=c("sp", "group"))
mean_disp_all$N_CELL<-ceiling(mean_disp_all$N_CELL/100)
mean_disp_all$mean_dist<-mean_disp_all$mean_dist/1000
#plot(mean_disp_all$N_CELL, mean_disp_all$continent+mean_disp_all$island)

mean_disp_all$mean_dist_raw<-mean_disp_all$mean_dist
#mean_disp_all[N_Reset_Cell==0]
mean_disp_all$mean_dist<-mean_disp_all$mean_dist_raw + mean_disp_all$N_New_Cell /
  (mean_disp_all$N_New_Cel+mean_disp_all$N_Reset_Cell+1)
mean_disp_all_se<-mean_disp_all[, .(mean_dist=mean(mean_dist),
                                    estimated_disp=mean(estimated_disp),
                                    N_CELL=mean(N_CELL)),
                                by=list(group, exposure, SSP, type)]


island_species<-mean_disp_all[type=="island"]
range(island_species$N_CELL)
range(island_species$estimated_disp)

mountain_species<-mean_disp_all[type=="mountain"]
range(mountain_species$N_CELL)
range(mountain_species$estimated_disp)

plain_species<-mean_disp_all[type=="plain"]
range(plain_species$N_CELL)
range(plain_species$estimated_disp)

cuts<-seq(0, 350, by=50)
disp_cuts<-seq(0, 120, by=10)
colnames<-c("sp", "N_CELL_cuts", "disp_cuts", "type", "group")

island_species$N_CELL_cuts<-cut2(island_species$N_CELL, cuts)
island_species$disp_cuts<-cut2(island_species$estimated_disp, disp_cuts)
island_species<-island_species[, ..colnames]
island_species<-unique(island_species)
island_species$N_CELL_cuts<-as.numeric(island_species$N_CELL_cuts)
table(island_species$N_CELL_cuts)
island_species$disp_cuts<-as.numeric(island_species$disp_cuts)
table(island_species$disp_cuts)

mountain_species$N_CELL_cuts<-cut2(mountain_species$N_CELL, cuts)
mountain_species$disp_cuts<-cut2(mountain_species$estimated_disp, disp_cuts)
mountain_species<-mountain_species[, ..colnames]
mountain_species<-unique(mountain_species)
mountain_species$N_CELL_cuts<-as.numeric(mountain_species$N_CELL_cuts)
table(mountain_species$N_CELL_cuts)
mountain_species$disp_cuts<-as.numeric(mountain_species$disp_cuts)
table(mountain_species$disp_cuts)

plain_species$N_CELL_cuts<-cut2(plain_species$N_CELL, cuts)
plain_species$disp_cuts<-cut2(plain_species$estimated_disp, disp_cuts)
plain_species<-plain_species[, ..colnames]
plain_species<-unique(plain_species)
plain_species$N_CELL_cuts<-as.numeric(plain_species$N_CELL_cuts)
table(plain_species$N_CELL_cuts)
plain_species$disp_cuts<-as.numeric(plain_species$disp_cuts)
table(plain_species$disp_cuts)

cuts_methods<-"dispersal"
cuts_methods<-"none"
cuts_methods<-"n_cell"

if (cuts_methods=="dispersal"){
  #resample
  table(mountain_species$disp_cuts)
  table(island_species$disp_cuts)
  table(plain_species$N_CELL_cuts)
  N=100
  all_cuts<-as.character(unique(mountain_species$disp_cuts))
  i=1
  island_species_sampled<-list()
  plain_species_sampled<-list()
  mountain_species_sampled<-list()
  for (i in c(4:13)){
    #N<-nrow(island_species[N_CELL_cuts==i])
    item<-mountain_species[disp_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    mountain_species_sampled[[i]]<-item
    
    item<-plain_species[disp_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    plain_species_sampled[[i]]<-item
    
    item<-island_species[disp_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    island_species_sampled[[i]]<-item
  }
  mountain_species_sampled<-rbindlist(mountain_species_sampled)
  table(mountain_species_sampled$N_CELL_cuts)
  plain_species_sampled<-rbindlist(plain_species_sampled)
  table(plain_species_sampled$N_CELL_cuts)
  island_species_sampled<-rbindlist(island_species_sampled)
  table(island_species_sampled$N_CELL_cuts)
  sampled_species<-rbindlist(list(mountain_species_sampled, plain_species_sampled, island_species_sampled))
}

if (cuts_methods=="n_cell"){
  #resample
  table(mountain_species$N_CELL_cuts)
  table(island_species$N_CELL_cuts)
  table(plain_species$N_CELL_cuts)
  N=10
  all_cuts<-as.character(unique(mountain_species$N_CELL_cuts))
  i=1
  island_species_sampled<-list()
  plain_species_sampled<-list()
  mountain_species_sampled<-list()
  for (i in c(1:7)){
    #N<-nrow(island_species[N_CELL_cuts==i])
    item<-mountain_species[N_CELL_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    mountain_species_sampled[[i]]<-item
    
    item<-plain_species[N_CELL_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    plain_species_sampled[[i]]<-item
    
    item<-island_species[N_CELL_cuts==i]
    if (nrow(item)==0){
      next()
    }
    if (nrow(item)>=N){
      item<-item[sample(nrow(item), N),]
    }else{
      item<-item[sample(nrow(item), N, replace=T),]
    }
    island_species_sampled[[i]]<-item
  }
  mountain_species_sampled<-rbindlist(mountain_species_sampled)
  table(mountain_species_sampled$N_CELL_cuts)
  plain_species_sampled<-rbindlist(plain_species_sampled)
  table(plain_species_sampled$N_CELL_cuts)
  island_species_sampled<-rbindlist(island_species_sampled)
  table(island_species_sampled$N_CELL_cuts)
  sampled_species<-rbindlist(list(mountain_species_sampled, plain_species_sampled, island_species_sampled))
}

if (cuts_methods=="none"){
  sampled_species<-rbindlist(list(island_species, mountain_species, plain_species))
}

mean_disp_all_se<-mean_disp_all[sp %in% sampled_species$sp, 
                                .(mean_dist=mean(mean_dist),
                                  estimated_disp=mean(estimated_disp),
                                  N_CELL=mean(N_CELL),
                                  sd_mean_dist=sd(mean_dist),
                                  sd_estimated_disp=sd(estimated_disp),
                                  sd_N_CELL=sd(N_CELL),
                                  N_New_Cell=mean(N_New_Cell),
                                  N_Reset_Cell=mean(N_Reset_Cell),
                                  sd_N_New_Cell=sd(N_New_Cell),
                                  sd_N_Reset_Cell=sd(N_Reset_Cell),
                                  mean_N_New_Cell=mean(mean_N_New_Cell),
                                  sd_mean_N_New_Cell=sd(mean_N_New_Cell)),
                                by=list(group, exposure, SSP, type)]
table(sampled_species$N_CELL_cuts)
table(sampled_species$disp_cuts)
#mean_disp_all_se[type=="plain"]$mean_dist<-mean_disp_all_se[type=="plain"]$mean_dist+30
mean_disp_all_se_n_cell<-mean_disp_all[sp %in% sampled_species$sp, 
                                .(mean_dist=mean(mean_dist),
                                  estimated_disp=mean(estimated_disp),
                                  sd_mean_dist=sd(mean_dist),
                                  sd_estimated_disp=sd(estimated_disp)),
                                by=list(exposure, SSP, type, N_CELL)]

ggplot(mean_disp_all_se_n_cell[N_CELL<=500])+geom_point(aes(x=N_CELL, y=mean_dist, color=type))+
  facet_grid(SSP~exposure)+theme_bw()

mean_disp_all_se$exposure<-ifelse(mean_disp_all_se$exposure==0, " no climate reslience", "climate reslience")
source("commonFuns/colors.r")
p<-ggplot(mean_disp_all_se)+
  geom_point(aes(x=type, y=mean_dist+mean_N_New_Cell*3.5, color=group), position=position_dodge(width=0.8))+
  #geom_point(aes(x=type, y=mean_dist, color=group), position=position_dodge(width=0.8))+
  geom_point(aes(x=type, y=estimated_disp, color=group), shape=2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=mean_dist+mean_N_New_Cell*3.5-sd_mean_dist, ymax=mean_dist+mean_N_New_Cell*3.5+sd_mean_dist, 
                    color=group, linetype=" mean dispersal distance"), 
                width=.2, position=position_dodge(width=0.8))+
  geom_errorbar(aes(x=type, ymin=estimated_disp-sd_estimated_disp, ymax=estimated_disp+sd_estimated_disp, color=group, 
                    linetype="estimate max natal dispersal distance"), 
                width=.2, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group", linetype="")+
  xlab("")+ylab("Dispersal distance (km)")+
  facet_grid(exposure~SSP)+theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/dispersal_distance.png", cuts_methods), width=10, height=6)
p<-ggplot(mean_disp_all_se)+
  geom_point(aes(x=type, y=N_New_Cell, color=group), position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=N_New_Cell-sd_N_New_Cell, ymax=N_New_Cell+sd_N_New_Cell, color=group), 
                width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Number of new cells")+
  facet_grid(exposure~SSP)+theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/number_of_new_cells.png", cuts_methods), width=10, height=6)

ggplot(mean_disp_all_se)+
  geom_point(aes(x=type, y=mean_N_New_Cell, color=group))+
  geom_errorbar(aes(x=type, ymin=mean_N_New_Cell-sd_mean_N_New_Cell, ymax=mean_N_New_Cell+sd_mean_N_New_Cell, color=group), 
                width=.1, position = "dodge")+
  facet_grid(SSP~exposure)+theme_bw()

p<-ggplot(mean_disp_all_se)+
  geom_point(aes(x=type, y=N_CELL, color=group), position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=N_CELL-sd_N_CELL, ymax=N_CELL+sd_N_CELL, color=group), 
                width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Number of cells")+
  facet_grid(exposure~SSP)+theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/number_of_cells.png", cuts_methods), width=10, height=6)

ggplot(mean_disp_all_se)+
  geom_point(aes(x=type, y=estimated_disp, color=group), position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=estimated_disp-sd_estimated_disp, ymax=estimated_disp+sd_estimated_disp, color=group), 
                width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("estimated max natal dispersal distance")+
  facet_grid(exposure~SSP)+theme_bw()


#extinct
extinct_sp<-readRDS("../../Figures/N_Extinction/sp_dis_all_0.rda")
extinct_sp<-extinct_sp[(year==2100)&(M==1)]
colnames<-c("sp", "GCM", "SSP", "N_type", "group")
extinct_sp<-extinct_sp[, ..colnames]
extinct_sp$exposure<-" no climate reslience"
extinct_sp<-extinct_sp[N_type=="EXTINCT"]

extinct_sp2<-readRDS("../../Figures/N_Extinction/sp_dis_all_5.rda")
extinct_sp2<-extinct_sp2[(year==2100)&(M==1)]
colnames<-c("sp", "GCM", "SSP", "N_type", "group")
extinct_sp2<-extinct_sp2[, ..colnames]
extinct_sp2$exposure<-"climate reslience"
extinct_sp2<-extinct_sp2[N_type=="EXTINCT"]

extinct_sp<-rbindlist(list(extinct_sp, extinct_sp2))

extinct_sp$type<-"mixed"
extinct_sp[sp %in% island_species$sp]$type<-"island"
extinct_sp[sp %in% mountain_species$sp]$type<-"mountain"
extinct_sp[sp %in% plain_species$sp]$type<-"plain"
table(extinct_sp$type)

extinct_sp_group<-extinct_sp[, .(N_Extinct=.N),
                          by=c("GCM", "SSP", "exposure", "type", "group")]
extinct_sp_group<-extinct_sp_group[type!="mixed"]
extinct_sp_group$N_Species<-0
extinct_sp_group[type=="island"]$N_Species<-nrow(island_species)
extinct_sp_group[type=="mountain"]$N_Species<-nrow(mountain_species)
extinct_sp_group[type=="plain"]$N_Species<-nrow(plain_species)
extinct_sp_group$extinct_proportion<-extinct_sp_group$N_Extinct/extinct_sp_group$N_Species

extinct_sp_group_se<-extinct_sp_group[, .(N_Extinct=mean(N_Extinct),
                                          sd_N_Extinct=sd(N_Extinct),
                                          extinct_proportion=mean(extinct_proportion),
                                          sd_extinct_proportion=sd(extinct_proportion)),
                                      by=c("SSP", "exposure", "type", "group")]



p<-ggplot(extinct_sp_group_se)+
  geom_point(aes(x=type, y=extinct_proportion, color=group), position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=extinct_proportion-sd_extinct_proportion, 
                    ymax=extinct_proportion+sd_extinct_proportion, color=group), 
                width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Extinct proportion")+
  facet_grid(exposure~SSP, scale="free")+theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/extinct_proportion_by_group.png", cuts_methods), width=10, height=6)

extinct_sp_group_se<-extinct_sp_group[, .(N_Extinct=mean(N_Extinct),
                                          sd_N_Extinct=sd(N_Extinct),
                                          extinct_proportion=mean(extinct_proportion),
                                          sd_extinct_proportion=sd(extinct_proportion)),
                                      by=c("SSP", "exposure", "type")]



p<-ggplot(extinct_sp_group_se)+
  geom_point(aes(x=type, y=extinct_proportion), position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=type, ymin=extinct_proportion-sd_extinct_proportion, 
                    ymax=extinct_proportion+sd_extinct_proportion), 
                width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Extinct proportion")+
  facet_grid(exposure~SSP, scale="free")+theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/extinct_proportion.png", cuts_methods), width=10, height=6)

#change of env
species_env_change<-readRDS("../../Objects/Island/species_env_change.rda")
species_env_change$type<-"mixed"
species_env_change[sp %in% island_species$sp]$type<-"island"
species_env_change[sp %in% mountain_species$sp]$type<-"mountain"
species_env_change[sp %in% plain_species$sp]$type<-"plain"
species_env_change<-species_env_change[type!="mixed"]
table(species_env_change$type)
species_env_change[group=="Bird"]$group<-"Birds"
species_env_change[group=="Mammal"]$group<-"Mammals"

species_env_change_se<-species_env_change[, .(range=mean(range/1000),
                                              sd_range=sd(range)),
                             by=c("SSP", "type", "group", "var")]


p<-ggplot(species_env_change_se)+
  geom_point(aes(x=type, y=range, color=group), position=position_dodge(width=0.5))+
  #geom_errorbar(aes(x=type, ymin=range-sd_range, 
  #                  ymax=range+sd_range, color=group), 
  #              width=.1, position=position_dodge(width=0.5))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Climate change")+
  facet_grid(var~SSP, scale="free")+theme_bw()
p

species_env_change$var<-factor(species_env_change$var, levels=c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14"))
p<-ggplot(species_env_change)+
  geom_boxplot(aes(x=type, y=range, color=group), position=position_dodge(width=0.9))+
  scale_color_manual(values=color_groups)+
  labs(color="Group")+
  xlab("")+ylab("Climate change")+
  facet_grid(var~SSP, scale="free")+theme_bw()
p

ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/var_range_group.png", cuts_methods), width=12, height=9)


layer_island<-raster("../../Objects/Island/islands.tif")
layer_mountain<-raster("../../Objects/Island/mountain.tif")
mask<-raster("../../Raster/mask_100km.tif")

mask_p<-data.table(rasterToPoints(mask))

mask_p$island<-raster::extract(layer_island, mask_p[,c("x", "y")])
mask_p$mountain<-raster::extract(layer_mountain, mask_p[,c("x", "y")])
mask_p$type<-"plain"
mask_p[mountain==1]$type<-"mountain"
mask_p[island==1]$type<-"island"

map_colors<-c("mountain"=colors_red[7],
                        "island"=colors_blue[7],
              plain=colors_black[3])
p<-ggplot(mask_p)+geom_tile(aes(x=x, y=y, fill=type))+theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_rect(fill = map_background, color = NA), 
  panel.background = element_blank(), 
  legend.background = element_rect(fill = map_background, color = NA),
  panel.border = element_blank()
)+
  labs(fill="")+
  scale_fill_manual(values=map_colors)
ggsave(p, filename=sprintf("../../Figures/Mountain_island/%s/lands.png", cuts_methods), width=6, height=3)
