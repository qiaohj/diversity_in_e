library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
exposure=0
rda<-readRDS(sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", 0))
rda$exposure<-0
rda2<-readRDS(sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", 5))
rda2$exposure<-5
rda<-rbindlist(list(rda, rda2), fill=T)
rda<-rda[year==2100]
table(rda$year)
distribution_all<-readRDS("../../Objects/N_cell_init_distribution.rda")
distribution_all$sp<-gsub(" ", "_", distribution_all$sp)
fullset<-merge(rda, distribution_all, by.x="sp", by.y="sp")

quantile(fullset[N_CELL==0]$N)
table(fullset[N_CELL==0]$N)


group_disp_mammals<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
group_disp_birds<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
colnames(group_disp_birds)[5]<-"Scientific"
group_disp<-rbindlist(list(group_disp_mammals, group_disp_birds), fill=T)
group_disp$sp<-gsub(" ", "_", group_disp$Scientific)
fullset_mreged<-merge(fullset, group_disp, by.x="sp", by.y="sp", all.x=T)
fullset_mreged[is.na(estimated_disp)]
quantile(fullset_mreged[N_CELL==0]$estimated_disp)
fullset_mreged[N_CELL==0&N==707]
fullset_mreged[N_CELL==0&estimated_disp==max(fullset_mreged[N_CELL==0]$estimated_disp)]
group_disp[sp=="Abeillia_abeillei"]
cols<-c("Scientific", "group.x", "estimated_disp", "N", "GCM", "SSP", "M", "exposure", "N_CELL")
fullset_mreged_item<-fullset_mreged[, ..cols]
#library(stringr)
#fullset_mreged_item$exposure<-as.numeric(str_split_fixed(fullset_mreged$TYPE, "_", 5)[,3])
big_sp<-unique(fullset_mreged_item[N_CELL==0])
table(big_sp$group.x)

saveRDS(big_sp, "../../Objects/extincted_species_list.rda")
