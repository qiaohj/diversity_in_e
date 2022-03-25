library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/colors.r")
df<-readRDS("../../Objects/mig_bird_log.rda")
stat<-df[YEAR==2021]
stat<-stat[, .(N=.N), by=list(dispersal, exposure, esm, ssp, res)]
df_se<-df[YEAR==2100]
df_se$Extincted_mig<-df_se$N_mig==0
df_se$Extincted_all<-df_se$N_all==0
df_sex<-df_se[, .(N=.N), by=list(dispersal, exposure, esm, ssp)]
max(df_sex$N)
df_se<-df_se[, .(N=.N), by=list(dispersal, exposure, esm, ssp, Extincted_mig, Extincted_all)]
df_se<-df_se[, .(N=mean(N), SD=sd(N)),
             by=list(dispersal, exposure, ssp, Extincted_mig, Extincted_all)]
df_se$All<-1050
df_se$per<-df_se$N/df_se$All
write.csv(df_se, "../../Figures/Area_change_mig/extinct_mig.csv", row.names = F)

df_max_mig<-df[, .(max_mig=max(N_mig), max_all=max(N_all)), by=list(dispersal, exposure, esm, ssp, sp)]
df_max_mig<-merge(df, df_max_mig, by=c("dispersal", "exposure", "esm", "ssp", "sp"))
df_2100<-df_max_mig[YEAR==2100]
colnames(df_2100)[c(7:8)]<-c("N_mig_2100", "N_all_2100")
df_2020<-df_max_mig[YEAR==2021]
colnames(df_2020)[c(7:8)]<-c("N_mig_2020", "N_all_2020")
df_compare<-merge(df_2020, df_2100, by=c("dispersal", "exposure", "esm", "ssp", "sp", "res", "max_mig", "max_all"), all=T)
df_compare[is.na(N_mig_2100)]$N_mig_2100<-0
df_compare[is.na(N_all_2100)]$N_all_2100<-0
df_compare$change_mig<-df_compare$N_mig_2100 - df_compare$N_mig_2020
df_compare$change_all<-df_compare$N_all_2100 - df_compare$N_all_2020
df_compare$change_mig_per<-1-df_compare$N_mig_2100/df_compare$max_mig
df_compare$change_all_per<-1-df_compare$N_all_2100/df_compare$max_all

hist(df_compare$change_mig_per)
hist(df_compare$change_all_per)

quantile_mig<-quantile(df_compare$change_mig_per, c(0.25, 0.75))
iqr_mig<-quantile_mig[2]-quantile_mig[1]
iqr_range_mig<-c(quantile_mig[1]-1.5*iqr_mig,
                 quantile_mig[2]+1.5*iqr_mig)

quantile_all<-quantile(df_compare$change_all_per, c(0.25, 0.75))
iqr_all<-quantile_all[2]-quantile_all[1]
iqr_range_all<-c(quantile_all[1]-1.5*iqr_all,
                 quantile_all[2]+1.5*iqr_all)

df_compare_filter<-df_compare[between(change_all_per, iqr_range_all[1], iqr_range_all[2])&
                         between(change_mig_per, iqr_range_mig[1], iqr_range_mig[2])]
range(df_compare$change_mig_per)
hist(df_compare_filter$change_mig_per)
df_compare_se<-df_compare_filter[, .(change_mig_per=mean(change_mig_per), change_all_per=mean(change_all_per),
                              sd_change_mig_per=sd(change_mig_per), sd_change_all_per=sd(change_all_per)),
                          by=list(dispersal, exposure, ssp)]

write.csv(df_compare_se, "../../Figures/Area_change_mig/area_change_compare.csv", row.names = F)


df_se<-df[YEAR==2100]
df_se$Extincted_mig<-df_se$N_mig==0
df_se$Extincted_all<-df_se$N_all==0
df_sex<-df_se[, .(N=.N), by=list(dispersal, exposure, esm, ssp,Extincted_mig)]
df_sex$per<-df_sex$N/1050
df_sexx<-df_sex[,.(per=mean(per*100),
                   sd_per=sd(per*100)), 
                by=list(dispersal, exposure, ssp, Extincted_mig)]

df_sey<-df_se[, .(N=.N), by=list(dispersal, exposure, esm, ssp,Extincted_all)]
df_sey$per<-df_sey$N/1050
df_sey<-df_sey[,.(per=mean(per*100),
                   sd_per=sd(per*100)), 
                by=list(dispersal, exposure, ssp, Extincted_all)]

write.csv(df_sexx, "../../Figures/Area_change_mig/extinction_per.csv", row.names = F)

df_sexx$dispersal<-ifelse(df_sexx$dispersal==0, "no dispersal", "with dispersal")
df_sexx$exposure<-ifelse(df_sexx$exposure==0, " no climate resilience", "climate resilience")
p<-ggplot(df_sexx[Extincted_mig==T], aes())+
  geom_bar(stat="identity", position="stack", 
           aes(y=per, x=ssp, fill=ssp))+
  geom_errorbar(position=position_dodge(0.1), width=0.1,
                aes(ymin=per-sd_per, 
                    ymax=per+sd_per, x=ssp)) +
  
  xlab("SSP scenario")+
  theme_bw()+
  facet_grid(exposure~dispersal)+
  scale_fill_manual(values=color_ssp)+
  ylab("%s of species that lose all suitable breeding habitat")
p

ggsave(p, filename="../../Figures/Area_change_mig/extinction_per.png", width=6, height=4)
