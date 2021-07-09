library(data.table)
library(ggplot2)
library(dplyr)

g<-"Mammals"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
df_all<-NULL
for (g in c("Birds", "Mammals")){
  df<-readRDS(sprintf("../../Objects/Species_property/%s_property_compared.rda", g))
  df$group<-g
  df_all<-bind(df_all, df)
}
df_all[which(df_all$range_bio1_sd_max>df_all$bio1_max), "range_bio1_sd_max"]<-
  df_all[which(df_all$range_bio1_sd_max>df_all$bio1_max), "bio1_max"]

df_all[which(df_all$range_bio1_sd_min<df_all$bio1_min), "range_bio1_sd_min"]<-
  df_all[which(df_all$range_bio1_sd_min<df_all$bio1_min), "bio1_min"]

df_all[which(df_all$range_bio12_sd_max>df_all$bio12_max), "range_bio12_sd_max"]<-
  df_all[which(df_all$range_bio12_sd_max>df_all$bio12_max), "bio12_max"]

df_all[which(df_all$range_bio12_sd_min<df_all$bio12_min), "range_bio12_sd_min"]<-
  df_all[which(df_all$range_bio12_sd_min<df_all$bio12_min), "bio12_min"]
df_all$diff_t_max<-df_all$range_bio1_sd_max-df_all$bio1_max
df_all$diff_t_min<-df_all$range_bio1_sd_min-df_all$bio1_min
df_all$diff_bio12_max<-df_all$range_bio12_sd_max-df_all$bio12_max
df_all$diff_bio12_min<-df_all$range_bio12_sd_min-df_all$bio12_min
df_all$nb_bio1_sd

df_N_CELL<-df_all%>%select(sp, group, N_CELL)
df_N_CELL$N_CELL<-df_N_CELL$N_CELL/100
colnames(df_N_CELL)[3]<-"V"
df_N_CELL$TYPE="Distribution range"

df_nb_bio1_sd<-df_all%>%select(sp, group, nb_bio1_sd)
colnames(df_nb_bio1_sd)[3]<-"V"
df_nb_bio1_sd$TYPE="Range of temperature"

df_nb_bio12_sd<-df_all%>%select(sp, group, nb_bio12_sd)
colnames(df_nb_bio12_sd)[3]<-"V"
df_nb_bio12_sd$TYPE="Range of precipitation"
df_g<-bind_rows(bind_rows(df_N_CELL, df_nb_bio1_sd), df_nb_bio12_sd)


p<-ggplot(df_nb_bio1_sd)+
  geom_density(aes(x = V/100, color=group, fill = group), alpha=0.3)+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="", color="", x="Niche breadth in temperature")

p

saveRDS(p, "../../Figures/NB_hist_combined/nb_temp_hist.rda")

ggsave(p, filename="../../Figures/niche_property/nb_hist_temp.png")
ggsave(p, filename="../../Figures/niche_property/nb_hist_temp.pdf")

p<-ggplot(df_nb_bio12_sd)+
  geom_density(aes(x = V/100, color=group, fill = group), alpha=0.3)+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="", color="", x="Niche breadth in precipitation", y="")
p
saveRDS(p, "../../Figures/NB_hist_combined/nb_prec_hist.rda")
ggsave(p, filename="../../Figures/niche_property/nb_hist_prec.png")
ggsave(p, filename="../../Figures/niche_property/nb_hist_prec.pdf")

p<-ggplot(df_N_CELL)+
  geom_density(aes(x = V, color=group, fill = group), alpha=0.3)+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  scale_x_log10()+
  theme_bw()+
  labs(fill="", color="", x="Range size", y="")
saveRDS(p, "../../Figures/NB_hist_combined/nb_range_size_hist.rda")
p
ggsave(p, filename="../../Figures/niche_property/nb_hist_rangesize.png")
ggsave(p, filename="../../Figures/niche_property/nb_hist_rangesize.pdf")

  
df_N_CELL<-df_all%>%select(sp, group, N_CELL)
colnames(df_N_CELL)[3]<-"V"
df_N_CELL$TYPE="Distribution range"

df_nb_bio1_sd<-df_all%>%select(sp, group, nb_bio1_sd)
colnames(df_nb_bio1_sd)[3]<-"V"
df_nb_bio1_sd$V<-df_nb_bio1_sd$V/100
df_nb_bio1_sd$TYPE="Range of temperature"

df_nb_bio12_sd<-df_all%>%select(sp, group, nb_bio12_sd)
colnames(df_nb_bio12_sd)[3]<-"V"
df_nb_bio12_sd$V<-df_nb_bio12_sd$V/100
df_nb_bio12_sd$TYPE="Range of precipitation"
df_g<-bind_rows(bind_rows(df_N_CELL, df_nb_bio1_sd), df_nb_bio12_sd)

p<-ggplot(df_g)+
  geom_density(aes(x = V, y=..count.., color=group, fill = group), position = "identity", alpha=0.3)+
  scale_x_log10()+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="Clade", color="Clade", y="Number of species")+
  facet_wrap(~TYPE, nrow=3, scale="free")+
  theme(axis.title.x=element_blank())
p
ggsave(p, filename="../../Figures/niche_property/hist_all.png")
ggsave(p, filename="../../Figures/niche_property/hist_all.pdf")
