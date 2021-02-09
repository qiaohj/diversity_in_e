library(data.table)
library(ggplot2)
library(dplyr)

g<-"Amphibians"
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
df_all<-NULL
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df<-readRDS(sprintf("../../Objects_Full_species/Species_property/%s_property.rda", g))
  df$group<-g
  df_all<-bind(df_all, df)
}
df_all[which(df_all$range_TEMP_sd_max>df_all$t_max_max), "range_TEMP_sd_max"]<-
  df_all[which(df_all$range_TEMP_sd_max>df_all$t_max_max), "t_max_max"]

df_all[which(df_all$range_TEMP_sd_min<df_all$t_min_min), "range_TEMP_sd_min"]<-
  df_all[which(df_all$range_TEMP_sd_min<df_all$t_min_min), "t_min_min"]

df_all[which(df_all$range_PR_sd_max>df_all$pr_max), "range_PR_sd_max"]<-
  df_all[which(df_all$range_PR_sd_max>df_all$pr_max), "pr_max"]

df_all[which(df_all$range_PR_sd_min<df_all$pr_min), "range_PR_sd_min"]<-
  df_all[which(df_all$range_PR_sd_min<df_all$pr_min), "pr_min"]
df_all$diff_t_max<-df_all$range_TEMP_sd_max-df_all$t_max_max
df_all$diff_t_min<-df_all$range_TEMP_sd_min-df_all$t_min_min
df_all$diff_pr_max<-df_all$range_PR_sd_max-df_all$pr_max
df_all$diff_pr_min<-df_all$range_PR_sd_min-df_all$pr_min
df_all$nb_TEMP_sd

p<-ggplot(df_all)+
  geom_density(aes(x = nb_TEMP_sd, color=group, fill = group), alpha=0.3)+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="Group", color="Group", x="Temperature", y="Density")
ggsave(p, filename="../../Figures_Full_species/niche_property/nb_hist_temp.png")
ggsave(p, filename="../../Figures_Full_species/niche_property/nb_hist_temp.pdf")

p<-ggplot(df_all)+
  geom_density(aes(x = nb_PR_sd, color=group, fill = group), alpha=0.3)+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="Group", color="Group", x="Precipitation", y="Density")
ggsave(p, filename="../../Figures_Full_species/niche_property/nb_hist_prec.png")
ggsave(p, filename="../../Figures_Full_species/niche_property/nb_hist_prec.pdf")

  
ttt<-2
df_all<-df_all%>%dplyr::filter(N_CELL>ttt)
df_N_CELL<-df_all%>%select(sp, group, N_CELL)
colnames(df_N_CELL)[3]<-"V"
df_N_CELL$TYPE="Distribution range"

df_nb_TEMP_sd<-df_all%>%select(sp, group, nb_TEMP_sd)
colnames(df_nb_TEMP_sd)[3]<-"V"
df_nb_TEMP_sd$TYPE="Range of temperature"

df_nb_PR_sd<-df_all%>%select(sp, group, nb_PR_sd)
colnames(df_nb_PR_sd)[3]<-"V"
df_nb_PR_sd$TYPE="Range of precipitation"
df_g<-bind_rows(bind_rows(df_N_CELL, df_nb_TEMP_sd), df_nb_PR_sd)

p<-ggplot(df_g)+
  geom_density(aes(x = V, y=..count.., color=group, fill = group), position = "identity", alpha=0.3)+
  scale_x_log10()+
  scale_fill_manual(values=color_groups)+
  scale_color_manual(values=color_groups)+
  theme_bw()+
  labs(fill="Group", color="Group", y="Number of species")+
  facet_wrap(~TYPE, nrow=3, scale="free")+
  theme(axis.title.x=element_blank())
p
ggsave(p, filename="../../Figures_Full_species/niche_property/hist_all.png")
ggsave(p, filename="../../Figures_Full_species/niche_property/hist_all.pdf")