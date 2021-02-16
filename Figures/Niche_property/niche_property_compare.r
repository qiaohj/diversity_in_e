df_all<-NULL
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  df<-readRDS(sprintf("../../Objects_Full_species/Species_property/%s_property_compared.rda", g))
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

ttt<-2
df_all<-df_all[N_CELL>ttt]


df_all$scaled_t_max_max<-scale(c(df_all$t_max_max, df_all$t_min_min, 
                                 df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$t_max_max)]
df_all$scaled_t_min_min<-scale(c(df_all$t_min_min, df_all$t_max_max, 
                                 df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$t_max_max)]
df_all$scaled_range_TEMP_sd_max<-scale(c(df_all$range_TEMP_sd_max, df_all$t_max_max, df_all$t_min_min, 
                                         df_all$range_TEMP_sd_min))[1:length(df_all$t_max_max)]
df_all$scaled_range_TEMP_sd_min<-scale(c(df_all$range_TEMP_sd_min, df_all$t_max_max, df_all$t_min_min, 
                                         df_all$range_TEMP_sd_max))[1:length(df_all$t_max_max)]

df_all$scaled_range_real_temp<-df_all$scaled_t_max_max-df_all$scaled_t_min_min
df_all$scaled_range_SD_temp<-df_all$scaled_range_TEMP_sd_max-df_all$scaled_range_TEMP_sd_min

df_all$scaled_pr_max<-scale(c(df_all$pr_max, df_all$pr_min, 
                              df_all$range_PR_sd_max, df_all$range_PR_sd_min))[1:length(df_all$pr_max)]
df_all$scaled_pr_min<-scale(c(df_all$pr_min, df_all$pr_max, 
                              df_all$range_PR_sd_max, df_all$range_PR_sd_min))[1:length(df_all$pr_max)]
df_all$scaled_range_PR_sd_max<-scale(c(df_all$range_PR_sd_max, df_all$pr_min, 
                                       df_all$pr_max, df_all$range_PR_sd_min))[1:length(df_all$pr_max)]
df_all$scaled_range_PR_sd_min<-scale(c(df_all$range_PR_sd_min, df_all$pr_min, df_all$pr_max, 
                                       df_all$range_PR_sd_max))[1:length(df_all$pr_max)]

plot(df_all$scaled_pr_max, df_all$scaled_range_PR_sd_max)
plot(df_all$scaled_pr_min, df_all$scaled_range_PR_sd_min)
plot(df_all$pr_max, df_all$range_PR_sd_max)
plot(df_all$pr_min, df_all$range_PR_sd_min)

df_all$scaled_range_real_prec<-df_all$scaled_pr_max-df_all$scaled_pr_min
df_all$scaled_range_SD_prec<-df_all$scaled_range_PR_sd_max-df_all$scaled_range_PR_sd_min

df_all$real_nb_size<-df_all$scaled_range_real_temp*df_all$scaled_range_real_prec
df_all$sd_nb_size<-df_all$scaled_range_SD_temp*df_all$scaled_range_SD_prec


p_v<-cor(df_all[target=="Niche_Models"]$sd_nb_size, df_all[target=="Niche_Models_1850_1970"]$sd_nb_size)
mean<-mean(df_all[target=="Niche_Models"]$sd_nb_size- df_all[target=="Niche_Models_1850_1970"]$sd_nb_size)

mean(df_all[target=="Niche_Models"]$range_TEMP_sd_max-df_all[target=="Niche_Models_1850_1970"]$range_TEMP_sd_max)
mean(df_all[target=="Niche_Models"]$range_TEMP_sd_min-df_all[target=="Niche_Models_1850_1970"]$range_TEMP_sd_min)
mean(df_all[target=="Niche_Models"]$range_PR_sd_max -df_all[target=="Niche_Models_1850_1970"]$range_PR_sd_max)
mean(df_all[target=="Niche_Models"]$range_PR_sd_min -df_all[target=="Niche_Models_1850_1970"]$range_PR_sd_min)

df_g<-data.frame(NB_1850=df_all[target=="Niche_Models"]$sd_nb_size,
                 NB_1970=df_all[target=="Niche_Models_1850_1970"]$sd_nb_size)
df_g$diff_nb_range<-df_g$NB_1850-df_g$NB_1970
p<-ggplot(df_g)+geom_point(aes(x=NB_1850, y=NB_1970, color=diff_nb_range))+
  #xlim(0, 11500)+
  #ylim(0, 8000)+
  geom_text(x=2, y=10, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f", 
                          p_v, mean))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species niche area 1850",
       y="Species niche area 1970",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area_1850_1970.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area_1850_1970.png")
