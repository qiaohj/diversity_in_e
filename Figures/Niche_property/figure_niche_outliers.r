library(data.table)
library(ggplot2)
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

df_all$range_real_temp<-df_all$t_max_max-df_all$t_min_min
df_all$range_real_prec<-df_all$pr_max-df_all$pr_min
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
plot(df_all$scaled_range_real_temp, df_all$scaled_range_SD_temp)

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
plot(df_all$scaled_range_real_prec, df_all$scaled_range_SD_prec)

hist(df_all$scaled_range_real_prec)
df_all$real_nb_size<-df_all$scaled_range_real_temp*df_all$scaled_range_real_prec
df_all$sd_nb_size<-df_all$scaled_range_SD_temp*df_all$scaled_range_SD_prec

#df_all$real_nb_size<-df_all$range_real_prec*df_all$range_real_temp
#df_all$sd_nb_size<-df_all$nb_PR_sd*df_all$nb_TEMP_sd
ttt<-2
df_all<-df_all[N_CELL>ttt]
plot(df_all$sd_nb_size, df_all$real_nb_size)

cor_df<-data.frame(p_t_max=cor(df_all$t_max_max, df_all$range_TEMP_sd_max, method="spearman"),
                   mean_diff_t_max=mean(df_all$range_TEMP_sd_max-df_all$t_max_max),
                   sd_diff_t_max=sd(df_all$range_TEMP_sd_max-df_all$t_max_max),
                   p_t_min=cor(df_all$t_min_min, df_all$range_TEMP_sd_min, method="spearman"),
                   mean_diff_t_min=mean(df_all$range_TEMP_sd_min-df_all$t_min_min),
                   sd_diff_t_min=sd(df_all$range_TEMP_sd_min-df_all$t_min_min),
                   p_pr_max=cor(df_all$pr_max, df_all$range_PR_sd_max, method="spearman"),
                   mean_diff_pr_max=mean(df_all$range_PR_sd_max-df_all$pr_max),
                   sd_diff_pr_max=sd(df_all$range_PR_sd_max-df_all$pr_max),
                   p_pr_min=cor(df_all$pr_min, df_all$range_PR_sd_min, method="spearman"),
                   mean_diff_pr_min=mean(df_all$range_PR_sd_min-df_all$pr_min),
                   sd_diff_pr_min=sd(df_all$range_PR_sd_min-df_all$pr_min),
                   p_range_t=cor(df_all$range_real_temp, df_all$nb_TEMP_sd, method="spearman"),
                   mean_diff_range_t=mean(df_all$nb_TEMP_sd-df_all$range_real_temp),
                   sd_diff_range_t=sd(df_all$nb_TEMP_sd-df_all$range_real_temp),
                   p_range_p=cor(df_all$range_real_pr, df_all$nb_PR_sd, method="spearman"),
                   mean_diff_range_p=mean(df_all$nb_PR_sd-df_all$range_real_prec),
                   sd_diff_range_p=sd(df_all$nb_PR_sd-df_all$range_real_prec),
                   p_nb=cor(df_all$sd_nb_size, df_all$real_nb_size, method="spearman"),
                   mean_diff_nb=mean(df_all$sd_nb_size-df_all$real_nb_size),
                   sd_diff_nb=sd(df_all$sd_nb_size-df_all$real_nb_size)
                   
                   )

write.csv(cor_df, "../../Figures_Full_species/niche_property/cor.csv")
df_all$diff_t_max<-df_all$range_TEMP_sd_max-df_all$t_max_max
df_all$diff_t_min<-df_all$range_TEMP_sd_min-df_all$t_min_min
df_all$diff_pr_max<-df_all$range_PR_sd_max-df_all$pr_max
df_all$diff_pr_min<-df_all$range_PR_sd_min-df_all$pr_min
df_all$diff_t_range<-df_all$nb_TEMP_sd-df_all$range_real_temp
df_all$diff_pr_range<-df_all$nb_PR_sd-df_all$range_real_prec
df_all$diff_nb_range<-df_all$sd_nb_size-df_all$real_nb_size

p<-ggplot(df_all)+geom_point(aes(x=range_TEMP_sd_max, y=t_max_max, color=diff_t_max))+
  xlim(0, 60)+
  ylim(0, 60)+
  geom_text(x=10, y=60, 
                label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                              cor_df$p_t_max, cor_df$mean_diff_t_max))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species Tmax (excluding outliers)",
       y="Species Tmax (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_max.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_max.png")


p<-ggplot(df_all)+geom_point(aes(x=range_TEMP_sd_min, y=t_min_min, color=diff_t_min))+
  xlim(-80, 30)+
  ylim(-80, 30)+
  geom_text(x=-60, y=20, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                          cor_df$p_t_min, cor_df$mean_diff_t_min))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species Tmin (excluding outliers)",
       y="Species Tmin (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_min.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_min.png")

p<-ggplot(df_all)+geom_point(aes(x=range_PR_sd_max, y=pr_max, color=diff_pr_max))+
  xlim(0, 8000)+
  ylim(0, 8000)+
  geom_text(x=1500, y=7500, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_pr_max, cor_df$mean_diff_pr_max))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species Pmax (excluding outliers)",
       y="Species Pmax (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_max.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_max.png")

p<-ggplot(df_all)+geom_point(aes(x=range_PR_sd_min, y=pr_min, color=diff_pr_min))+
  xlim(0, 6000)+
  ylim(0, 6000)+
  geom_text(x=1000, y=5000, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_pr_min, cor_df$mean_diff_pr_min))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species Pmin (excluding outliers)",
       y="Species Pmin (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_min.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_min.png")

p<-ggplot(df_all)+geom_point(aes(x=nb_TEMP_sd, y=range_real_temp, color=diff_t_range))+
  xlim(0, 200)+
  ylim(0, 150)+
  geom_text(x=30, y=150, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                          cor_df$p_range_t, cor_df$mean_diff_range_t))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species temperature breadth (excluding outliers)",
       y="Species temperature breadth (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_breadth.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_breadth.png")

p<-ggplot(df_all)+geom_point(aes(x=nb_PR_sd, y=range_real_prec, color=diff_pr_range))+
  xlim(0, 11500)+
  ylim(0, 8000)+
  geom_text(x=2000, y=8000, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_range_p, cor_df$mean_diff_range_p))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species precipitation breadth (excluding outliers)",
       y="Species precipitation breadth (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Prec_breadth.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Prec_breadth.png")

p<-ggplot(df_all)+geom_point(aes(x=sd_nb_size, y=real_nb_size, color=diff_nb_range))+
  #xlim(0, 11500)+
  #ylim(0, 8000)+
  geom_text(x=2, y=28, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f", 
                          cor_df$p_nb, cor_df$mean_diff_nb))+
  scale_color_gradient(low="blue", high="red")+
  theme_bw()+
  labs(x="Species niche area (excluding outliers)",
       y="Species niche area (including outliers)",
       color="Differ")
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area.png")
