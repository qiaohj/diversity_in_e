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

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df_all$diff_t_max<-df_all$range_TEMP_sd_max-df_all$t_max_max
df_all$diff_t_min<-df_all$range_TEMP_sd_min-df_all$t_min_min
df_all$diff_pr_max<-df_all$range_PR_sd_max-df_all$pr_max
df_all$diff_pr_min<-df_all$range_PR_sd_min-df_all$pr_min
df_all$diff_t_range<-df_all$nb_TEMP_sd-df_all$range_real_temp
df_all$diff_pr_range<-df_all$nb_PR_sd-df_all$range_real_prec
df_all$diff_nb_range<-df_all$sd_nb_size-df_all$real_nb_size

df_all$TEMP_MAX_density <- get_density(df_all$range_TEMP_sd_max, df_all$t_max_max, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=range_TEMP_sd_max, y=t_max_max, color=TEMP_MAX_density))+
  xlim(0, 60)+
  ylim(0, 60)+
  geom_text(x=10, y=60, 
                label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                              cor_df$p_t_max, cor_df$mean_diff_t_max))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$TEMP_MAX_density, 0.2))+
  theme_bw()+
  labs(x="Species Tmax (excluding outliers)",
       y="Species Tmax (including outliers)",
       color="Density")+
  theme(legend.position = "none")
p
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_max.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_max.png")

df_all$TEMP_MIN_density <- get_density(df_all$range_TEMP_sd_min, df_all$diff_t_min, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=range_TEMP_sd_min, y=t_min_min, color=TEMP_MIN_density))+
  xlim(-80, 30)+
  ylim(-80, 30)+
  geom_text(x=-60, y=20, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                          cor_df$p_t_min, cor_df$mean_diff_t_min))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$TEMP_MIN_density, 0.2))+
  theme_bw()+
  labs(x="Species Tmin (excluding outliers)",
       y="Species Tmin (including outliers)",
       color="Density")+
  theme(legend.position = "none")
p
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_min.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_min.png")

df_all$PR_MAX_density <- get_density(df_all$range_PR_sd_max, df_all$pr_max, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=range_PR_sd_max, y=pr_max, color=PR_MAX_density))+
  xlim(0, 8000)+
  ylim(0, 8000)+
  geom_text(x=1500, y=7500, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_pr_max, cor_df$mean_diff_pr_max))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$PR_MAX_density, 0.2))+
  theme_bw()+
  labs(x="Species Pmax (excluding outliers)",
       y="Species Pmax (including outliers)",
       color="Density")+
  theme(legend.position = "none")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_max.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_max.png")

df_all$PR_MIN_density <- get_density(df_all$range_PR_sd_min, df_all$pr_min, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=range_PR_sd_min, y=pr_min, color=PR_MIN_density))+
  xlim(0, 6000)+
  ylim(0, 6000)+
  geom_text(x=1000, y=5000, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_pr_min, cor_df$mean_diff_pr_min))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$PR_MIN_density, 0.2))+
  theme_bw()+
  labs(x="Species Pmin (excluding outliers)",
       y="Species Pmin (including outliers)",
       color="Density")+
  theme(legend.position = "none")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_min.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/PR_min.png")

df_all$nb_TEMP_density <- get_density(df_all$nb_TEMP_sd, df_all$range_real_temp, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=nb_TEMP_sd, y=range_real_temp, color=nb_TEMP_density))+
  xlim(0, 200)+
  ylim(0, 150)+
  geom_text(x=30, y=150, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                          cor_df$p_range_t, cor_df$mean_diff_range_t))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$nb_TEMP_density, 0.2))+
  theme_bw()+
  labs(x="Species temperature breadth (excluding outliers)",
       y="Species temperature breadth (including outliers)",
       color="Density")+
  theme(legend.position = "none")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_breadth.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Temp_breadth.png")

df_all$nb_PR_density <- get_density(df_all$nb_PR_sd, df_all$range_real_prec, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=nb_PR_sd, y=range_real_prec, color=nb_PR_density))+
  xlim(0, 11500)+
  ylim(0, 8000)+
  geom_text(x=2000, y=8000, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3fmm", 
                          cor_df$p_range_p, cor_df$mean_diff_range_p))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$nb_PR_density, 0.2))+
  theme_bw()+
  labs(x="Species precipitation breadth (excluding outliers)",
       y="Species precipitation breadth (including outliers)",
       color="Density")+
  theme(legend.position = "none")
ggsave(p, filename="../../Figures_Full_species/niche_property/Prec_breadth.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Prec_breadth.png")

df_all$nb_size_density <- get_density(df_all$sd_nb_size, df_all$real_nb_size, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=sd_nb_size, y=real_nb_size, color=nb_size_density))+
  #xlim(0, 11500)+
  #ylim(0, 8000)+
  geom_text(x=2, y=22, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f", 
                          cor_df$p_nb, cor_df$mean_diff_nb))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$nb_size_density, 0.01))+
  theme_bw()+
  labs(x="Species niche area (excluding outliers)",
       y="Species niche area (including outliers)",
       color="Density")+
  theme(legend.position = "none")
p
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area.pdf")
ggsave(p, filename="../../Figures_Full_species/niche_property/Niche_area.png")



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

df_g$nb_size_density <- get_density(df_g$NB_1850, df_g$NB_1970, n = 100)


p<-ggplot(df_all)+geom_point(aes(x=sd_nb_size, y=real_nb_size, color=nb_size_density))+
  #xlim(0, 11500)+
  #ylim(0, 8000)+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1.5, alpha=0.5)+
  geom_text(x=3, y=22, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f", 
                          cor_df$p_nb, cor_df$mean_diff_nb))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$nb_size_density, 0.01))+
  theme_bw()+
  labs(x="Species niche area (excluding outliers)",
       y="Species niche area (including outliers)",
       color="Density")+
  theme(legend.position = "none")
p

p2<-ggplot(df_g)+geom_point(aes(x=NB_1850, y=NB_1970, color=nb_size_density))+
  #xlim(0, 11500)+
  #ylim(0, 8000)+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1.5, alpha=0.5)+
  geom_text(x=2.5, y=10, 
            label=sprintf("ρ=%.3f, mean(x-y)=%.3f", 
                          p_v, mean))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_g$nb_size_density, 0.01))+
  theme_bw()+
  labs(x="Species niche area 1850 - 2020",
       y="Species niche area 1970 - 2020",
       color="Density")+
  theme(legend.position = "none")
p2
ggsave(p2, filename="../../Figures_Full_species/niche_property/Niche_area_1850_1970.pdf")
ggsave(p2, filename="../../Figures_Full_species/niche_property/Niche_area_1850_1970.png")

pp<-ggarrange(p, p2, nrow = 1, ncol=2, labels=c("(a)", "(b)"))
pp
ggsave(pp, filename="../../Figures_Full_species/niche_property/Niche_area_combined.pdf",width=10, height=4)
ggsave(pp, filename="../../Figures_Full_species/niche_property/Niche_area_combined.png",width=10, height=4)
