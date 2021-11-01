library(data.table)
library(ggplot2)
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
df_all[which(df_all$range_bio1_sd_min>df_all$bio1_min), "range_bio1_sd_min"]<-
  df_all[which(df_all$range_bio1_sd_min>df_all$bio1_min), "bio1_min"]

df_all[which(df_all$range_bio5_sd_max>df_all$bio5_max), "range_bio5_sd_max"]<-
  df_all[which(df_all$range_bio5_sd_max>df_all$bio5_max), "bio5_max"]
df_all[which(df_all$range_bio5_sd_min>df_all$bio5_min), "range_bio5_sd_min"]<-
  df_all[which(df_all$range_bio5_sd_min>df_all$bio5_min), "bio5_min"]

df_all[which(df_all$range_bio6_sd_max>df_all$bio6_max), "range_bio6_sd_max"]<-
  df_all[which(df_all$range_bio6_sd_max>df_all$bio6_max), "bio6_max"]
df_all[which(df_all$range_bio6_sd_min>df_all$bio6_min), "range_bio6_sd_min"]<-
  df_all[which(df_all$range_bio6_sd_min>df_all$bio6_min), "bio6_min"]

df_all[which(df_all$range_bio12_sd_max>df_all$bio12_max), "range_bio12_sd_max"]<-
  df_all[which(df_all$range_bio12_sd_max>df_all$bio12_max), "bio12_max"]
df_all[which(df_all$range_bio12_sd_min>df_all$bio12_min), "range_bio12_sd_min"]<-
  df_all[which(df_all$range_bio12_sd_min>df_all$bio12_min), "bio12_min"]

df_all[which(df_all$range_bio13_sd_max>df_all$bio13_max), "range_bio13_sd_max"]<-
  df_all[which(df_all$range_bio13_sd_max>df_all$bio13_max), "bio13_max"]
df_all[which(df_all$range_bio13_sd_min>df_all$bio13_min), "range_bio13_sd_min"]<-
  df_all[which(df_all$range_bio13_sd_min>df_all$bio13_min), "bio13_min"]

df_all[which(df_all$range_bio14_sd_max>df_all$bio14_max), "range_bio14_sd_max"]<-
  df_all[which(df_all$range_bio14_sd_max>df_all$bio14_max), "bio14_max"]
df_all[which(df_all$range_bio14_sd_min>df_all$bio14_min), "range_bio14_sd_min"]<-
  df_all[which(df_all$range_bio14_sd_min>df_all$bio14_min), "bio14_min"]

df_all$range_real_bio1<-df_all$bio1_max-df_all$bio1_min
df_all$range_real_bio5<-df_all$bio5_max-df_all$bio5_min
df_all$range_real_bio6<-df_all$bio6_max-df_all$bio6_min
df_all$range_real_bio12<-df_all$bio12_max-df_all$bio12_min
df_all$range_real_bio13<-df_all$bio13_max-df_all$bio13_min
df_all$range_real_bio14<-df_all$bio14_max-df_all$bio14_min
df_all<-df_all[(!is.infinite(bio1_max))&
                 (!is.infinite(bio5_max))&
                 (!is.infinite(bio6_max))&
                 (!is.infinite(bio12_max))&
                 (!is.infinite(bio13_max))&
                 (!is.infinite(bio14_max))&
                 (!is.infinite(bio1_min))&
                 (!is.infinite(bio5_min))&
                 (!is.infinite(bio6_min))&
                 (!is.infinite(bio12_min))&
                 (!is.infinite(bio13_min))&
                 (!is.infinite(bio14_min))]

df_all$scaled_bio1_max<-scale(c(df_all$bio1_max, df_all$bio1_min, 
                                 df_all$range_bio1_sd_max, 
                                df_all$range_bio1_sd_min))[1:length(df_all$bio1_max)]

df_all$scaled_bio1_min<-scale(c(df_all$bio1_min, df_all$bio1_max, 
                                 df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio1_max)]
df_all$scaled_range_bio1_sd_max<-scale(c(df_all$range_bio1_sd_max, df_all$bio1_max, df_all$bio1_min, 
                                  df_all$range_bio1_sd_min))[1:length(df_all$bio1_max)]
df_all$scaled_range_bio1_sd_min<-scale(c(df_all$range_bio1_sd_min, df_all$bio1_max, df_all$bio1_min, 
                                         df_all$range_bio1_sd_max))[1:length(df_all$bio1_max)]

df_all$scaled_range_real_bio1<-df_all$scaled_bio1_max-df_all$scaled_bio1_min
df_all$scaled_range_SD_bio1<-df_all$scaled_range_bio1_sd_max-df_all$scaled_range_bio1_sd_min
plot(df_all$scaled_range_real_bio1, df_all$scaled_range_SD_bio1)

df_all$scaled_bio5_max<-scale(c(df_all$bio5_max, df_all$bio5_min, 
                                df_all$range_bio5_sd_max, 
                                df_all$range_bio5_sd_min))[1:length(df_all$bio5_max)]
df_all$scaled_bio5_min<-scale(c(df_all$bio5_min, df_all$bio5_max, 
                                df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio5_max)]
df_all$scaled_range_bio5_sd_max<-scale(c(df_all$range_bio5_sd_max, df_all$bio5_max, df_all$bio5_min, 
                                         df_all$range_bio5_sd_min))[1:length(df_all$bio5_max)]
df_all$scaled_range_bio5_sd_min<-scale(c(df_all$range_bio5_sd_min, df_all$bio5_max, df_all$bio5_min, 
                                         df_all$range_bio5_sd_max))[1:length(df_all$bio5_max)]

df_all$scaled_range_real_bio5<-df_all$scaled_bio5_max-df_all$scaled_bio5_min
df_all$scaled_range_SD_bio5<-df_all$scaled_range_bio5_sd_max-df_all$scaled_range_bio5_sd_min
plot(df_all$scaled_range_real_bio5, df_all$scaled_range_SD_bio5)

df_all$scaled_bio6_max<-scale(c(df_all$bio6_max, df_all$bio6_min, 
                                df_all$range_bio6_sd_max, 
                                df_all$range_bio6_sd_min))[1:length(df_all$bio6_max)]
df_all$scaled_bio6_min<-scale(c(df_all$bio6_min, df_all$bio6_max, 
                                df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio6_max)]
df_all$scaled_range_bio6_sd_max<-scale(c(df_all$range_bio6_sd_max, df_all$bio6_max, df_all$bio6_min, 
                                         df_all$range_bio6_sd_min))[1:length(df_all$bio6_max)]
df_all$scaled_range_bio6_sd_min<-scale(c(df_all$range_bio6_sd_min, df_all$bio6_max, df_all$bio6_min, 
                                         df_all$range_bio6_sd_max))[1:length(df_all$bio6_max)]

df_all$scaled_range_real_bio6<-df_all$scaled_bio6_max-df_all$scaled_bio6_min
df_all$scaled_range_SD_bio6<-df_all$scaled_range_bio6_sd_max-df_all$scaled_range_bio6_sd_min
plot(df_all$scaled_range_real_bio6, df_all$scaled_range_SD_bio6)

df_all$scaled_bio12_max<-scale(c(df_all$bio12_max, df_all$bio12_min, 
                                df_all$range_bio12_sd_max, 
                                df_all$range_bio12_sd_min))[1:length(df_all$bio12_max)]
df_all$scaled_bio12_min<-scale(c(df_all$bio12_min, df_all$bio12_max, 
                                df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio12_max)]
df_all$scaled_range_bio12_sd_max<-scale(c(df_all$range_bio12_sd_max, df_all$bio12_max, df_all$bio12_min, 
                                         df_all$range_bio12_sd_min))[1:length(df_all$bio12_max)]
df_all$scaled_range_bio12_sd_min<-scale(c(df_all$range_bio12_sd_min, df_all$bio12_max, df_all$bio12_min, 
                                         df_all$range_bio12_sd_max))[1:length(df_all$bio12_max)]

df_all$scaled_range_real_bio12<-df_all$scaled_bio12_max-df_all$scaled_bio12_min
df_all$scaled_range_SD_bio12<-df_all$scaled_range_bio12_sd_max-df_all$scaled_range_bio12_sd_min
plot(df_all$scaled_range_real_bio12, df_all$scaled_range_SD_bio12)

df_all$scaled_bio13_max<-scale(c(df_all$bio13_max, df_all$bio13_min, 
                                df_all$range_bio13_sd_max, 
                                df_all$range_bio13_sd_min))[1:length(df_all$bio13_max)]
df_all$scaled_bio13_min<-scale(c(df_all$bio13_min, df_all$bio13_max, 
                                df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio13_max)]
df_all$scaled_range_bio13_sd_max<-scale(c(df_all$range_bio13_sd_max, df_all$bio13_max, df_all$bio13_min, 
                                         df_all$range_bio13_sd_min))[1:length(df_all$bio13_max)]
df_all$scaled_range_bio13_sd_min<-scale(c(df_all$range_bio13_sd_min, df_all$bio13_max, df_all$bio13_min, 
                                         df_all$range_bio13_sd_max))[1:length(df_all$bio13_max)]

df_all$scaled_range_real_bio13<-df_all$scaled_bio13_max-df_all$scaled_bio13_min
df_all$scaled_range_SD_bio13<-df_all$scaled_range_bio13_sd_max-df_all$scaled_range_bio13_sd_min
plot(df_all$scaled_range_real_bio13, df_all$scaled_range_SD_bio13)

df_all$scaled_bio14_max<-scale(c(df_all$bio14_max, df_all$bio14_min, 
                                df_all$range_bio14_sd_max, 
                                df_all$range_bio14_sd_min))[1:length(df_all$bio14_max)]
df_all$scaled_bio14_min<-scale(c(df_all$bio14_min, df_all$bio14_max, 
                                df_all$range_TEMP_sd_max, df_all$range_TEMP_sd_min))[1:length(df_all$bio14_max)]
df_all$scaled_range_bio14_sd_max<-scale(c(df_all$range_bio14_sd_max, df_all$bio14_max, df_all$bio14_min, 
                                         df_all$range_bio14_sd_min))[1:length(df_all$bio14_max)]
df_all$scaled_range_bio14_sd_min<-scale(c(df_all$range_bio14_sd_min, df_all$bio14_max, df_all$bio14_min, 
                                         df_all$range_bio14_sd_max))[1:length(df_all$bio14_max)]

df_all$scaled_range_real_bio14<-df_all$scaled_bio14_max-df_all$scaled_bio14_min
df_all$scaled_range_SD_bio14<-df_all$scaled_range_bio14_sd_max-df_all$scaled_range_bio14_sd_min
plot(df_all$scaled_range_real_bio14, df_all$scaled_range_SD_bio14)

df_all$real_nb_size<-df_all$scaled_range_real_bio1*df_all$scaled_range_real_bio12
df_all$sd_nb_size<-df_all$scaled_range_SD_bio1*df_all$scaled_range_SD_bio12

plot(df_all$sd_nb_size, df_all$real_nb_size)

cor_df<-data.frame(p_bio1_max=cor(df_all$bio1_max, df_all$range_bio1_sd_max, method="spearman"),
                   mean_diff_bio1_max=mean(df_all$range_bio1_sd_max-df_all$bio1_max),
                   sd_diff_bio1_max=sd(df_all$range_bio1_sd_max-df_all$bio1_max),
                   p_bio5_max=cor(df_all$bio5_max, df_all$range_bio5_sd_max, method="spearman"),
                   mean_diff_bio5_max=mean(df_all$range_bio5_sd_max-df_all$bio5_max),
                   sd_diff_bio5_max=sd(df_all$range_bio5_sd_max-df_all$bio5_max),
                   p_bio6_max=cor(df_all$bio6_max, df_all$range_bio6_sd_max, method="spearman"),
                   mean_diff_bio6_max=mean(df_all$range_bio6_sd_max-df_all$bio6_max),
                   sd_diff_bio6_max=sd(df_all$range_bio6_sd_max-df_all$bio6_max),
                   p_bio12_max=cor(df_all$bio12_max, df_all$range_bio12_sd_max, method="spearman"),
                   mean_diff_bio12_max=mean(df_all$range_bio12_sd_max-df_all$bio12_max),
                   sd_diff_bio12_max=sd(df_all$range_bio12_sd_max-df_all$bio12_max),
                   p_bio13_max=cor(df_all$bio13_max, df_all$range_bio13_sd_max, method="spearman"),
                   mean_diff_bio13_max=mean(df_all$range_bio13_sd_max-df_all$bio13_max),
                   sd_diff_bio13_max=sd(df_all$range_bio13_sd_max-df_all$bio13_max),
                   p_bio14_max=cor(df_all$bio14_max, df_all$range_bio14_sd_max, method="spearman"),
                   mean_diff_bio14_max=mean(df_all$range_bio14_sd_max-df_all$bio14_max),
                   sd_diff_bio14_max=sd(df_all$range_bio14_sd_max-df_all$bio14_max),
                   
                   p_bio1_min=cor(df_all$bio1_min, df_all$range_bio1_sd_min, method="spearman"),
                   mean_diff_bio1_min=mean(df_all$range_bio1_sd_min-df_all$bio1_min),
                   sd_diff_bio1_min=sd(df_all$range_bio1_sd_min-df_all$bio1_min),
                   p_bio5_min=cor(df_all$bio5_min, df_all$range_bio5_sd_min, method="spearman"),
                   mean_diff_bio5_min=mean(df_all$range_bio5_sd_min-df_all$bio5_min),
                   sd_diff_bio5_min=sd(df_all$range_bio5_sd_min-df_all$bio5_min),
                   p_bio6_min=cor(df_all$bio6_min, df_all$range_bio6_sd_min, method="spearman"),
                   mean_diff_bio6_min=mean(df_all$range_bio6_sd_min-df_all$bio6_min),
                   sd_diff_bio6_min=sd(df_all$range_bio6_sd_min-df_all$bio6_min),
                   p_bio12_min=cor(df_all$bio12_min, df_all$range_bio12_sd_min, method="spearman"),
                   mean_diff_bio12_min=mean(df_all$range_bio12_sd_min-df_all$bio12_min),
                   sd_diff_bio12_min=sd(df_all$range_bio12_sd_min-df_all$bio12_min),
                   p_bio13_min=cor(df_all$bio13_min, df_all$range_bio13_sd_min, method="spearman"),
                   mean_diff_bio13_min=mean(df_all$range_bio13_sd_min-df_all$bio13_min),
                   sd_diff_bio13_min=sd(df_all$range_bio13_sd_min-df_all$bio13_min),
                   p_bio14_min=cor(df_all$bio14_min, df_all$range_bio14_sd_min, method="spearman"),
                   mean_diff_bio14_min=mean(df_all$range_bio14_sd_min-df_all$bio14_min),
                   sd_diff_bio14_min=sd(df_all$range_bio14_sd_min-df_all$bio14_min),
                   
                   p_range_bio1=cor(df_all$range_real_bio1, df_all$nb_bio1_sd, method="spearman"),
                   mean_diff_range_bio1=mean(df_all$nb_bio1_sd-df_all$range_real_bio1),
                   sd_diff_range_bio1=sd(df_all$nb_bio1_sd-df_all$range_real_bio1),
                   
                   p_range_bio5=cor(df_all$range_real_bio5, df_all$nb_bio5_sd, method="spearman"),
                   mean_diff_range_bio5=mean(df_all$nb_bio5_sd-df_all$range_real_bio5),
                   sd_diff_range_bio5=sd(df_all$nb_bio5_sd-df_all$range_real_bio5),
                   
                   p_range_bio6=cor(df_all$range_real_bio6, df_all$nb_bio6_sd, method="spearman"),
                   mean_diff_range_bio6=mean(df_all$nb_bio6_sd-df_all$range_real_bio6),
                   sd_diff_range_bio6=sd(df_all$nb_bio6_sd-df_all$range_real_bio6),
                   
                   p_range_bio12=cor(df_all$range_real_bio12, df_all$nb_bio12_sd, method="spearman"),
                   mean_diff_range_bio12=mean(df_all$nb_bio12_sd-df_all$range_real_bio12),
                   sd_diff_range_bio12=sd(df_all$nb_bio12_sd-df_all$range_real_bio12),
                   
                   p_range_bio13=cor(df_all$range_real_bio13, df_all$nb_bio13_sd, method="spearman"),
                   mean_diff_range_bio13=mean(df_all$nb_bio13_sd-df_all$range_real_bio13),
                   sd_diff_range_bio13=sd(df_all$nb_bio13_sd-df_all$range_real_bio13),
                   
                   p_range_bio14=cor(df_all$range_real_bio14, df_all$nb_bio14_sd, method="spearman"),
                   mean_diff_range_bio14=mean(df_all$nb_bio14_sd-df_all$range_real_bio14),
                   sd_diff_range_bio14=sd(df_all$nb_bio14_sd-df_all$range_real_bio14),
                   
                   p_nb=cor(df_all$sd_nb_size, df_all$real_nb_size, method="spearman"),
                   mean_diff_nb=mean(df_all$sd_nb_size-df_all$real_nb_size),
                   sd_diff_nb=sd(df_all$sd_nb_size-df_all$real_nb_size)
                   
                   )

write.csv(cor_df, "../../Figures/niche_property/cor.csv")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df_all$diff_bio1_max<-df_all$range_bio1_sd_max-df_all$bio1_max
df_all$diff_bio1_min<-df_all$range_bio1_sd_min-df_all$bio1_min
df_all$diff_bio5_max<-df_all$range_bio5_sd_max-df_all$bio5_max
df_all$diff_bio5_min<-df_all$range_bio5_sd_min-df_all$bio5_min
df_all$diff_bio6_max<-df_all$range_bio6_sd_max-df_all$bio6_max
df_all$diff_bio6_min<-df_all$range_bio6_sd_min-df_all$bio6_min
df_all$diff_bio12_max<-df_all$range_bio12_sd_max-df_all$bio12_max
df_all$diff_bio12_min<-df_all$range_bio12_sd_min-df_all$bio12_min
df_all$diff_bio13_max<-df_all$range_bio13_sd_max-df_all$bio13_max
df_all$diff_bio13_min<-df_all$range_bio13_sd_min-df_all$bio13_min
df_all$diff_bio14_max<-df_all$range_bio14_sd_max-df_all$bio14_max
df_all$diff_bio14_min<-df_all$range_bio14_sd_min-df_all$bio14_min

df_all$diff_bio1_range<-df_all$nb_bio1_sd-df_all$range_real_bio1
df_all$diff_bio5_range<-df_all$nb_bio5_sd-df_all$range_real_bio5
df_all$diff_bio6_range<-df_all$nb_bio6_sd-df_all$range_real_bio6
df_all$diff_bio12_range<-df_all$nb_bio12_sd-df_all$range_real_bio12
df_all$diff_bio13_range<-df_all$nb_bio13_sd-df_all$range_real_bio13
df_all$diff_bio14_range<-df_all$nb_bio14_sd-df_all$range_real_bio14

df_all$diff_nb_range<-df_all$sd_nb_size-df_all$real_nb_size

df_all$bio1_max_density <- get_density(df_all$range_bio1_sd_max, df_all$bio1_max, n = 100)

p<-ggplot(df_all)+geom_point(aes(x=range_bio1_sd_max/100, y=bio1_max/100, color=bio1_max_density))+
  xlim(0, 40)+
  ylim(0, 40)+
  geom_text(x=10, y=60, 
                label=sprintf("ρ=%.3f, mean(x-y)=%.3f°C", 
                              cor_df$p_bio1_max, cor_df$mean_diff_bio1_max))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(df_all$bio1_max_density, 0.2))+
  theme_bw()+
  labs(x="Species bio1 (excluding outliers)",
       y="Species bio1 (including outliers)",
       color="Density")+
  theme(legend.position = "none")
p
ggsave(p, filename="../../Figures/niche_property/bio1_max.pdf")
ggsave(p, filename="../../Figures/niche_property/bio1_max.png")


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
ggsave(p, filename="../../Figures/niche_property/Niche_area.pdf")
ggsave(p, filename="../../Figures/niche_property/Niche_area.png")



p_v<-cor(df_all[target=="1850"]$sd_nb_size, df_all[target=="1970"]$sd_nb_size)
mean<-mean(df_all[target=="1850"]$sd_nb_size- df_all[target=="1970"]$sd_nb_size)


df_g<-data.frame(NB_1850=df_all[target=="1850"]$sd_nb_size,
                 NB_1970=df_all[target=="1970"]$sd_nb_size)
df_g$diff_nb_range<-df_g$NB_1850-df_g$NB_1970

df_g$nb_size_density <- get_density(df_g$NB_1850, df_g$NB_1970, n = 100)
df_all$nb_size_d

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
  geom_text(x=3, y=10, 
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
ggsave(p2, filename="../../Figures/niche_property/Niche_area_1850_1970.pdf")
ggsave(p2, filename="../../Figures/niche_property/Niche_area_1850_1970.png")
library(ggpubr)
pp<-ggarrange(p, p2, nrow = 1, ncol=2, labels=c("A", "B"))
pp
ggsave(pp, filename="../../Figures/niche_property/Niche_area_combined.pdf",width=10, height=4)
ggsave(pp, filename="../../Figures/niche_property/Niche_area_combined.png",width=10, height=4)
