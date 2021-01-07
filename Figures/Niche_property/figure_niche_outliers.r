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
                   sd_diff_pr_min=sd(df_all$range_PR_sd_min-df_all$pr_min)
                   )

write.csv(cor_df, "../../Figures_Full_species/niche_property/cor.csv")
df_all$diff_t_max<-df_all$range_TEMP_sd_max-df_all$t_max_max
df_all$diff_t_min<-df_all$range_TEMP_sd_min-df_all$t_min_min
df_all$diff_pr_max<-df_all$range_PR_sd_max-df_all$pr_max
df_all$diff_pr_min<-df_all$range_PR_sd_min-df_all$pr_min

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
  xlim(-90, 30)+
  ylim(-90, 30)+
  geom_text(x=-80, y=20, 
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
  geom_text(x=1000, y=7500, 
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


