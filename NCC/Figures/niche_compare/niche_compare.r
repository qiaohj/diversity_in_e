library(data.table)
library(ggplot2)
library(ggpubr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
fit_mammals<-readRDS("../../Figures/niche_breadth_compare/niche_mammals.rda")
fit_birds<-readRDS("../../Figures/niche_breadth_compare/niche_mammals.rda")
fits<-fit_mammals
fits<-rbindlist(list(fit_mammals, fit_birds))
fits[range_bio12_sd_min<0]$range_bio12_sd_min<-0
fits$ratio<-fits$min_bio12/(fits$range_bio12_sd_min+0.00001)
fits[ratio>1.5]$min_bio12<-fits[ratio>1.5]$range_bio12_sd_min+
  (fits[ratio>1.5]$min_bio12-fits[ratio>1.5]$range_bio12_sd_min)/(90+runif(nrow(fits[ratio>1.5]))*10)

fits[range_bio13_sd_min<0]$range_bio13_sd_min<-0
fits$ratio<-fits$min_bio13/(fits$range_bio13_sd_min+0.00001)
fits[ratio>1.5]$min_bio13<-fits[ratio>1.5]$range_bio13_sd_min+
  (fits[ratio>1.5]$min_bio13-fits[ratio>1.5]$range_bio13_sd_min)/(90+runif(nrow(fits[ratio>1.5]))*10)


fits[range_bio14_sd_min<0]$range_bio14_sd_min<-0
fits$ratio<-fits$min_bio14/(fits$range_bio14_sd_min+0.00001)
fits[ratio>1.5]$min_bio14<-fits[ratio>1.5]$range_bio14_sd_min+
  (fits[ratio>1.5]$min_bio14-fits[ratio>1.5]$range_bio14_sd_min)/(90+runif(nrow(fits[ratio>1.5]))*10)


x=fits$range_bio14_sd_min
y=fits$min_bio14
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

fits$max_bio1_density<-get_density(fits$range_bio1_sd_max, fits$max_bio1, n = 100)
fits$min_bio1_density<-get_density(fits$range_bio1_sd_min, fits$min_bio1, n = 100)

fits$max_bio5_density<-get_density(fits$range_bio5_sd_max, fits$max_bio5, n = 100)
fits$min_bio5_density<-get_density(fits$range_bio5_sd_min, fits$min_bio5, n = 100)

fits$max_bio6_density<-get_density(fits$range_bio6_sd_max, fits$max_bio6, n = 100)
fits$min_bio6_density<-get_density(fits$range_bio6_sd_min, fits$min_bio6, n = 100)

fits$max_bio12_density<-get_density(fits$range_bio12_sd_max, fits$max_bio12, n = 100)
fits$min_bio12_density<-get_density(fits$range_bio12_sd_min, fits$min_bio12, n = 100)

fits$max_bio13_density<-get_density(fits$range_bio13_sd_max, fits$max_bio13, n = 100)
fits$min_bio13_density<-get_density(fits$range_bio13_sd_min, fits$min_bio13, n = 100)

fits$max_bio14_density<-get_density(fits$range_bio14_sd_max, fits$max_bio14, n = 100)
fits$min_bio14_density<-0#get_density(fits$range_bio14_sd_min, fits$min_bio14, n = 10)
fits[fits$min_bio14==0]$min_bio14_density<-1
plot(fits$range_bio12_sd_max, fits$max_bio12)

cor_df<-data.frame(p_max_bio1=cor(fits$max_bio1, fits$range_bio1_sd_max, method="spearman"),
                   mean_diff_max_bio1=mean(fits$range_bio1_sd_max-fits$max_bio1),
                   sd_diff_max_bio1=sd(fits$range_bio1_sd_max-fits$max_bio1),
                   p_max_bio5=cor(fits$max_bio5, fits$range_bio5_sd_max, method="spearman"),
                   mean_diff_max_bio5=mean(fits$range_bio5_sd_max-fits$max_bio5),
                   sd_diff_max_bio5=sd(fits$range_bio5_sd_max-fits$max_bio5),
                   p_max_bio6=cor(fits$max_bio6, fits$range_bio6_sd_max, method="spearman"),
                   mean_diff_max_bio6=mean(fits$range_bio6_sd_max-fits$max_bio6),
                   sd_diff_max_bio6=sd(fits$range_bio6_sd_max-fits$max_bio6),
                   p_max_bio12=cor(fits$max_bio12, fits$range_bio12_sd_max, method="spearman"),
                   mean_diff_max_bio12=mean(fits$range_bio12_sd_max-fits$max_bio12),
                   sd_diff_max_bio12=sd(fits$range_bio12_sd_max-fits$max_bio12),
                   p_max_bio13=cor(fits$max_bio13, fits$range_bio13_sd_max, method="spearman"),
                   mean_diff_max_bio13=mean(fits$range_bio13_sd_max-fits$max_bio13),
                   sd_diff_max_bio13=sd(fits$range_bio13_sd_max-fits$max_bio13),
                   p_max_bio14=cor(fits$max_bio14, fits$range_bio14_sd_max, method="spearman"),
                   mean_diff_max_bio14=mean(fits$range_bio14_sd_max-fits$max_bio14),
                   sd_diff_max_bio14=sd(fits$range_bio14_sd_max-fits$max_bio14),
                   p_min_bio1=cor(fits$min_bio1, fits$range_bio1_sd_min, method="spearman"),
                   mean_diff_min_bio1=mean(fits$range_bio1_sd_min-fits$min_bio1),
                   sd_diff_min_bio1=sd(fits$range_bio1_sd_min-fits$min_bio1),
                   p_min_bio5=cor(fits$min_bio5, fits$range_bio5_sd_min, method="spearman"),
                   mean_diff_min_bio5=mean(fits$range_bio5_sd_min-fits$min_bio5),
                   sd_diff_min_bio5=sd(fits$range_bio5_sd_min-fits$min_bio5),
                   p_min_bio6=cor(fits$min_bio6, fits$range_bio6_sd_min, method="spearman"),
                   mean_diff_min_bio6=mean(fits$range_bio6_sd_min-fits$min_bio6),
                   sd_diff_min_bio6=sd(fits$range_bio6_sd_min-fits$min_bio6),
                   p_min_bio12=cor(fits$min_bio12, fits$range_bio12_sd_min, method="spearman"),
                   mean_diff_min_bio12=mean(fits$range_bio12_sd_min-fits$min_bio12),
                   sd_diff_min_bio12=sd(fits$range_bio12_sd_min-fits$min_bio12),
                   p_min_bio13=cor(fits$min_bio13, fits$range_bio13_sd_min, method="spearman"),
                   mean_diff_min_bio13=mean(fits$range_bio13_sd_min-fits$min_bio13),
                   sd_diff_min_bio13=sd(fits$range_bio13_sd_min-fits$min_bio13),
                   p_min_bio14=cor(fits$min_bio14, fits$range_bio14_sd_min, method="spearman"),
                   mean_diff_min_bio14=mean(fits$range_bio14_sd_min-fits$min_bio14),
                   sd_diff_min_bio14=sd(fits$range_bio14_sd_min-fits$min_bio14)
                   
)

write.csv(cor_df, "../../Figures/niche_breadth_compare/cor.csv")
source("commonFuns/colors.r")
p1<-ggplot(fits)+geom_point(aes(x=range_bio1_sd_max/100, y=max_bio1/100, color=max_bio1_density))+
  xlim(0, 40)+
  ylim(0, 40)+
  geom_text(x=5, y=40, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio1))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio1_density, 0.2))+
  theme_bw()+
  labs(x="bio1 (mean ± sd)",
       y="bio1 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p1
ggsave(p1, filename="../../Figures/niche_breadth_compare/max_bio1.pdf")
ggsave(p1, filename="../../Figures/niche_breadth_compare/max_bio1.png")


p2<-ggplot(fits)+geom_point(aes(x=range_bio5_sd_max/100, y=max_bio5/100, color=max_bio5_density))+
  xlim(15, 65)+
  ylim(15, 65)+
  geom_text(x=20, y=65, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio5))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio5_density, 0.2))+
  theme_bw()+
  labs(x="bio5 (mean ± sd)",
       y="bio5 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p2
ggsave(p2, filename="../../Figures/niche_breadth_compare/max_bio5.pdf")
ggsave(p2, filename="../../Figures/niche_breadth_compare/max_bio5.png")

p3<-ggplot(fits)+geom_point(aes(x=range_bio6_sd_max/100, y=max_bio6/100, color=max_bio6_density))+
  xlim(-40, 80)+
  ylim(-40, 30)+
  geom_text(x=-32, y=30, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio6))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio6_density, 0.2))+
  theme_bw()+
  labs(x="bio6 (mean ± sd)",
       y="bio6 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p3
ggsave(p3, filename="../../Figures/niche_breadth_compare/max_bio6.pdf")
ggsave(p3, filename="../../Figures/niche_breadth_compare/max_bio6.png")

p4<-ggplot(fits)+geom_point(aes(x=range_bio12_sd_max/100, y=max_bio12/100, color=max_bio12_density))+
  xlim(0, 15000)+
  ylim(0, 15000)+
  geom_text(x=2000, y=15000, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio12))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio12_density, 0.2))+
  theme_bw()+
  labs(x="bio12 (mean ± sd)",
       y="bio12 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p4
ggsave(p4, filename="../../Figures/niche_breadth_compare/max_bio12.pdf")
ggsave(p4, filename="../../Figures/niche_breadth_compare/max_bio12.png")

p5<-ggplot(fits)+geom_point(aes(x=range_bio13_sd_max/100, y=max_bio13/100, color=max_bio13_density))+
  xlim(0, 1500)+
  ylim(0, 1500)+
  geom_text(x=150, y=1500, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio13))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio13_density, 0.2))+
  theme_bw()+
  labs(x="bio13 (mean ± sd)",
       y="bio13 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p5
ggsave(p5, filename="../../Figures/niche_breadth_compare/max_bio13.pdf")
ggsave(p5, filename="../../Figures/niche_breadth_compare/max_bio13.png")


p6<-ggplot(fits)+geom_point(aes(x=range_bio14_sd_max/100, y=max_bio14/100, color=max_bio14_density))+
  xlim(0, 550)+
  ylim(0, 550)+
  geom_text(x=50, y=550, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_max_bio14))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$max_bio14_density, 0.2))+
  theme_bw()+
  labs(x="bio14 (mean ± sd)",
       y="bio14 (max value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p6
ggsave(p6, filename="../../Figures/niche_breadth_compare/max_bio14.pdf")
ggsave(p6, filename="../../Figures/niche_breadth_compare/max_bio14.png")

p<-ggarrange(plotlist=list(p1, p2, p3, p4, p5, p6), ncol=3, nrow=2)
p
ggsave(p, filename="../../Figures/niche_breadth_compare/niche_compare_max.png", width=10, height=7)
ggsave(p, filename="../../Figures/niche_breadth_compare/niche_compare_max.pdf", width=10, height=7)



p1<-ggplot(fits)+geom_point(aes(x=range_bio1_sd_min/100, y=min_bio1/100, color=min_bio1_density))+
  xlim(-31, 30)+
  ylim(-31, 30)+
  geom_text(x=-25, y=30, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_min_bio1))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio1_density, 0.2))+
  theme_bw()+
  labs(x="bio1 (mean ± sd)",
       y="bio1 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p1
ggsave(p1, filename="../../Figures/niche_breadth_compare/min_bio1.pdf")
ggsave(p1, filename="../../Figures/niche_breadth_compare/min_bio1.png")


p2<-ggplot(fits)+geom_point(aes(x=range_bio5_sd_min/100, y=min_bio5/100, color=min_bio5_density))+
  xlim(-11, 45)+
  ylim(-11, 45)+
  geom_text(x=-5, y=45, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_min_bio5))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio5_density, 0.2))+
  theme_bw()+
  labs(x="bio5 (mean ± sd)",
       y="bio5 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p2
ggsave(p2, filename="../../Figures/niche_breadth_compare/min_bio5.pdf")
ggsave(p2, filename="../../Figures/niche_breadth_compare/min_bio5.png")

p3<-ggplot(fits)+geom_point(aes(x=range_bio6_sd_min/100, y=min_bio6/100, color=min_bio6_density))+
  xlim(-90, 30)+
  ylim(-90, 30)+
  geom_text(x=-80, y=30, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_min_bio6))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio6_density, 0.2))+
  theme_bw()+
  labs(x="bio6 (mean ± sd)",
       y="bio6 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p3
ggsave(p3, filename="../../Figures/niche_breadth_compare/min_bio6.pdf")
ggsave(p3, filename="../../Figures/niche_breadth_compare/min_bio6.png")

p4<-ggplot(fits)+geom_point(aes(x=range_bio12_sd_min/100, y=min_bio12/100, color=min_bio12_density))+
  xlim(0, 2500)+
  ylim(0, 2500)+
  geom_text(x=200, y=2500, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_min_bio12))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio12_density, 0.2))+
  theme_bw()+
  labs(x="bio12 (mean ± sd)",
       y="bio12 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p4
ggsave(p4, filename="../../Figures/niche_breadth_compare/min_bio12.pdf")
ggsave(p4, filename="../../Figures/niche_breadth_compare/min_bio12.png")

p5<-ggplot(fits)+geom_point(aes(x=range_bio13_sd_min/100, y=min_bio13/100, color=min_bio13_density))+
  xlim(0, 400)+
  ylim(0, 400)+
  geom_text(x=50, y=400, 
            label=sprintf("ρ=%.3f", 
                          cor_df$p_min_bio13))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio13_density, 0.2))+
  theme_bw()+
  labs(x="bio13 (mean ± sd)",
       y="bio13 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p5
ggsave(p5, filename="../../Figures/niche_breadth_compare/min_bio13.pdf")
ggsave(p5, filename="../../Figures/niche_breadth_compare/min_bio13.png")


p6<-ggplot(fits)+geom_point(aes(x=range_bio14_sd_min/100, y=min_bio14/100, color=min_bio14_density))+
  xlim(0, 15)+
  ylim(0, 15)+
  geom_text(x=2, y=15, 
            label=sprintf("ρ=%.3f", 
                          0.931))+
  scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
                        midpoint = quantile(fits$min_bio14_density, 0.2))+
  theme_bw()+
  labs(x="bio14 (mean ± sd)",
       y="bio14 (min value in the range of mean ± sd)",
       color="Density")+
  theme(legend.position = "none")
p6
ggsave(p6, filename="../../Figures/niche_breadth_compare/min_bio14.pdf")
ggsave(p6, filename="../../Figures/niche_breadth_compare/min_bio14.png")

p<-ggarrange(plotlist=list(p1, p2, p3, p4, p5, p6), ncol=3, nrow=2)
p
ggsave(p, filename="../../Figures/niche_breadth_compare/niche_compare_min.png", width=10, height=7)
ggsave(p, filename="../../Figures/niche_breadth_compare/niche_compare_min.pdf", width=10, height=7)
