library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/colors.r")
source("commonFuns/functions.r")
df<-readRDS("../../Objects/mean_env_year.rda")
df[which(df$VAR!="pr"),]$V<-df[which(df$VAR!="pr"),]$V*10-273.16
df_se_2020<-df%>%dplyr::filter(Y<=2020)%>%
  dplyr::group_by(Y, GCM, VAR)%>%
  dplyr::summarise(annul_prec=mean(V),
                   annul_max_temp=mean(V),
                   annual_min_temp=mean(V))
df_se_2020$SSP<-"SSP_ALL"
df_se_2100<-df%>%dplyr::filter(Y>2020)%>%
  dplyr::group_by(Y, GCM, SSP, VAR)%>%
  dplyr::summarise(annul_prec=sum(V),
                   annul_max_temp=mean(V),
                   annual_min_temp=mean(V))
df_se<-bind_rows(df_se_2020, df_se_2100)

df_se_se<-df_se%>%dplyr::group_by(Y, VAR, SSP)%>%
  dplyr::summarise(max_v=max(annul_prec),
                   mean_v=mean(annul_prec),
                   min_v=min(annul_prec))



when_extinct_1<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", 0))
names(when_extinct_1)[1]<-"Group"
when_extinct_1$label<-ifelse((when_extinct_1$dispersal==0), 
                             " no dispersal, no climate resilience",
                             "with dispersal, no climate resilience")
when_extinct_1$exposure<-" no climate resilience"
when_extinct_1$da<-ifelse((when_extinct_1$dispersal==0), 
                          "no dispersal",
                          "with dispersal")

when_extinct_5<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/when_extinct_final.rda", 5))
names(when_extinct_5)[1]<-"Group"
when_extinct_5$label<-ifelse((when_extinct_5$dispersal==0), 
                             " no dispersal, climate resilience",
                             "with dispersal, climate resilience")
when_extinct_5$exposure<-"climate resilience"
when_extinct_5$da<-ifelse((when_extinct_5$dispersal==0), 
                          "no dispersal",
                          "with dispersal")



when_extinct<-bind_rows(when_extinct_1, when_extinct_5)

when_extinct_se<-when_extinct%>%dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(n_sp=sum(mean_n_sp))
when_extinct_se$label<-paste(when_extinct_se$exposure, ", ", when_extinct_se$da, sep="")

dataset1<-df_se_se%>%filter((VAR=="tasmax")&(between(Y, 2021, 2100)))
#dataset2<-when_extinct_se%>%
#  dplyr:: filter(label %in% c(" no climate resilience, no dispersal", "climate resilience, with dispersal"))
dataset2<-when_extinct_se
max_y1<-max(dataset1$max_v)
min_y1<-min(dataset1$min_v)
max_y2<-max(dataset2$n_sp)
min_y2<-min(dataset2$n_sp)
scale<-(max_y2-min_y2)/(max_y1-min_y1)
intercept<-min_y1

dataset2$n_sp_scale<-dataset2$n_sp/scale+intercept
range(dataset2$n_sp_scale)

extinct_propotion<-read.csv("../../Figures/when_where_extinction_all/proportion_by_year.csv",
                            stringsAsFactors = F)
extinct_propotion<-extinct_propotion%>%dplyr::filter(Group=="ALL")
extinct_propotion_with_sp<-inner_join(extinct_propotion, dataset2, 
                              by=c("SSP", "exposure", "da", "extinct_year"))
labels<-extinct_propotion_with_sp%>%dplyr::filter(SSP=="SSP585")
#extinct_propotion_with_sp<-extinct_propotion_with_sp%>%filter(extinct_threshold>0.2)
dataset2%>%dplyr::filter((SSP=="SSP585")&(exposure=="climate resilience")&(extinct_year==2095))
p2<-ggplot()+
  geom_ribbon(data=dataset1, aes(x=Y, ymin=min_v, ymax=max_v, fill=SSP), alpha=0.2)+
  #geom_line(data=dataset1, aes(x=Y, y=mean_v, color=SSP))+
  geom_line(data=dataset2, aes(x=extinct_year, y=n_sp/scale+intercept, 
                               color=SSP, linetype=label),
            size=1)+
  geom_segment(data=extinct_propotion_with_sp, 
               aes(x = extinct_year, y=min_y1, 
                   xend = extinct_year, yend=n_sp_scale, 
                   color=SSP), linetype=2)+
  geom_segment(data=extinct_propotion_with_sp, 
               aes(x = extinct_year, y=n_sp_scale, 
                   xend = 2020, yend=n_sp_scale, 
                   color=SSP), linetype=2)+
  scale_color_manual(values=color_ssp)+
  scale_fill_manual(values=color_ssp)+
  scale_x_continuous(breaks = c(seq(2020, 2100, by=10)))+
  scale_y_continuous(
    name = "Annual maximum temperature",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~(.-intercept)*scale, 
                        name="Average number of extinctions across ESMs",
                        breaks=waiver(),
                        labels=waiver())
  )+
  geom_text(data=labels%>%dplyr::filter(label==" no climate resilience, no dispersal"), aes(x=2020, y=n_sp_scale, 
                             label=sprintf(">%d%%", extinct_threshold*100)),
            hjust=0.1, vjust=-0.8)+
  geom_text(data=extinct_propotion_with_sp, 
            aes(x=extinct_year-1.1, y=n_sp_scale+0.12, label=extinct_year), size=3)+
  
  xlab("Year")+
  scale_linetype_manual(values=linetype_label)+
  #scale_linetype_manual(values=linetype_label)+
  labs(linetype="", color="", fill="")+
  theme_bw()+
  #theme(legend.position = "bottom")
  theme(legend.position = "bottom", legend.box="vertical")
p2
ggsave(p2, filename="../../Figures/Combined_Figure/Figure_env_year.pdf", width=10, height=6)
ggsave(p2, filename="../../Figures/Combined_Figure/Figure_env_year.png", width=10, height=6)

p2<-ggplot()+
  geom_ribbon(data=dataset1, aes(x=Y, ymin=min_v, ymax=max_v, fill=SSP), alpha=0.2)+
  #geom_line(data=dataset1, aes(x=Y, y=mean_v, color=SSP))+
  geom_line(data=dataset2%>%dplyr::filter(label %in% c(" no climate resilience, no dispersal", "climate resilience, with dispersal")), 
            aes(x=extinct_year, y=n_sp/scale+intercept, color=SSP, linetype=label), size=1)+
  geom_segment(data=extinct_propotion_with_sp%>%dplyr::filter(label %in% c(" no climate resilience, no dispersal", "climate resilience, with dispersal")), 
               aes(x = extinct_year, y=min_y1, 
                   xend = extinct_year, yend=n_sp_scale, 
                   color=SSP), linetype=2)+
  geom_segment(data=extinct_propotion_with_sp%>%dplyr::filter(label %in% c(" no climate resilience, no dispersal", "climate resilience, with dispersal")), 
               aes(x = extinct_year, y=n_sp_scale, 
                   xend = 2020, yend=n_sp_scale, 
                   color=SSP), linetype=2)+
  scale_color_manual(values=color_ssp)+
  scale_fill_manual(values=color_ssp)+
  scale_x_continuous(breaks = c(seq(2020, 2100, by=10)))+
  scale_y_continuous(
    name = "Annual maximum temperature",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~(.-intercept)*scale, 
                        name="Average number of extinctions across ESMs",
                        breaks=waiver(),
                        labels=waiver())
  )+
  geom_text(data=labels%>%dplyr::filter(label==" no climate resilience, no dispersal"), aes(x=2020, y=n_sp_scale, 
                                                                                 label=sprintf(">%d%%", extinct_threshold*100)),
            hjust=0.1, vjust=-0.8)+
  geom_text(data=extinct_propotion_with_sp%>%dplyr::filter(label %in% c(" no climate resilience, no dispersal", "climate resilience, with dispersal")),
            aes(x=extinct_year-1.1, y=n_sp_scale+0.12, label=extinct_year), size=3)+
  
  xlab("Year")+
  #scale_linetype_manual(values=linetype_label)+
  labs(linetype="", color="", fill="")+
  theme_bw()+
  theme(legend.position = "bottom")
  #theme(legend.position = "bottom", legend.box="vertical")
p2
ggsave(p2, filename="../../Figures/Combined_Figure/Figure_env_year_few.pdf", width=12, height=6)
ggsave(p2, filename="../../Figures/Combined_Figure/Figure_env_year_few.png", width=12, height=6)

write.csv(extinct_propotion_with_sp, "../../Figures/Combined_Figure/extinct_propotion_with_sp.csv")
