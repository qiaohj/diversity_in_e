library(dplyr)
library(data.table)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/colors.r")
source("commonFuns/functions.r")
ttt=2
if (F){
  df_all<-NULL
  for (threshold in c(1, 5)){
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      df<-readRDS(sprintf("../../Objects_Full_species/dispersal_trait/dispersal_trait_exposure_%d_%s.rda", threshold, g))
      df_all<-bind(df_all, df)
    }
  }
  saveRDS(df_all, "../../Objects_Full_species/dispersal_trait/all.rda")
}
df_all<-readRDS("../../Objects_Full_species/dispersal_trait/all.rda")

df_all$status<-ifelse(df_all$extinct_year==2101, "extant", "extinct")
df_all$exposure<-ifelse(df_all$threshold==1, " no exposure", "5-year exposure")
df_all$da<-ifelse(df_all$dispersal==0, "no dispersal", "with dispersal")

df_se<-df_all%>%dplyr::group_by(type, SSP, da, exposure, status, group)%>%
  dplyr::summarise(alt=mean(mean_alt),
                   sd_alt=sd(mean_alt),
                   CI_alt=CI(mean_alt)[1]-CI(mean_alt)[2],
                   y=mean(mean_abs_y),
                   sd_y=sd(mean_abs_y),
                   CI_y=CI(mean_abs_y)[1]-CI(mean_abs_y)[2])

df_se$type<-factor(df_se$type, levels=c("start", "end"))
p1<-ggplot(df_se%>%dplyr::filter(da=="with dispersal"))+
  geom_point(aes(x=type, y=alt, color=group), size=2)+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Elevation (m)", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


legend<-g_legend(p1)

df_se$y<-df_se$y/1000
df_se$CI_y<-df_se$CI_y/1000
p2<-ggplot(df_se%>%dplyr::filter(da=="with dispersal"))+
  geom_point(aes(x=type, y=y, color=group))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Latitude (km)", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank())

pp<-ggarrange(p1, p2, nrow=2, ncol=1, common.legend=T, legend.grob=legend, legend="right")
pp2<-annotate_figure(pp,
                     bottom = text_grob("Time spots", size = 10)
)

write.csv(df_se, sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d.csv", "ALL", ttt), row.names=F)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d.png", "ALL", ttt), width=10, height=8)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d.pdf", "ALL", ttt), width=10, height=8)





#Figure XXX

df_se<-df_all%>%dplyr::group_by(type, SSP, da, exposure, group)%>%
  dplyr::summarise(alt=mean(mean_alt),
                   sd_alt=sd(mean_alt),
                   CI_alt=CI(mean_alt)[1]-CI(mean_alt)[2],
                   y=mean(mean_abs_y),
                   sd_y=sd(mean_abs_y),
                   CI_y=CI(mean_abs_y)[1]-CI(mean_abs_y)[2])

df_se$type<-factor(df_se$type, levels=c("start", "end"))
df_se$y<-df_se$y/1000
df_se$CI_y<-df_se$CI_y/1000

p1<-ggplot(df_se%>%dplyr::filter((SSP=="SSP119")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=alt, color=group), size=2)+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  #facet_wrap(~SSP, scale="free", strip.position="right", nrow=2)+
  labs(x="Time spots", y="Elevation (m)", color="Group", linetype="Exposure")+
  theme_bw()+
  scale_x_discrete(expand=c(0.1,0), drop=FALSE)+
  guides(color = guide_legend(nrow = 2, byrow = T)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.direction = "vertical", legend.box = "horizontal")
legend<-g_legend(p1)


p2<-ggplot(df_se%>%dplyr::filter((SSP=="SSP119")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=y, color=group))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  #facet_wrap(~SSP, scale="free", strip.position="right", nrow=2)+
  labs(x="Time spots", y="Latitude (km)", color="Group", linetype="Exposure")+
  theme_bw()+
  scale_x_discrete(expand=c(0.1,0), drop=FALSE)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0, "lines"))

pp<-ggarrange(p1, p2, nrow=1, ncol=2, common.legend=T, legend.grob=legend, legend="top")
pp2=pp
#pp2<-annotate_figure(pp,
#                     bottom = text_grob("Begin to end", size = 10)
#)

#write.csv(raw_path_final_se_g, sprintf("../../Figures_Full_species/Top_Figure_all/gradient_%s_ttt_%d.csv", "ALL", ttt), row.names=F)
ggsave(pp2, filename=
         sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_combined.png", "ALL", ttt), 
       width=3.5, height=4)
ggsave(pp2, filename=
         sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_combined.pdf", "ALL", ttt), 
       width=4, height=6)



#Figure YYYYY

df_se<-df_all%>%dplyr::group_by(type, SSP, da, status, exposure)%>%
  dplyr::summarise(alt=mean(mean_alt),
                   sd_alt=sd(mean_alt),
                   CI_alt=CI(mean_alt)[1]-CI(mean_alt)[2],
                   y=mean(mean_abs_y),
                   sd_y=sd(mean_abs_y),
                   CI_y=CI(mean_abs_y)[1]-CI(mean_abs_y)[2])

df_se$type<-factor(df_se$type, levels=c("start", "end"))
df_se$y<-df_se$y/1000
df_se$CI_y<-df_se$CI_y/1000

p1<-ggplot(df_se%>%dplyr::filter((da=="with dispersal")))+
  geom_point(aes(x=type, y=alt, color=status), size=2)+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=status), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=status, linetype=factor(SSP)))+
  scale_color_manual(values=color_survive)+
  facet_wrap(~exposure, scale="free", strip.position="right", nrow=2)+
  labs(x="Time spots", y="Elevation (m)", color="Status", linetype="SSP")+
  theme_bw()+
  guides(color = guide_legend(nrow = 3, byrow = T)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.direction = "vertical", legend.box = "horizontal")
legend<-g_legend(p1)


p2<-ggplot(df_se%>%dplyr::filter((da=="with dispersal")))+
  geom_point(aes(x=type, y=y, color=status))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=status), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=status, linetype=factor(SSP)))+
  scale_color_manual(values=color_survive)+
  facet_wrap(~exposure, scale="free", strip.position="right", nrow=2)+
  labs(x="Time spots", y="Latitude (km)", color="Status", linetype="SSP")+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"))

pp<-ggarrange(p1, p2, nrow=2, ncol=1, common.legend=T, legend.grob=legend, legend="top")
pp2=pp
#pp2<-annotate_figure(pp,
#                     bottom = text_grob("Begin to end", size = 10)
#)

#write.csv(raw_path_final_se_g, sprintf("../../Figures_Full_species/Top_Figure_all/gradient_%s_ttt_%d.csv", "ALL", ttt), row.names=F)
ggsave(pp2, filename=
         sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_combined_all_sp.png", "ALL", ttt), 
       width=4, height=8)
ggsave(pp2, filename=
         sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_combined_all_sp.pdf", "ALL", ttt), 
       width=4, height=8)



#Separated by tropic and n/s temperate
tropic_ll<-data.table(lon=0, lat=c(23.43654, -23.43654))
mask_ll<-raster("../../Raster/Continent.tif")
mask_eck4<-raster("../../Raster/Continent_ect4.tif")
p_tropic_all<-SpatialPointsDataFrame(tropic_ll, tropic_ll, proj4string=crs(mask_ll))
p_tropic_all_eck4<-spTransform(p_tropic_all, crs(mask_eck4))
lat_threshold<-p_tropic_all_eck4@coords[,2]

df_all_start<-df_all%>%dplyr::filter(type=="start")

df_all_start$is_tropic<-"tropics"
df_all_start[(df_all_start$mean_y>lat_threshold[1]), "is_tropic"]<-"north temperate"
df_all_start[(df_all_start$mean_y<lat_threshold[2]), "is_tropic"]<-"south temperate"
df_all_start<-unique(df_all_start[, c("sp", "group", "is_tropic")])
df_all_with_tropic<-left_join(df_all, df_all_start, by=c("sp", "group"))

ggplot(df_all_start[sample(nrow(df_all_start), 1000),])+
  geom_point(aes(x=mean_y, y=mean_y, color=factor(is_tropic)))

df_se<-df_all_with_tropic%>%dplyr::group_by(type, SSP, da, status, exposure, is_tropic, group)%>%
  dplyr::summarise(alt=mean(mean_alt),
                   sd_alt=sd(mean_alt),
                   CI_alt=CI(mean_alt)[1]-CI(mean_alt)[2],
                   y=mean(mean_abs_y),
                   sd_y=sd(mean_abs_y),
                   CI_y=CI(mean_abs_y)[1]-CI(mean_abs_y)[2])

df_se$type<-factor(df_se$type, levels=c("start", "end"))
df_se$y<-df_se$y/1000
df_se$CI_y<-df_se$CI_y/1000

p1<-ggplot(df_se%>%filter((is_tropic=="tropics")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=y, color=group))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Latitude (km) in tropics", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
legend<-g_legend(p1)

p2<-ggplot(df_se%>%filter((is_tropic=="north temperate")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=y, color=group))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Latitude (km) in north temperate", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3<-ggplot(df_se%>%filter((is_tropic=="south temperate")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=y, color=group))+
  geom_errorbar(aes(x=type, y=y, ymin=y-CI_y, ymax=y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Latitude (km) in south temperate", color="Group", linetype="Exposure")+
  theme_bw()+
  theme()


pp<-ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend=T, legend.grob=legend, legend="right")
pp2<-annotate_figure(pp,
                     bottom = text_grob("Time spots", size = 10)
)

#write.csv(raw_path_final_se_g, sprintf("../../Figures_Full_species/Top_Figure_all/gradient_%s_ttt_%d.csv", "ALL", ttt), row.names=F)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_by_tropics_ns.png", 
                             "ALL", ttt), width=10, height=10)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_by_tropics_ns.pdf", 
                             "ALL", ttt), width=12, height=8)


p1<-ggplot(df_se%>%filter((is_tropic=="tropics")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=alt, color=group))+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Elevation (m) in tropics", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
legend<-g_legend(p1)

p2<-ggplot(df_se%>%filter((is_tropic=="north temperate")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=alt, color=group))+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Elevation (m) in north temperate", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3<-ggplot(df_se%>%filter((is_tropic=="south temperate")&(da=="with dispersal")))+
  geom_point(aes(x=type, y=alt, color=group))+
  geom_errorbar(aes(x=type, y=alt, ymin=alt-CI_alt, ymax=alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(type), y=alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(status~SSP, scale="free")+
  labs(x="Time spots", y="Elevation (m) in south temperate", color="Group", linetype="Exposure")+
  theme_bw()+
  theme()


pp<-ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend=T, legend.grob=legend, legend="right")
pp2<-annotate_figure(pp,
                     bottom = text_grob("Time spots", size = 10)
)

#write.csv(raw_path_final_se_g, sprintf("../../Figures_Full_species/Top_Figure_all/gradient_%s_ttt_%d.csv", "ALL", ttt), row.names=F)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_by_tropics_ns_alt.png", 
                             "ALL", ttt), width=10, height=10)
ggsave(pp2, filename=sprintf("../../Figures_Full_species/dispersal_traits/gradient_%s_ttt_%d_by_tropics_ns_alt.pdf", 
                             "ALL", ttt), width=12, height=8)

