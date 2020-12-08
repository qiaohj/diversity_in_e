library(raster)
library(dplyr)
library(concaveman)
library(sf)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
alt<-raster("../../Raster/ALT/alt_eck4")


no_na<-!is.na(values(mask))

source("commonFuns/functions.r")
source("commonFuns/colors.r")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

if (F){
  
  j=9
  threshold=1
  raw_path_final_all<-NULL
  for (SSP_i in SSPs){
    for (threshold in c(1, 5)){
      for (GCM_i in GCMs){
        print(paste(SSP_i, GCM_i, threshold))
        raw_path<-readRDS(sprintf("../../Figures/Top_Figure_%d/raw_path_%s_%s.rda", threshold, GCM_i, SSP_i))
        raw_path<-data.frame(raw_path)
        raw_path$alt<-raster::extract(alt, raw_path[, c("gravity_x", "gravity_y")])
        saveRDS(raw_path, sprintf("../../Figures/Top_Figure_%d/raw_path_%s_%s_with_alt.rda", threshold, GCM_i, SSP_i))
        
        raw_path_min_max_year<-raw_path%>%
          dplyr::group_by(group, sp, survive, continent_i)%>%
          dplyr::summarise(MAX_YEAR=max(YEAR),
                           MIN_YEAR=min(YEAR))
        
        raw_path_2<-inner_join(raw_path_min_max_year, raw_path, 
                               by=c("group", "sp", "survive", "continent_i", "MAX_YEAR"="YEAR"))
        colnames(raw_path_2)[7:9]<-c("end_x", "end_y", "end_alt")
        
        raw_path_1<-inner_join(raw_path_min_max_year, raw_path, 
                               by=c("group", "sp", "survive", "continent_i", "MIN_YEAR"="YEAR"))
        colnames(raw_path_1)[7:9]<-c("start_x", "start_y", "start_alt")
        
        raw_path_final<-inner_join(raw_path_1, raw_path_2, 
                                   by=c("group", "sp", "survive", "continent_i", "MAX_YEAR", "MIN_YEAR"))
        
        raw_path_final<-raw_path_final%>%dplyr::filter(!is.na(start_alt)&
                                                         !is.na(end_alt))
        
        raw_path_final$start_agent<-raw_path_final$start_alt+abs(raw_path_final$start_y)/1000
        raw_path_final$end_agent<-raw_path_final$end_alt+abs(raw_path_final$end_y)/1000
        
        raw_path_final$alt_differ<-raw_path_final$end_alt-raw_path_final$start_agent
        raw_path_final$y_differ<-abs(raw_path_final$end_y)-abs(raw_path_final$start_agent)
        raw_path_final$agent_differ<-raw_path_final$end_agent-raw_path_final$start_agent
        
        raw_path_final$GCM<-GCM_i
        raw_path_final$SSP<-SSP_i
        raw_path_final$threshold<-threshold
        raw_path_final_all<-bind_dplyr(raw_path_final_all, raw_path_final)
      }
    }
  }
  saveRDS(raw_path_final_all, "../../Figures/Top_Figure_all/Data/raw_path_final_all.rda")
}

raw_path_final_all<-readRDS("../../Figures/Top_Figure_all/Data/raw_path_final_all.rda")


raw_path_final_se<-raw_path_final_all%>%dplyr::group_by(group, survive, SSP, threshold)%>%
  dplyr::summarise(mean_start_alt=mean(start_alt, na.rm=T),
                   sd_start_alt=sd(start_alt, na.rm=T),
                   CI_start_alt=CI(start_alt)[1]-CI(start_alt)[2],
                   mean_end_alt=mean(end_alt, na.rm=T),
                   sd_end_alt=sd(end_alt, na.rm=T),
                   CI_end_alt=CI(end_alt)[1]-CI(end_alt)[2],
                   mean_start_y=mean(abs(start_y), na.rm=T),
                   sd_start_y=sd(abs(start_y), na.rm=T),
                   CI_start_y=CI(start_y)[1]-CI(start_y)[2],
                   mean_end_y=mean(abs(end_y), na.rm=T),
                   sd_end_y=sd(abs(end_y), na.rm=T),
                   CI_end_y=CI(end_y)[1]-CI(end_y)[2],
                   mean_start_agent=mean(start_agent, na.rm=T),
                   sd_start_agent=sd(start_agent, na.rm=T),
                   CI_start_agent=CI(start_agent)[1]-CI(start_agent)[2],
                   mean_end_agent=mean(end_agent, na.rm=T),
                   sd_end_agent=sd(end_agent, na.rm=T),
                   CI_end_agent=CI(end_agent)[1]-CI(end_agent)[2])

raw_path_final_se_1<-raw_path_final_se[, c("group", "survive", "SSP", "threshold",
                                           "mean_start_alt", "sd_start_alt", "CI_start_alt",
                                           "mean_start_y", "sd_start_y", "CI_start_y",
                                           "mean_start_agent", "sd_start_agent", "CI_start_agent")]
colnames(raw_path_final_se_1)<-c("group", "survive","SSP", "threshold",
                                 "mean_alt", "sd_alt", "CI_alt",
                                 "mean_y", "sd_y", "CI_y",
                                 "mean_agent", "sd_agent", "CI_agent")
raw_path_final_se_1$year="Start"

raw_path_final_se_2<-raw_path_final_se[, c("group", "survive","SSP", "threshold",
                                           "mean_end_alt", "sd_end_alt", "CI_end_alt",
                                           "mean_end_y", "sd_end_y", "CI_end_y",
                                           "mean_end_agent", "sd_end_agent", "CI_end_agent")]
colnames(raw_path_final_se_2)<-c("group", "survive","SSP", "threshold",
                                 "mean_alt", "sd_alt", "CI_alt",
                                 "mean_y", "sd_y", "CI_y",
                                 "mean_agent", "sd_agent", "CI_agent")
raw_path_final_se_2$year="End"

raw_path_final_se_g<-rbind(raw_path_final_se_1, raw_path_final_se_2)
raw_path_final_se_g$year<-factor(raw_path_final_se_g$year, levels = c("Start", "End"))
raw_path_final_se_g$exposure<-" no exposure"
raw_path_final_se_g[which(raw_path_final_se_g$threshold==5),]$exposure<-"5-year exposure"
#raw_path_final_se_g$exposure<-factor(raw_path_final_se_g$exposure, levels=c("no exposure", "5-year exposure"))

p1<-ggplot(raw_path_final_se_g)+
  geom_point(aes(x=year, y=mean_alt, color=group), size=2)+
  geom_errorbar(aes(x=year, y=mean_alt, ymin=mean_alt-CI_alt, ymax=mean_alt+CI_alt, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(year), y=mean_alt, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(survive~SSP, scale="free")+
  labs(x="Time spots", y="Elevation (m)", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
legend<-g_legend(p1)

raw_path_final_se_g$mean_y<-raw_path_final_se_g$mean_y/1000
raw_path_final_se_g$CI_y<-raw_path_final_se_g$CI_y/1000
p2<-ggplot(raw_path_final_se_g)+
  geom_point(aes(x=year, y=mean_y, color=group))+
  geom_errorbar(aes(x=year, y=mean_y, ymin=mean_y-CI_y, ymax=mean_y+CI_y, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(year), y=mean_y, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(survive~SSP, scale="free")+
  labs(x="Time spots", y="Latitude (km)", color="Group", linetype="Exposure")+
  theme_bw()+
  theme(axis.title.x=element_blank())

pp<-ggarrange(p1, p2, nrow=2, ncol=1, common.legend=T, legend.grob=legend, legend="right")
pp2<-annotate_figure(pp,
                     bottom = text_grob("Time spots", size = 10)
)
ggsave(pp2, filename=sprintf("../../Figures/Top_Figure_all/gradient_%s.png", "ALL"), width=10, height=8)
ggsave(pp2, filename=sprintf("../../Figures/Top_Figure_all/gradient_%s.pdf", "ALL"), width=10, height=8)

p3<-ggplot(raw_path_final_se_g)+
  geom_point(aes(x=year, y=mean_agent, color=group))+
  geom_errorbar(aes(x=year, y=mean_agent, ymin=mean_agent-CI_agent, ymax=mean_agent+CI_agent, color=group), 
                width = 0.1, position = "dodge2")+
  geom_line(aes(x=as.numeric(year), y=mean_agent, color=group, linetype=factor(exposure)))+
  scale_color_manual(values=color_groups)+
  facet_grid(survive~SSP, scale="free")+
  labs(x="Time spots", y="Elevation + Latitude", color="Group", linetype="Exposure")+
  theme_bw()


raw_path_continent_se<-raw_path_final_all%>%dplyr::group_by(group, survive, continent_i, SSP, threshold)%>%
  dplyr::summarise(mean_start_alt=mean(start_alt, na.rm=T),
                   sd_start_alt=sd(start_alt, na.rm=T),
                   CI_start_alt=CI(start_alt)[1]-CI(start_alt)[2],
                   mean_end_alt=mean(end_alt, na.rm=T),
                   sd_end_alt=sd(end_alt, na.rm=T),
                   CI_end_alt=CI(end_alt)[1]-CI(end_alt)[2],
                   mean_start_y=mean(abs(start_y), na.rm=T),
                   sd_start_y=sd(abs(start_y), na.rm=T),
                   CI_start_y=CI(start_y)[1]-CI(start_y)[2],
                   mean_end_y=mean(abs(end_y), na.rm=T),
                   sd_end_y=sd(abs(end_y), na.rm=T),
                   CI_end_y=CI(end_y)[1]-CI(end_y)[2],
                   mean_start_agent=mean(start_agent, na.rm=T),
                   sd_start_agent=sd(start_agent, na.rm=T),
                   CI_start_agent=CI(start_agent)[1]-CI(start_agent)[2],
                   mean_end_agent=mean(end_agent, na.rm=T),
                   sd_end_agent=sd(end_agent, na.rm=T),
                   CI_end_agent=CI(end_agent)[1]-CI(end_agent)[2])

raw_path_continent_se$continent_label<-""
raw_path_continent_se[which(raw_path_continent_se$continent_i==1), "continent_label"]<-"Africa"
raw_path_continent_se[which(raw_path_continent_se$continent_i==2), "continent_label"]<-"Euroasia"
raw_path_continent_se[which(raw_path_continent_se$continent_i==3), "continent_label"]<-"Australia"
raw_path_continent_se[which(raw_path_continent_se$continent_i==5), "continent_label"]<-"North America"
raw_path_continent_se[which(raw_path_continent_se$continent_i==6), "continent_label"]<-"South America"

raw_path_continent_se_1<-raw_path_continent_se[, c("group", "survive", "SSP", "threshold", "continent_i", "continent_label",
                                           "mean_start_alt", "sd_start_alt", "CI_start_alt",
                                           "mean_start_y", "sd_start_y", "CI_start_y",
                                           "mean_start_agent", "sd_start_agent", "CI_start_agent")]

colnames(raw_path_continent_se_1)<-c("group", "survive","SSP", "threshold", "continent_i", "continent_label",
                                 "mean_alt", "sd_alt", "CI_alt",
                                 "mean_y", "sd_y", "CI_y",
                                 "mean_agent", "sd_agent", "CI_agent")
raw_path_continent_se_1$year="Start"

raw_path_continent_se_2<-raw_path_continent_se[, c("group", "survive","SSP", "threshold", "continent_i", "continent_label",
                                           "mean_end_alt", "sd_end_alt", "CI_end_alt",
                                           "mean_end_y", "sd_end_y", "CI_end_y",
                                           "mean_end_agent", "sd_end_agent", "CI_end_agent")]
colnames(raw_path_continent_se_2)<-c("group", "survive","SSP", "threshold", "continent_i", "continent_label",
                                 "mean_alt", "sd_alt", "CI_alt",
                                 "mean_y", "sd_y", "CI_y",
                                 "mean_agent", "sd_agent", "CI_agent")
raw_path_continent_se_2$year="End"

raw_path_continent_se_g<-rbind(raw_path_continent_se_1, raw_path_continent_se_2)
raw_path_continent_se_g$year<-factor(raw_path_continent_se_g$year, levels = c("Start", "End"))
raw_path_continent_se_g$exposure<-" no exposure"
raw_path_continent_se_g[which(raw_path_continent_se_g$threshold==5),]$exposure<-"5-year exposure"
co<-"Africa"
for (co in unique(raw_path_continent_se_g$continent_label)){
  print(co)
  item<-raw_path_continent_se_g%>%dplyr::filter(continent_label==co)
  p1<-ggplot(item)+
    geom_point(aes(x=year, y=mean_alt, color=group))+
    geom_errorbar(aes(x=year, y=mean_alt, ymin=mean_alt-CI_alt, ymax=mean_alt+CI_alt, color=group), 
                  width = 0.1, position = "dodge2")+
    geom_line(aes(x=as.numeric(year), y=mean_alt, color=group, linetype=factor(exposure)))+
    scale_color_manual(values=color_groups)+
    facet_grid(survive~SSP, scale="free")+
    labs(x="Time spots", y="Elevation (m)", color="Group", linetype="Exposure")+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  legend<-g_legend(p1)
  
  item$mean_y<-item$mean_y/1000
  item$CI_y<-item$CI_y/1000
  p2<-ggplot(item)+
    geom_point(aes(x=year, y=mean_y, color=group))+
    geom_errorbar(aes(x=year, y=mean_y, ymin=mean_y-CI_y, ymax=mean_y+CI_y, color=group), 
                  width = 0.1, position = "dodge2")+
    geom_line(aes(x=as.numeric(year), y=mean_y, color=group, linetype=factor(exposure)))+
    scale_color_manual(values=color_groups)+
    facet_grid(survive~SSP, scale="free")+
    labs(x="Time spots", y="Latitude (km)", color="Group", linetype="Exposure")+
    theme_bw()+
    theme(axis.title.x=element_blank())
  
  pp<-ggarrange(p1, p2, nrow=2, ncol=1, common.legend=T, legend.grob=legend, legend="right")
  pp2<-annotate_figure(pp,
                  top = text_grob(co, face = "bold", size = 14),
                  bottom = text_grob("Time spots", size = 10)
  )
  ggsave(pp2, filename=sprintf("../../Figures/Top_Figure_all/gradient_%s.png", co), width=10, height=8)
  ggsave(pp2, filename=sprintf("../../Figures/Top_Figure_all/gradient_%s.pdf", co), width=10, height=8)
}

if (F){
  
  library(ggplot2)
  raw_path$agent<-raw_path$alt+abs(raw_path$gravity_y)/1000
  
  raw_path$next_year<-raw_path$YEAR+20
  raw_path_with_next<-left_join(raw_path, raw_path, by=c("continent_i", "group", "sp", "survive", "next_year"="YEAR"))
  raw_path_with_next$alt_differ<-raw_path_with_next$alt.y-raw_path_with_next$alt.x
  raw_path_with_next$y_differ<-abs(raw_path_with_next$gravity_y.y)-abs(raw_path_with_next$gravity_y.x)
  raw_path_with_next$agent_differ<-raw_path_with_next$agent.y-raw_path_with_next$agent.x
  
  
  raw_path_se<-raw_path_with_next%>%dplyr::group_by(group, survive, YEAR, next_year)%>%
    dplyr::summarise(mean_alt=mean(alt.x, na.rm=T),
                     sd_alt=sd(alt.x, na.rm=T),
                     mean_y=mean(abs(gravity_y.x), na.rm=T),
                     sd_y=sd(abs(gravity_y.x), na.rm=T),
                     sd_agent=sd(agent.x, na.rm=T),
                     mean_agent=mean(agent.x, na.rm=T),
                     mean_alt_differ=mean(alt_differ, na.rm=T),
                     sd_alt_differ=sd(alt_differ, na.rm=T),
                     mean_y_differ=mean(y_differ, na.rm=T),
                     sd_y_differ=sd(y_differ, na.rm=T),
                     mean_agent_differ=mean(agent_differ, na.rm=T),
                     sd_agent_differ=sd(agent_differ, na.rm=T))
  
  raw_path_se<-raw_path_se%>%dplyr::filter(!is.nan(mean_alt_differ))
  ggplot(raw_path_se, aes(x=next_year, y=mean_alt, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_wrap(~survive, scale="free", nrow=2)
  
  ggplot(raw_path_se, aes(x=next_year, y=mean_alt_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_wrap(~survive, scale="free", nrow=2)
  
  ggplot(raw_path_se, aes(x=next_year, y=mean_y_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_wrap(~survive, scale="free", nrow=2)
  
  ggplot(raw_path_se, aes(x=next_year, y=mean_agent_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_wrap(~survive, scale="free", nrow=2)
  
  raw_path_continent_se<-raw_path_with_next%>%dplyr::group_by(group, survive, YEAR, next_year, continent_i)%>%
    dplyr::summarise(mean_alt=mean(alt.x, na.rm=T),
                     sd_alt=sd(alt.x, na.rm=T),
                     mean_y=mean(abs(gravity_y.x), na.rm=T),
                     sd_y=sd(abs(gravity_y.x), na.rm=T),
                     sd_agent=sd(agent.x, na.rm=T),
                     mean_agent=mean(agent.x, na.rm=T),
                     mean_alt_differ=mean(alt_differ, na.rm=T),
                     sd_alt_differ=sd(alt_differ, na.rm=T),
                     mean_y_differ=mean(y_differ, na.rm=T),
                     sd_y_differ=sd(y_differ, na.rm=T),
                     mean_agent_differ=mean(agent_differ, na.rm=T),
                     sd_agent_differ=sd(agent_differ, na.rm=T))
  
  raw_path_continent_se$continent_label<-""
  raw_path_continent_se[which(raw_path_continent_se$continent_i==1), "continent_label"]<-"Africa"
  raw_path_continent_se[which(raw_path_continent_se$continent_i==2), "continent_label"]<-"Euroasia"
  raw_path_continent_se[which(raw_path_continent_se$continent_i==3), "continent_label"]<-"Australia"
  raw_path_continent_se[which(raw_path_continent_se$continent_i==5), "continent_label"]<-"North America"
  raw_path_continent_se[which(raw_path_continent_se$continent_i==6), "continent_label"]<-"South America"
  raw_path_continent_se<-raw_path_continent_se%>%dplyr::filter(!is.nan(mean_alt_differ))
  
  ggplot(raw_path_continent_se, aes(x=next_year, y=mean_alt, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_grid(continent_label~survive, scale="free")
  
  ggplot(raw_path_continent_se, aes(x=next_year, y=mean_alt_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_grid(continent_label~survive, scale="free")
  
  ggplot(raw_path_continent_se, aes(x=next_year, y=mean_y_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_grid(continent_label~survive, scale="free")
  
  ggplot(raw_path_continent_se, aes(x=next_year, y=mean_agent_differ, color=factor(group)))+
    geom_line()+
    geom_smooth(method="lm")+
    facet_grid(continent_label~survive, scale="free")
}