library(dplyr)
library(ggplot2)
library(raster)
library(Rmisc)
library(ggpubr)

source("commonFuns/colors.r")
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
df<-readRDS("../../Objects/why_extinct/why_extinct.rda")
head(df)
df$exposure<-ifelse(df$exposure==0, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==0, "no dispersal", "with dispersal")

#factor by year
df_se_temp<-df%>%dplyr::filter((!is_bio1)|(!is_bio5)|(!is_bio6))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp$causation<-"Temperature"
df_se_prec<-df%>%dplyr::filter((!is_bio12)|(!is_bio13)|(!is_bio14))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec$causation<-"Precipitation"

df_se<-bind_rows(df_se_prec, df_se_temp)

df_se_prec%>%ungroup()%>%dplyr::filter((N>1000))
p1<-ggplot(df_se)+geom_line(aes(x=extinct_year, y=N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("Year")+
  ylab("Number of extinction events")+
  labs(color="Causation")+
  scale_y_sqrt(labels=seq(0, 45, by=5)^2, breaks=seq(0, 45, by=5)^2)+
  theme_bw()
p1
ggsave(p1, filename="../../Figures/why_extinct/by_year.pdf")

#factor by lat
range(df$y)
df$lat<-round(df$y/10e5)
df$lon<-round(df$x/10e5)

mask<-raster("../../Raster/mask_100km.tif")
mask_p<-data.frame(rasterToPoints(mask))
mask_p$lat<-round(mask_p$y/10e5)
mask_p$lon<-round(mask_p$x/10e5)

lat_N<-data.frame(table(mask_p$lat))
colnames(lat_N)<-c("lat", "lat_N")
lat_N$lat<-as.numeric(as.character(lat_N$lat))

lon_N<-data.frame(table(mask_p$lon))
colnames(lon_N)<-c("lon", "lon_N")
lon_N$lon<-as.numeric(as.character(lon_N$lon))


df<-inner_join(df, lat_N, by="lat")
df<-inner_join(df, lon_N, by="lon")

df_se_temp<-df%>%dplyr::filter((!is_bio1)|(!is_bio5)|(!is_bio6))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp<-df_se_temp%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_temp[is.na(df_se_temp)]<-0

df_se_temp$causation<-"Temperature"

df_se_prec<-df%>%dplyr::filter((!is_bio12)|(!is_bio13)|(!is_bio14))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec<-df_se_prec%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_prec[is.na(df_se_prec)]<-0
df_se_prec$causation<-"Precipitation"

df_se<-bind_rows(df_se_prec, df_se_temp)

p2<-ggplot(df_se)+
  geom_errorbar(aes(x=lat, ymin=mean_N-sd_N, ymax=mean_N+sd_N, color=causation, linetype=SSP))+
  geom_line(aes(x=lat, y=mean_N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  scale_fill_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("Mean number of extinction events")+
  labs(color="Causation")+
  theme_bw()
p2
ggsave(p2, filename="../../Figures/why_extinct/by_lat.pdf")

p3<-ggplot(df_se)+
  geom_errorbar(aes(x=lat, ymin=mean_N/lat_N-sd_N/lat_N, ymax=mean_N/lat_N+sd_N/lat_N), alpha=0.2)+
  geom_line(aes(x=lat, y=mean_N/lat_N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("Mean number of extinction events per grid cell")+
  labs(color="Causation")+
  theme_bw()
plot(df_se$lat, df_se$lat_N)
p3
ggsave(p3, filename="../../Figures/why_extinct/by_lat_fixed_by_area.pdf", width=8, height=4)
ggsave(p3, filename="../../Figures/why_extinct/by_lat_fixed_by_area.png", width=8, height=4)


legends<-get_legend(p1)
p1_formatted<-p1+theme(legend.position = "none", strip.background.x = element_blank(),
                       strip.text.x = element_blank())
p3_formatted<-p3+theme(legend.position = "none")

pp<-ggarrange(p3_formatted, p1_formatted, common.legend = T, legend = "right", 
          legend.grob = legends, labels=c("(a)", "(b)"), nrow=2, ncol=1,
          label.x = 0.06, label.y=0.93)

pp<-ggarrange(p3_formatted, p1_formatted, common.legend = T, legend = "right", 
              legend.grob = legends, nrow=2, ncol=1)

pp

ggsave(pp, filename="../../Figures/why_extinct/combined_why_extinct.png", width=10, height=8)
ggsave(pp, filename="../../Figures/why_extinct/combined_why_extinct.pdf", width=10, height=8)


df<-readRDS("../../Objects/why_extinct/why_extinct.rda")
head(df)
df$exposure<-ifelse(df$exposure==0, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==0, "no dispersal", "with dispersal")

#factor by year
df_se_temp_upper<-df%>%dplyr::filter((upper_bio1)|(upper_bio5)|(upper_bio6))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp_upper$causation<-"Temperature (upper limit)"

df_se_temp_lower<-df%>%dplyr::filter((lower_bio1)|(lower_bio5)|(lower_bio6))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp_lower$causation<-"Temperature (lower limit)"

df_se_prec_upper<-df%>%dplyr::filter((upper_bio12)|(upper_bio13)|(upper_bio14))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec_upper$causation<-"Precipitation (upper limit)"
df_se_prec_lower<-df%>%dplyr::filter((lower_bio12)|(lower_bio13)|(lower_bio14))%>%
  dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec_lower$causation<-"Precipitation (lower limit)"

df_se<-bind_rows(bind_rows(df_se_temp_upper, df_se_temp_lower), 
                 bind_rows(df_se_prec_upper, df_se_prec_lower))

p1<-ggplot(df_se)+geom_line(aes(x=extinct_year, y=N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation_ul)+
  facet_grid(exposure~da, scale="free")+
  xlab("Year")+
  ylab("Number of extinction events")+
  labs(color="Causation")+
  scale_y_sqrt(labels=seq(0, 45, by=5)^2, breaks=seq(0, 45, by=5)^2)+
  theme_bw()
p1
ggsave(p1, filename="../../Figures/why_extinct/by_year_ul.pdf")

p1_2<-ggplot(df_se)+geom_line(aes(x=extinct_year, y=N, linetype=SSP))+
  #scale_color_manual(values=color_causation_ul)+
  facet_grid(exposure+causation~da, scale="free")+
  xlab("Year")+
  ylab("Number of extinction events")+
  labs(color="Causation")+
  scale_y_sqrt(labels=seq(0, 45, by=2)^2, breaks=seq(0, 45, by=2)^2)+
  theme_bw()
p1_2
ggsave(p1_2, filename="../../Figures/why_extinct/by_year_ul_splitted.pdf", width=10, height=15)
ggsave(p1_2, filename="../../Figures/why_extinct/by_year_ul_splitted.png", width=10, height=15)


#factor by lat
range(df$y)
df$lat<-round(df$y/10e5)
df$lon<-round(df$x/10e5)

mask<-raster("../../Raster/mask_100km.tif")
mask_p<-data.frame(rasterToPoints(mask))
mask_p$lat<-round(mask_p$y/10e5)
mask_p$lon<-round(mask_p$x/10e5)

lat_N<-data.frame(table(mask_p$lat))
colnames(lat_N)<-c("lat", "lat_N")
lat_N$lat<-as.numeric(as.character(lat_N$lat))

lon_N<-data.frame(table(mask_p$lon))
colnames(lon_N)<-c("lon", "lon_N")
lon_N$lon<-as.numeric(as.character(lon_N$lon))


df<-inner_join(df, lat_N, by="lat")
df<-inner_join(df, lon_N, by="lon")

df_se_temp_upper<-df%>%dplyr::filter((upper_bio1)|(upper_bio5)|(upper_bio6))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp_upper<-df_se_temp_upper%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_temp_upper[is.na(df_se_temp_upper)]<-0

df_se_temp_upper$causation<-"Temperature (upper limit)"

df_se_temp_lower<-df%>%dplyr::filter((lower_bio1)|(lower_bio5)|(lower_bio6))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp_lower<-df_se_temp_lower%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_temp_lower[is.na(df_se_temp_lower)]<-0

df_se_temp_lower$causation<-"Temperature (lower limit)"

df_se_prec_upper<-df%>%dplyr::filter((upper_bio1)|(upper_bio5)|(upper_bio6))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec_upper<-df_se_prec_upper%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_prec_upper[is.na(df_se_prec_upper)]<-0

df_se_prec_upper$causation<-"Precipitation (upper limit)"

df_se_prec_lower<-df%>%dplyr::filter((lower_bio1)|(lower_bio5)|(lower_bio6))%>%
  dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec_lower<-df_se_prec_lower%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_prec_lower[is.na(df_se_prec_lower)]<-0

df_se_prec_lower$causation<-"Precipitation (lower limit)"



df_se<-bind_rows(bind_rows(df_se_temp_upper, df_se_temp_lower), 
                 bind_rows(df_se_prec_upper, df_se_prec_lower))


p2<-ggplot(df_se)+
  geom_errorbar(aes(x=lat, ymin=mean_N-sd_N, ymax=mean_N+sd_N, color=causation, linetype=SSP))+
  geom_line(aes(x=lat, y=mean_N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation_ul)+
  scale_fill_manual(values=color_causation_ul)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("Mean number of extinction events")+
  labs(color="Causation")+
  scale_y_sqrt(labels=seq(0, 45, by=5)^2, breaks=seq(0, 45, by=5)^2)+
  theme_bw()
p2
ggsave(p2, filename="../../Figures/why_extinct/by_lat_ul.pdf")

p3<-ggplot(df_se)+
  geom_errorbar(aes(x=lat, ymin=mean_N/lat_N-sd_N/lat_N, ymax=mean_N/lat_N+sd_N/lat_N), alpha=0.2)+
  geom_line(aes(x=lat, y=mean_N/lat_N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation_ul)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("Mean number of extinction events per grid cell")+
  labs(color="Causation")+
  theme_bw()
plot(df_se$lat, df_se$lat_N)
p3
ggsave(p3, filename="../../Figures/why_extinct/by_lat_fixed_by_area_ul.pdf", width=8, height=4)
ggsave(p3, filename="../../Figures/why_extinct/by_lat_fixed_by_area_ul.png", width=8, height=4)


legends<-get_legend(p1)
p1_formatted<-p1+theme(legend.position = "none", strip.background.x = element_blank(),
                       strip.text.x = element_blank())
p3_formatted<-p3+theme(legend.position = "none")

pp<-ggarrange(p3_formatted, p1_formatted, common.legend = T, legend = "right", 
              legend.grob = legends, labels=c("(a)", "(b)"), nrow=2, ncol=1,
              label.x = 0.06, label.y=0.93)

pp<-ggarrange(p3_formatted, p1_formatted, common.legend = T, legend = "right", 
              legend.grob = legends, nrow=2, ncol=1)

pp

ggsave(pp, filename="../../Figures/why_extinct/combined_why_extinct_ul.png", width=10, height=8)
ggsave(pp, filename="../../Figures/why_extinct/combined_why_extinct_ul.pdf", width=10, height=8)

