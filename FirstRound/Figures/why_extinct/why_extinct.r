library(dplyr)
library(ggplot2)
library(raster)
library(Rmisc)

source("commonFuns/colors.r")
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
df<-readRDS("../../Objects_Full_species/why_extinct/why_extinct.rda")
head(df)
df$exposure<-ifelse(df$threshold==1, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==0, "no dispersal", "with dispersal")

#factor by year
df_se_temp<-df%>%dplyr::filter(!is_temp_in)%>%dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp$causation<-"Temperature"
df_se_prec<-df%>%dplyr::filter(!is_prec_in)%>%dplyr::group_by(SSP, extinct_year, da, exposure)%>%
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
ggsave(p1, filename="../../Figures_Full_species/why_extinct/by_year.pdf")

#factor by lat
range(df$y)
df$lat<-round(df$y/10e5)
df$lon<-round(df$x/10e5)

mask<-raster("../../Raster/mask_index.tif")
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

df_se_temp<-df%>%dplyr::filter(!is_temp_in)%>%dplyr::group_by(SSP, GCM, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp<-df_se_temp%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[2]-CI(N)[3])
df_se_temp[is.na(df_se_temp)]<-0

df_se_temp$causation<-"Temperature"

df_se_prec<-df%>%dplyr::filter(!is_prec_in)%>%
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
ggsave(p2, filename="../../Figures_Full_species/why_extinct/by_lat.pdf")

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
ggsave(p3, filename="../../Figures_Full_species/why_extinct/by_lat_fixed_by_area.pdf", width=8, height=4)
ggsave(p3, filename="../../Figures_Full_species/why_extinct/by_lat_fixed_by_area.png", width=8, height=4)


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
ggsave(pp, filename="../../Figures_Full_species/why_extinct/combined_why_extinct.png", width=10, height=8)
ggsave(pp, filename="../../Figures_Full_species/why_extinct/combined_why_extinct.pdf", width=10, height=8)
