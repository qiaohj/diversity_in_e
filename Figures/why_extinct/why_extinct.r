library(dplyr)
library(ggplot2)
library(raster)
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
p<-ggplot(df_se)+geom_line(aes(x=extinct_year, y=N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("Year")+
  ylab("N extinct events")+
  labs(color="Causation")+
  theme_bw()
ggsave(p, filename="../../Figures_Full_species/why_extinct/by_year.pdf")

#factor by lat
range(df$y)
df$lat<-round(df$y/10e5)
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
mask_p$lat<-round(mask_p$y/10e5)
lat_N<-data.frame(table(mask_p$lat))
colnames(lat_N)<-c("lat", "lat_N")
lat_N$lat<-as.numeric(as.character(lat_N$lat))

df<-inner_join(df, lat_N, by="lat")
df_se_temp<-df%>%dplyr::filter(!is_temp_in)%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp$causation<-"Temperature"
df_se_prec<-df%>%dplyr::filter(!is_prec_in)%>%dplyr::group_by(SSP, lat, lat_N, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec$causation<-"Precipitation"

df_se<-bind_rows(df_se_prec, df_se_temp)

p<-ggplot(df_se)+geom_line(aes(x=lat, y=N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("N extinct events")+
  labs(color="Causation")+
  theme_bw()
ggsave(p, filename="../../Figures_Full_species/why_extinct/by_lat.pdf")

p<-ggplot(df_se)+geom_line(aes(x=lat, y=N/lat_N, color=causation, linetype=SSP))+
  scale_color_manual(values=color_causation)+
  facet_grid(exposure~da, scale="free")+
  xlab("latitudinal band (in 100km)")+
  ylab("mean N extinct events")+
  labs(color="Causation")+
  theme_bw()
plot(df_se$lat, df_se$lat_N)
ggsave(p, filename="../../Figures_Full_species/why_extinct/by_lat_fixed_by_area.pdf")
