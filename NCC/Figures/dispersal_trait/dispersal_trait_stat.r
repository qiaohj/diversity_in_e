library(dplyr)
library(data.table)
library(ggplot2)
library(Rmisc)
source("commonFuns/colors.r")
df_all<-readRDS("../../Objects/dispersal_trait/all.rda")
df_all$status<-ifelse(df_all$extinct_year==2101, "extant", "extinct")
df_all$exposure<-ifelse(df_all$exposure==0, " no exposure", "5-year exposure")
df_all$da<-ifelse(df_all$dispersal==0, "no dispersal", "with dispersal")
df_all_start<-df_all[type=="start"]
df_all_end<-df_all[type=="end"]
colnames(df_all_start)[c(1:5)]<-paste(colnames(df_all_start)[c(1:5)], "start", sep="_")
colnames(df_all_end)[c(1:5)]<-paste(colnames(df_all_end)[c(1:5)], "end", sep="_")

tropic_ll<-data.table(lon=0, lat=c(23.43654, -23.43654))
mask_ll<-raster("../../Raster/Continent.tif")
mask_eck4<-raster("../../Raster/Continent_ect4.tif")
p_tropic_all<-SpatialPointsDataFrame(tropic_ll, tropic_ll, proj4string=crs(mask_ll))
p_tropic_all_eck4<-spTransform(p_tropic_all, crs(mask_eck4))
lat_exposure<-p_tropic_all_eck4@coords[,2]


df_all_merged<-inner_join(df_all_start, df_all_end, 
                          by=colnames(df_all_start)[c(6:ncol(df_all_start))])
dim(df_all_merged)
dim(df_all_start)

df_all_merged$change_alt<-ifelse(df_all_merged$mean_alt_start>df_all_merged$mean_alt_end,
                                 "to lower elevation", "to upper elevation")
df_all_merged$change_lat<-ifelse(df_all_merged$mean_abs_y_start>df_all_merged$mean_abs_y_end,
                                 "to lower latitude", "to upper latitude")

df_all_merged$is_tropic<-"tropics"
df_all_merged[(df_all_merged$mean_y_start>lat_exposure[1]), "is_tropic"]<-"north temperate"
df_all_merged[(df_all_merged$mean_y_start<lat_exposure[2]), "is_tropic"]<-"south temperate"

df_all_merged$label<-"neither higher elevation nor higher latitude"
df_all_merged[which(df_all_merged$change_alt=="to upper elevation"|df_all_merged$change_lat=="to upper latitude"),
                 "label"]<-"higher elevation and/or higher latitude"

df_se<-df_all_merged%>%group_by(SSP, GCM, label, da, exposure, status)%>%
  dplyr::summarise(N=n())
unique(df_se$exposure)
df_se<-df_se%>%dplyr::filter(da=="with dispersal")


df_se_across_GCM<-df_se%>%dplyr::group_by(SSP, label, da, status)%>%
  dplyr::summarise(mean_N=mean(N),
                   sd_N=sd(N),
                   CI_N=CI(N)[1]-CI(N)[2])


p<-ggplot(df_se_across_GCM)+
  geom_errorbar(aes(x=SSP, ymin=mean_N/1000-sd_N/1000, ymax=mean_N/1000+sd_N/1000, 
                    color=label), width=0.1)+
  geom_point(aes(x=SSP, y=mean_N/1000, color=label))+
  facet_wrap(~status, scale="free_y", nrow=2, strip.position="right")+
  scale_color_manual(values=color_dipsersal_type)+
  theme_bw()+
  ylab("Number of species (×1000)")+
  labs(color="")+
  theme(legend.position="top",
        legend.direction = "vertical",
        axis.line=element_line())
p
ggsave(p, filename = "../../Figures/dispersal_traits/dispersal_stat.png", width=3, height=4)
SSS<-"SSP119"
for (SSS in unique(df_se_across_GCM$SSP)){
  p<-ggplot(df_se_across_GCM%>%dplyr::filter(SSP==SSS))+
    geom_errorbar(aes(x=exposure, ymin=mean_N/1000-sd_N/1000, ymax=mean_N/1000+sd_N/1000, color=label), width=0.1)+
    geom_point(aes(x=exposure, y=mean_N/1000, color=label))+
    facet_wrap(~status, scale="free_y", nrow=2, strip.position="right")+
    scale_color_manual(values=color_dipsersal_type)+
    theme_bw()+
    ylab("Number of species (×1000)")+
    labs(color="")+
    xlab("")+
    theme(legend.position="top",
          legend.direction = "vertical",
          axis.line=element_line())
  p
  ggsave(p, filename = sprintf("../../Figures/dispersal_traits/dispersal_stat_%s.png", SSS), width=3, height=4)
}
hist(df_se_across_GCM$mean_N, col=factor(df_se_across_GCM$status))
write.table(df_se_across_GCM, "../../Figures/dispersal_traits/dispersal_stat.csv", sep=",", row.names = F)
ggplot(df_se_across_GCM)+geom_point(aes(x=SSP, y=mean_N, color=factor(status), shape=factor(exposure)))+
  facet_grid(change_alt~change_lat)
  
df_se_across_GCM$label<-paste(df_se_across_GCM$change_alt, df_se_across_GCM$change_lat)

ggplot(df_se_across_GCM)+geom_point(aes(x=SSP, y=mean_N, color=factor(status)))+
  facet_grid(exposure~label)

df_se<-df_all_merged%>%group_by(SSP, GCM, change_alt, change_lat, da, exposure, status, is_tropic)%>%
  dplyr::summarise(N=n())

df_se<-df_se%>%dplyr::filter(da=="with dispersal")


df_se_across_GCM<-df_se%>%dplyr::group_by(SSP, change_alt, change_lat, da, exposure, status, is_tropic)%>%
  dplyr::summarise(mean_N=mean(N))

ggplot(df_se_across_GCM%>%dplyr::filter(exposure=="5-year exposure"))+
  geom_point(aes(x=SSP, y=mean_N, color=factor(is_tropic), shape=factor(status)))+
  facet_grid(change_alt~change_lat, scale="free")


  
