library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(scales)
library(rgdal)
library(rgeos)
library(Rmisc)
library(geometry)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")

tttt=2

n_ext_final<-readRDS("../../Figures_Full_species/when_where_extinction_all/n_ext_final_with_ttt.rda")
n_ext_final_1<-n_ext_final%>%dplyr::filter((dispersal==0)&
                                             (ttt==tttt)&
                                             (SSP=="SSP585")&
                                             (threshold==1))
n_ext_final_2<-n_ext_final%>%dplyr::filter((dispersal==1)&
                                             (ttt==tttt)&
                                             (SSP=="SSP585")&
                                             (threshold==5))
n_ext_final<-bind_rows(n_ext_final_1, n_ext_final_2)
n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
n_ext_final$exposure<-ifelse(n_ext_final$threshold==1, " no exposure", "5-year exposure")
n_ext_final$da<-ifelse(n_ext_final$dispersal==0, "no dispersal", "with dispersal")
mask<-raster("../../Raster/mask_index.tif")

mask_p<-data.frame(rasterToPoints(mask))
n_ext_final$label<-paste(n_ext_final$SSP, n_ext_final$exposure, n_ext_final$da, sep=", ")
n_ext_final$label<-factor(n_ext_final$label, 
                             levels = c("SSP585,  no exposure, no dispersal",
                                        "SSP585, 5-year exposure, with dispersal"))

if (T){
keyspots<-read.csv("../../Objects_Full_species/keyspots/keyspots.csv", stringsAsFactors = F)
cord.dec = SpatialPoints(keyspots[, c("lon", "lat")], proj4string=CRS("+proj=longlat"))
cord.eck4 <- spTransform(cord.dec, crs(mask))
plot(mask)
plot(cord.eck4, add=T)
keyspots$x<-cord.eck4@coords[,1]
keyspots$y<-cord.eck4@coords[,2]

pc100km <- gBuffer( cord.eck4, width=100*10e3, byid=TRUE )
plot(pc100km, add=T, col="red")
pc100km <- SpatialPolygonsDataFrame(pc100km, data=keyspots)

mask_p_sp<-SpatialPointsDataFrame(mask_p[, c("x", "y")], data=mask_p, proj4string=crs(mask))

mask_p_sp$keysport<-NA
for (i in c(1:nrow(pc100km))){
  overed<-!is.na(sp::over(mask_p_sp, pc100km[i,])[,1])
  mask_p_sp[overed, "keysport"]<-pc100km[i,]$name
}
table(mask_p_sp$keysport)

pc100km@data$id = rownames(pc100km@data)
pc100km.points = fortify(pc100km, region="id")
pc100km.df = join(pc100km.points[], pc100km@data[, c("name", "id")], by="id")

p<-ggplot()+
  geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  geom_tile(data=n_ext_final, aes(x=x, y=y, fill=sum_V))+
  geom_polygon(data=pc100km.df, aes(x=long, y=lat, color=factor(name)), fill=NA)+
  facet_wrap(~label)+
  scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                      limits=c(0, 100), oob=squish,
                      breaks=c(0, 20, 40, 60, 80, 100),
                      labels=c("0", "20", "40", "60", "80", 
                               sprintf(">100, up to %d", round(max(n_ext_final$sum_V)))))+
  labs(fill = "Number of extinctions", color="Key spots")+
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.background = element_rect(fill = map_background, color = NA), 
    panel.background = element_blank(), 
    #legend.background = element_rect(fill = map_background, color = NA),
    #panel.border = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    strip.background = element_blank()
  )
p

head(n_ext_final)
mask_p_sp_df<-mask_p_sp@data%>%dplyr::filter(!is.na(keysport))
head(mask_p_sp)

n_ext_final_with_keyspots<-left_join(n_ext_final, mask_p_sp_df[, c("mask_index", "keysport")], by="mask_index")


p2<-ggplot()+
  stat_density(data=n_ext_final_with_keyspots, aes(x=sum_V, y=..count..), fill="grey50", alpha=0.2)+
  stat_density(data=n_ext_final_with_keyspots%>%filter(!is.na(keysport)), 
               aes(x=sum_V, y=..count.., fill=factor(keysport)), position = "identity", alpha=0.3)+
  stat_density(data=n_ext_final_with_keyspots%>%filter(!is.na(keysport)), 
               aes(x=sum_V, y=..count.., color=factor(keysport)), position = "identity", fill=NA)+
  scale_x_log10()+
  scale_y_sqrt(breaks=seq(0, 100, by=20)^2)+
  #scale_y_continuous()+
  theme_bw()+
  labs(x="Number of extinctions", y="")+
  facet_wrap(~label, scale="free")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.background = element_rect(fill = map_background, color = NA), 
    panel.background = element_blank(),
    axis.title.y=element_blank()
  )

p2
}
g_legend<-get_legend(p)
p3<-ggarrange(p, p2, nrow=2, ncol=1, common.legend=T, legend = "right", legend.grob = g_legend)
ggsave(p3, filename="../../Figures_Full_species/Combined_Figure/keyspots_extinction.pdf", width=12, height=8)
