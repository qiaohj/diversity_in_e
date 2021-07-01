library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(ggpubr)
setwd("Y:/Script/diversity_in_e")

source("colors.R")
source("functions.r")
source("genCircle.R")
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2), na.rm = T)
}

group<-"Amphibians"


example_sp<-"Dendropsophus_walfordi"
target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
all_v<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))

mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
#plot(mask)
target<-sprintf("%s/dispersal", target_folder)
dispersal_df<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1_1"))

year=2015
g2<-ggplot()+geom_raster(data=mask_p, aes(x=x, y=y), fill=colors_black[3])
centers<-NULL
for (year in c(2015:2099)){
  print(year)
  dis1<-dispersal_df%>%dplyr::filter(YEAR==year)
  center<-c(x=mean(dis1$x), y=mean(dis1$y), year=year)
  p1<-left_join(mask_p, dis1, by=c("x", "y", "mask_index"))
  g2<-g2+geom_raster(data=p1 %>% dplyr::filter(!is.na(YEAR)), aes(x=x, y=y), fill=colors_black[8], alpha=(year-2015)/85)
  centers<-bind(centers, center)
}
g2<-g2+geom_point(data=centers, aes(x=x, y=y, color=year))

g2<-g2+xlim(c(-8000000, -3000000))+
  ylim(c(-3300000, 1500000))+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")
g2<-g2+ggtitle("(a)")
#g2
#ggsave(g2, filename=sprintf("../../Figures/Methods/dispersal_center.png", year), width = 6, height=5)

g3<-ggplot()+geom_point(data=centers, aes(x=x, y=y, color=year))+
  geom_path(data=centers, aes(x=x, y=y, color=year))
keyyears<-round(quantile(centers$year, seq(0, 1, 0.2)))
keycenters<-centers%>%dplyr::filter(year %in% keyyears)
g3<-g3+geom_point(data=keycenters, aes(x=x, y=y), color=colors_red[8], size=1.5)+
  geom_path(data=keycenters, aes(x=x, y=y), color=colors_red[5])+
  geom_text(data=keycenters, aes(x=x, y=y, label=year),hjust=0, vjust=0, color=colors_red[5])+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")
g3<-g3+ggtitle("(b)")
plot.new()
sm_path<-data.frame(xspline(keycenters$x, keycenters$y, shape=1, draw=F, repEnds=T))
sm_path$YEAR<-seq(keyyears[1], keyyears[length(keyyears)], by=(keyyears[length(keyyears)]-keyyears[1])/(nrow(sm_path)-1))
sm_path$alpha<-((sm_path$YEAR-2020)/80)^5
g4<-ggplot()+geom_point(data=keycenters, aes(x=x, y=y), color=colors_red[8], size=1.5)+
  geom_path(data=keycenters, aes(x=x, y=y), color=colors_red[5])+
  geom_path(data=sm_path, aes(x=x, y=y), color=colors_purple[9])+
  geom_text(data=keycenters, aes(x=x, y=y, label=year),hjust=0, vjust=0, color=colors_red[5])+
  #scale_alpha_continuous()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")
g4<-g4+ggtitle("(c)")
g4

library(png)
library(grid)
img <- readPNG("../../Figures/Top_Figure/Top_Figure_ALL_UKESM1_SSP585.png")
g <- rasterGrob(img, interpolate=TRUE)

g1<-qplot(1:1000, 1:1000, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")
g1<-g1+ggtitle("(d)")
g1
g<-ggarrange(ggarrange(g2, g3, g4, nrow=1),
          g1, nrow=2, heights=c(2, 6))

ggsave(g, filename="../../Figures/Methods/dispersal_pathway.png", height=10, width=10)

#ffmpeg -r 1 -start_number 2015 -i %04d.png -y ../dispersal_diagram.gif
#ffmpeg -r 2 -start_number 2015 -i %04d.png -y ../dispersal_diagram.mp4
