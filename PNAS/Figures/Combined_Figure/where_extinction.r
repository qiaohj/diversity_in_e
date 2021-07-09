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

SSPi<-"SSP245"
n_ext_final<-readRDS("../../Figures/when_where_extinction_all/n_ext_final.rda")
n_ext_final_1<-n_ext_final%>%dplyr::filter((dispersal==0)&
                                             (SSP==SSPi)&
                                             (exposure==0))
n_ext_final_2<-n_ext_final%>%dplyr::filter((dispersal==1)&
                                             (SSP==SSPi)&
                                             (exposure==5))
n_ext_final<-bind_rows(n_ext_final_1, n_ext_final_2)
n_ext_final<-n_ext_final%>%dplyr::filter(sum_V>0)
n_ext_final$exposure<-ifelse(n_ext_final$exposure==0, " no exposure", "5-year exposure")
n_ext_final$da<-ifelse(n_ext_final$dispersal==0, "no dispersal", "with dispersal")
mask<-raster("../../Raster/mask_100km_plot.tif")

mask_p<-data.frame(rasterToPoints(mask))
colnames(mask_p)[3]<-"mask_100km"
n_ext_final$label<-paste(n_ext_final$SSP, n_ext_final$exposure, n_ext_final$da, sep=", ")
n_ext_final$label<-factor(n_ext_final$label, 
                          levels = c(sprintf("%s,  no exposure, no dispersal", SSPi),
                                     sprintf("%s, 5-year exposure, with dispersal", SSPi)))

if (T){
  keyspots<-read.csv("../../Objects/keyspots/keyspots.csv", stringsAsFactors = F)
  cord.dec = SpatialPoints(keyspots[, c("lon", "lat")], proj4string=CRS("+proj=longlat"))
  cord.eck4 <- spTransform(cord.dec, crs(mask))
  plot(mask)
  plot(cord.eck4, add=T)
  keyspots$x<-cord.eck4@coords[,1]
  keyspots$y<-cord.eck4@coords[,2]
  
  pc100km <- gBuffer( cord.eck4, width=100*10e3, byid=TRUE )
  plot(pc100km, add=T, col="red")
  pc100km <- SpatialPolygonsDataFrame(pc100km, data=keyspots)
  
  mask_p_sp<-SpatialPointsDataFrame(mask_p[, c("x", "y")], 
                                    data=mask_p, proj4string=crs(mask))
  
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
    geom_polygon(data=pc100km.df, aes(x=long, y=lat, color=name), fill=NA)+
    facet_wrap(~label)+
    
    scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                        limits=c(0, 80), oob=squish,
                        breaks=c(0, 20, 40, 60, 80),
                        labels=c("0", "20", "40", "60", 
                                 sprintf(">80, up to %d", round(max(n_ext_final$sum_V)))))+
    scale_color_manual(values=color_keyspots)+
    labs(fill = "Number of extinctions", color="Key locations")+
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
      #panel.border = element_rect(colour = "black", fill=NA),
      #strip.background = element_blank()
    )
  p
  
  head(n_ext_final)
  mask_p_sp_df<-mask_p_sp@data%>%dplyr::filter(!is.na(keysport))
  head(mask_p_sp)
  
  plot(mask_p_sp$x, mask_p_sp$y)
  n_ext_final_with_keyspots<-left_join(n_ext_final, 
                                       mask_p_sp_df[, c("mask_100km", "keysport")], 
                                       by="mask_100km")
  N_keys<-unique(n_ext_final_with_keyspots$keysport)
  N_keys<-N_keys[!is.na(N_keys)]
  N_keys<-as.character(N_keys)
  table(n_ext_final_with_keyspots$dispersal, n_ext_final_with_keyspots$exposure)
  p2<-list()
  ll<- unique(n_ext_final_with_keyspots$label)[1]
  for (ll in unique(n_ext_final_with_keyspots$label)){
    item_df<-n_ext_final_with_keyspots%>%dplyr::filter(label==ll)
    binwidth<-0.1
    p2_item<-ggplot()+
      geom_density(data=item_df, aes(x=sum_V, y=..density.. * nrow(item_df) * binwidth), 
                   bw=binwidth, fill="grey50", color="grey50", alpha=0.2)
    p2_item<-ggplot()
    i=1
    for (i in c(1:length(N_keys))){
      item_df_sub<-item_df%>%dplyr::filter(keysport==N_keys[i])
      p2_item<-p2_item+geom_density(data=item_df_sub, 
                                    aes(x=sum_V, y=..density.. * nrow(item_df_sub) * binwidth, fill=keysport, color=keysport), 
                                    bw=binwidth, position = "identity", alpha=0.3)+
        geom_density(data=item_df_sub, 
                     aes(x=sum_V, y=..density.. * nrow(item_df_sub) * binwidth, color=keysport), 
                     bw=binwidth, position = "identity", fill=NA)
    }
    p2_item<-p2_item+
      scale_x_log10()+
      #scale_x_sqrt(breaks=seq(0, 14, by=2)^2)+
      #scale_y_log10()+
      #scale_y_sqrt(breaks=seq(0, 20, by=5)^2, labels=seq(0, 20, by=5)^2)+
      #scale_y_continuous()+
      theme_bw()+
      scale_color_manual(values=color_keyspots)+
      scale_fill_manual(values=color_keyspots)+
      labs(x="Number of extinction events (log transformed)", y="Number of cells (sq root transformed)")+
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.background = element_rect(fill = map_background, color = NA), 
        panel.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
        #axis.title.y=element_blank()
      )
    p2_insert<-ggplot()+geom_density(data=item_df, aes(x=sum_V, y=..density.. * nrow(item_df) * binwidth), 
                                         bw=binwidth, fill="grey50", color="grey50", alpha=0.2)+
      scale_x_log10()+
      #scale_x_sqrt(breaks=seq(0, 14, by=2)^2)+
      #scale_y_log10()+
      #scale_y_sqrt(breaks=seq(0, 20, by=5)^2, labels=seq(0, 20, by=5)^2)+
      #scale_y_continuous()+
      theme_bw()+
      scale_color_manual(values=color_keyspots)+
      scale_fill_manual(values=color_keyspots)+
      labs(x="Number of extinction events (log transformed)", y="Number of cells (sq root transformed)")+
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.background = element_rect(fill = map_background, color = NA), 
        panel.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
        #axis.title.y=element_blank()
      )
    ggsave(p2_insert, filename=sprintf("../../Figures/Combined_Figure/keyspots_extinction_%s_%s.png", 
                                       SSPi, gsub(",", "", ll)), width=6, height=3)
    
    p2[[ll]]<-p2_item
  }
}
p2_c<-ggarrange(plotlist=p2, nrow=1, ncol=2)

p2_c<-annotate_figure(p2_c,
                bottom = text_grob("Number of extinction events (log transformed)", 
                                   hjust = 2.2, x = 1, size = 10),
                left = text_grob("Number of cells (sq root transformed)", rot = 90, size = 10)
)

g_legend<-get_legend(p)
p3<-ggarrange(p, p2_c, nrow=2, ncol=1, common.legend=T, legend = "right", legend.grob = g_legend)
p3
ggsave(p3, filename=sprintf("../../Figures/Combined_Figure/keyspots_extinction_%s.pdf", SSPi), width=12, height=6)
ggsave(p3, filename=sprintf("../../Figures/Combined_Figure/keyspots_extinction_%s.png", SSPi), width=12, height=6)

