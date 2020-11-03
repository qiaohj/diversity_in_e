setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sp)
library(Rmisc)
library(ggpubr)
f<-"Diversity"
g<-"Amphibians"

for (f in c("Diversity", "Diversity_with_human")){
  for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
    print(paste(f, g))
    df<-readRDS(sprintf("../../Figures/Species_gain_loss/%s_%s.rda", f, g))
    if (F){
      df_sample<-df[sample(nrow(df), 1000),]
      plot(df_sample$slope, df_sample$mean_gain_loss)
    }
    
    source("functions.r")
    
    keyspots<-read.csv("../../Objects/Species_exposure/Keyspots/keyspots.csv", stringsAsFactors = F)
    mask<-raster("../../Raster/mask_index.tif")
    plot(mask)
    cord.dec = SpatialPoints(keyspots[, c("lon", "lat")], proj4string=CRS("+proj=longlat"))
    cord.eck4 <- spTransform(cord.dec, crs(mask))
    plot(cord.eck4, add=T)
    keyspots$x<-cord.eck4@coords[,1]
    keyspots$y<-cord.eck4@coords[,2]
    
    
    
    pc100km<-readOGR("../../Shape/pc100km", "pc100km") 
    plot(pc100km, add=T)
    points<-SpatialPoints(df[, c("x", "y")], proj4string=crs(mask))
    crs(points)<-crs(pc100km)
    i=1
    
    diversity<-readRDS(sprintf("../../Figures/%s/%s/%s.rda", f, g, g))
    diversity_se<-diversity%>%dplyr::group_by(group, SSP, M, N, year, type, name)%>%
      dplyr::summarise(mean_mean=mean(mean))
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    diversity_item<-diversity_se%>%dplyr::filter(type=="species.richness")
    df[which((df$M==0)&(df$mean_gain_loss>0)), "mean_gain_loss"]<-0
    df$gain_loss<-"GAIN"
    df[which(df$mean_gain_loss<0), ]$gain_loss<-"LOSS"
    df[which(df$mean_gain_loss==0), ]$gain_loss<-"UNCHANGED"
    
    all_se<-df%>%dplyr::group_by(type, group, SSP, M, N, year, gain_loss)%>%
      dplyr::summarise(mean_v_mean=mean(mean_v, na.rm=T),
                       mean_v_sd=sd(mean_v, na.rm=T),
                       mean_v_CI=CI(mean_v)[2]-CI(mean_v)[3],
                       mean_gain_loss_mean=mean(mean_gain_loss, na.rm=T),
                       mean_gain_loss_sd=sd(mean_gain_loss, na.rm=T),
                       mean_gain_loss_CI=CI(mean_gain_loss)[2]-CI(mean_gain_loss)[3],
                       N_pixel=n())
    
    p1<-ggplot(all_se, aes(x=year, y=mean_gain_loss_mean, color=factor(gain_loss)))+
      #geom_point(aes(shape=factor(SSP)))+
      geom_line(aes(linetype=factor(SSP)))+
      #geom_line(data=diversity_item, aes(y=mean_mean, linetype=factor(SSP)), color="black")+
      facet_wrap(~M, ncol=1, scale="free")+theme_bw()
    legend<-g_legend(p1)
    p1<-p1+theme(legend.position="none")
    p2<-ggplot(all_se, aes(x=year, y=N_pixel, color=factor(gain_loss)))+
      #geom_point(aes(shape=factor(SSP)))+
      geom_line(aes(linetype=factor(SSP)))+
      #geom_line(data=diversity_item, aes(y=mean_mean, linetype=factor(SSP)), color="black")+
      facet_wrap(~M, ncol=1, scale="free")+theme_bw()+
      theme(legend.position="none")
    p<-ggarrange(p1, p2, nrow=1, legend)
    p<-annotate_figure(p,
                       top = text_grob(paste(f, g, "Global"), color = "black", size = 14)
    )
    ggsave(p, file=sprintf("../../Figures/Species_gain_loss/gain_loss_by_year/%s_%s_%s.png", f, g, "Global"),
           width=12, height=10)
    
    for (i in c(1:nrow(keyspots))){
      
      keyspot<-keyspots[i,]
      print(paste(f, g, keyspot$name))
      over_points<-sp::over(points, pc100km[which(pc100km$name==keyspot$name),])
      over_points_result<-df[!is.na(over_points$name),]
      
      
      over_points_result_se<-over_points_result%>%dplyr::group_by(type, group, SSP, M, N, year, gain_loss)%>%
        dplyr::summarise(mean_v_mean=mean(mean_v, na.rm=T),
                         mean_v_sd=sd(mean_v, na.rm=T),
                         mean_v_CI=CI(mean_v)[2]-CI(mean_v)[3],
                         mean_gain_loss_mean=mean(mean_gain_loss, na.rm=T),
                         mean_gain_loss_sd=sd(mean_gain_loss, na.rm=T),
                         mean_gain_loss_CI=CI(mean_gain_loss)[2]-CI(mean_gain_loss)[3],
                         N_pixel=n())
      over_points_result_se$name<-keyspot$name
      
      over_points_result_se%>%dplyr::filter((SSP=="SSP119")&(M==1)&(year==2015))
      diversity_item<-diversity_se%>%dplyr::filter((type=="species.richness")&(name==keyspot$name))
      
      p1<-ggplot(over_points_result_se, aes(x=year, y=mean_gain_loss_mean, color=factor(gain_loss)))+
        #geom_point(aes(shape=factor(SSP)))+
        geom_line(aes(linetype=factor(SSP)))+
        #geom_line(data=diversity_item, aes(y=mean_mean, linetype=factor(SSP)), color="black")+
        facet_wrap(~M, ncol=1, scale="free")+theme_bw()
      legend<-g_legend(p1)
      p1<-p1+theme(legend.position="none")
      p2<-ggplot(over_points_result_se, aes(x=year, y=N_pixel, color=factor(gain_loss)))+
        #geom_point(aes(shape=factor(SSP)))+
        geom_line(aes(linetype=factor(SSP)))+
        #geom_line(data=diversity_item, aes(y=mean_mean, linetype=factor(SSP)), color="black")+
        facet_wrap(~M, ncol=1, scale="free")+theme_bw()+
        theme(legend.position="none")
      p<-ggarrange(p1, p2, nrow=1, legend)
      p<-annotate_figure(p,
                         top = text_grob(paste(f, g, keyspot$name), color = "black", size = 14)
      )
      ggsave(p, file=sprintf("../../Figures/Species_gain_loss/gain_loss_by_year/%s_%s_%s.png", f, g, keyspot$name),
             width=12, height=10)
    }
  }
}