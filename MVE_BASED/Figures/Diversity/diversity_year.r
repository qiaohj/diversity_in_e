setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(ggplot2)
library(dplyr)
g<-"Amphibians"
f<-"Diversity"
source("functions.r")

for (f in c("Diversity", "Diversity_with_human")){
  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    print(sprintf("../../Figures/%s/%s/%s.rda", f, g, g))
    df_o<-readRDS(sprintf("../../Figures/%s/%s/%s.rda", f, g, g))
    for (metric in unique(df_o$type)){
      df<-df_o%>%filter(type==metric)

      
      df_se<-df%>%dplyr::group_by(name, type, year, M, N, SSP, group)%>%
        dplyr::summarise(mean_v=mean(mean), 
                         sd_v=mean(sd),
                         CI_v=mean(CI))
      f_l<-"With HFP"
      if (f=="Diversity"){
        f_l<-"Without HFP"
      }
      p<-ggplot(df_se, aes(x=year, y=mean_v, color=factor(M)))+
        #geom_point()+
        geom_line(aes(linetype=factor(SSP)))+
        theme_bw()+
        labs(title=paste(g, metric, f_l))+
        facet_wrap(~name, scale="free", ncol=2)
      ggsave(p, file=sprintf("../../Figures/%s/%s_%s.png", f, g, metric))
    }
  }
}


r1<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity/Amphibians/species.richness/EC-Earth3-Veg_SSP585_0_1/2014.tif")
r2<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity/Amphibians/species.richness/EC-Earth3-Veg_SSP585_0_1/2100.tif")
r3<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity/Amphibians/species.richness/EC-Earth3-Veg_SSP585_1_1/2014.tif")
r4<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity/Amphibians/species.richness/EC-Earth3-Veg_SSP585_1_1/2100.tif")
r5<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity_with_human/Amphibians/species.richness/EC-Earth3-Veg_SSP585_1_1/2014.tif")
r6<-raster("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/Diversity_with_human/Amphibians/species.richness/EC-Earth3-Veg_SSP585_1_1/2100.tif")

plot(r1)
plot(r3)
plot(r2)
plot(r2-r1)
plot(r3-r1)

plot(r4-r3)
plot(r5-r4)
plot(r6-r5)

t<-r4-r3
t1<-t
library(RColorBrewer)
loss <- colorRampPalette(c("white","blue"))
gain <- colorRampPalette(c("white","red"))

t1<-t
values(t1)[which(values(t1)>=0)]<-NA
t2<-t
values(t2)[which(values(t2)<0)]<-NA
plot(t1, col=loss(100))
plot(t2, col=gain(100), add=T)
#writeRaster(t, "~/Downloads/test.tif", overwrite=T)

t<-r-r3
t1<-t
library(RColorBrewer)
loss <- colorRampPalette(c("white","blue"))
gain <- colorRampPalette(c("white","red"))

t1<-t
values(t1)[which(values(t1)>=0)]<-NA
t2<-t
values(t2)[which(values(t2)<0)]<-NA
plot(t1, col=loss(100))
plot(t2, col=gain(100), add=T)
