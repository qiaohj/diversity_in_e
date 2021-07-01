library(raster)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
GCMs<-c(NA, "bc", "cc", "gs", "hd", "he", "ip", "mc", "mg", "mi", "mr", "no")
rcps<-c("26", "45", "60", "85")
GCM<-GCMs[2]
rcp<-rcps[1]
all_p<-NULL
for (GCM in GCMs){
  for (rcp in rcps){
    print(paste(GCM, rcp))
    if (is.na(GCM)){
      bio1<-raster("../Raster/Bioclim/Present/bio1.tif")
      bio12<-raster("../Raster/Bioclim/Present/bio12.tif")
    }else{
      bio1<-raster(sprintf("../Raster/Bioclim/Future/%s%sbi70/%s%sbi701.tif", GCM, rcp, GCM, rcp))
      bio12<-raster(sprintf("../Raster/Bioclim/Future/%s%sbi70/%s%sbi7012.tif", GCM, rcp, GCM, rcp))
    }
    p<-as_tibble(rasterToPoints(bio1))
    colnames(p)[3]<-"bio1"
    p$bio12<-raster::extract(bio12, p[, c("x", "y")])
    if (is.na(GCM)){ 
      p$GCM<-"Present"
      p$rcp<-rcp
    }else{
      p$GCM<-GCM
      p$rcp<-rcp
    }
    if (!is.null(p)){
      if (is.null(all_p)){
        all_p<-p
      }else{
        all_p<-bind_rows(all_p, p)
      }
    }
  }
}
saveRDS(all_p, "../Object/all_p_bio1_bio12.rda")

all_p<-readRDS("../Object/all_p_bio1_bio12.rda")
all_p[!complete.cases(all_p),]
all_p_se<-all_p%>%dplyr::group_by(GCM, rcp, y)%>%
  dplyr::summarise(mean_bio1=mean(bio1),
                   mean_bio12=mean(bio12))

n <- length(GCMs)-1
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors<-c(sample(col_vector, n), "black")

p_bio1<-ggplot(all_p_se, aes(x=y, color=factor(GCM)))+
  geom_line(aes(y=mean_bio1))+
  theme_bw()+
  scale_color_manual(name = "GCM", labels = c(GCMs[-1], "Present"), values=colors)+
  facet_wrap(~rcp, nrow=2)
p_bio12<-ggplot(all_p_se, aes(x=y, color=factor(GCM)))+
  geom_line(aes(y=mean_bio12))+
  theme_bw()+
  scale_color_manual(name = "GCM", labels = c(GCMs[-1], "Present"), values=colors)+
  facet_wrap(~rcp, nrow=2)
pp<-ggarrange(p_bio1, p_bio12, nrow=1, common.legend=T, legend="right")
pp2<-annotate_figure(pp,
                top = text_grob("Annual Mean Temperature (bio1) and Annual Precipitation (bio12)", 
                                face = "bold", size = 14))

ggsave(pp2, file="../Figures/future_evn_profile.png", width=10, height = 6)                
length(unique(all_p_se$y))
