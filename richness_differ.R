library(raster)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
GCMs<-c(NA, "bc", "cc", "gs", "hd", "he", "ip", "mc", "mg", "mi", "mr", "no")
rcps<-c("26", "45", "60", "85")
groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")
i=2
GCM<-GCMs[2]
rcp<-rcps[1]
all_p<-NULL

for (i in c(1:length(groups))){
  for (GCM in GCMs){
    for (rcp in rcps){
      print(paste(GCM, rcp, groups[i]))
      if (is.na(GCM)){
        richness_present<-raster(sprintf("../Raster/Niche_E/%s_present.tif", groups[i]))
        p<-NULL
        p_present<-as_tibble(rasterToPoints(richness_present))
        colnames(p_present)[3]<-"richness_present"
      }else{
        richness<-raster(sprintf("../Raster/Niche_E/%s_%s%s.tif", groups[i], GCM, rcp))
        p<-as_tibble(rasterToPoints(richness))
        colnames(p)[3]<-"richness_future"
        
        p<-full_join(p, p_present, by=c("x", "y"))
        p$GCM<-GCM
        p$rcp<-rcp
        p$group<-groups[i]
        if (is.null(all_p)){
          all_p<-p
        }else{
          all_p<-bind_rows(all_p, p)
        }
      }
    }
  }
}
dim(all_p)
saveRDS(all_p, "../Object/future_richness.rda")

head(all_p)
tail(all_p)
all_p<-readRDS("../Object/future_richness.rda")

all_p[is.na(all_p$richness_future), "richness_future"]<-0
all_p[is.na(all_p$richness_present), "richness_present"]<-0
all_p[!complete.cases(all_p),]

all_p_se<-all_p%>%dplyr::group_by(GCM, rcp, y, group)%>%
  dplyr::summarise(mean_richness_future=mean(richness_future),
                   mean_richness_present=mean(richness_present))

if (F){
  n <- length(GCMs)-1
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors<-c(sample(col_vector, n), "black")
  saveRDS(colors, "../Object/colors.rda")
}
colors<-readRDS("../Object/colors.rda")

p_richness<-ggplot(all_p_se, aes(x=y, color=factor(GCM)))+
  geom_line(aes(y=mean_richness_future))+
  geom_line(aes(y=mean_richness_present), color="black")+
  theme_bw()+
  #scale_color_manual(name = "GCM", labels = c(GCMs[-1], "Present"), values=colors)+
  facet_wrap(~rcp+group, ncol=4, scales="free")

ggsave(p_richness, file="../Figures/future_richness.png", width=10, height = 6)                


all_p_se<-all_p%>%dplyr::group_by(rcp, x, y, group)%>%
  dplyr::summarise(mean_richness_future=mean(richness_future),
                   mean_richness_present=mean(richness_present))

all_p_se$differ<-all_p_se$mean_richness_future-all_p_se$mean_richness_present

mask<-raster("../Raster/mask.tif")
for (i in c(1:length(groups))){
    for (rcp in rcps){
      print(paste(groups[i], rcp))
      item<-all_p_se %>% dplyr::filter((rcp==rcp)&(group==groups[i]))
      r_richness_differ<-rasterFromXYZ(item[, c("x", "y", "differ")], crs=crs(mask), res=res(mask))
      NAvalue(r_richness_differ)<--9999
      writeRaster(r_richness_differ, NAflag=-9999,
                  sprintf("../Raster/Richness_differ/%s_%s.tif", rcp, groups[i]), 
                  overwrite=T)
    }
}




