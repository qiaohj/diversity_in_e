library(raster)
library(dplyr)
library(RColorBrewer)
library(Rmisc)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("functions.r")
SSPs<-c("SSP119", "SSP245", "SSP585")
GCMs<-c("EC-Earth3-Veg", "UKESM1", "MRI-ESM2-0")
dispersals<-data.frame(M=c(0:5), N=1)
years<-c(2014:2100)
args = commandArgs(trailingOnly=TRUE)
f<-args[1]
g<-args[2]
if (is.na(f)){
  f<-"Diversity"
}
if (is.na(g)){
  g<-"Amphibians"
}
SSP_i<-SSPs[3]
GCM_i<-GCMs[1]
dis_i<-2
y=2014
y=2100
slope<-raster("../../Raster/ALT/slope_eck4.tif")
alt<-raster("../../Raster/ALT/alt_eck4.tif")
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))


loss <- colorRampPalette(c("white","blue"))
gain <- colorRampPalette(c("white","red"))
result<-NULL
#for (f in c("Diversity", "Diversity_with_human")){
#  for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
    for (SSP_i in SSPs){
      for (dis_i in c(1:nrow(dispersals))){
        dis<-dispersals[dis_i,]
        for (y in years){
          item<-data.frame(type=f, group=g, SSP=SSP_i, M=dis$M, N=dis$N, year=y)
          print(paste(f, g, SSP_i, dis$M, dis$N, y))
          points<-NULL
          for (GCM_i in GCMs){
            #print(GCM_i)
            r<-raster(sprintf("../../Figures/%s/%s/species.richness/%s_%s_%d_%d/%d.tif", f, g, GCM_i, SSP_i, dis$M, dis$N, y))
            if (y==2014){
              r_raw<-r
              p_raw<-data.frame(rasterToPoints(r_raw))
              colnames(p_raw)<-c("x", "y", "v_raw")
              p_raw$slope<-extract(slope, p_raw[, c("x", "y")])
              p_raw$alt<-extract(alt, p_raw[, c("x", "y")])
            }else{
              p_r<-data.frame(rasterToPoints(r))
              colnames(p_r)<-c("x", "y", "v")
              p_r<-left_join(p_raw, p_r, by=c("x", "y"))
              p_r[is.na(p_r)]<-0
              p_r$gain_loss<-p_r$v-p_r$v_raw
              p_r$gain_loss_ratio<-p_r$gain_loss/p_r$v_raw
              p_r$GCM<-GCM_i
              points<-bind(points, p_r)
              if (F){
                r_temp<-mask
                p_temp<-mask_p
                p_temp<-left_join(p_temp, p_r, by=c("x", "y"))
                values(r_temp)[!is.na(values(r_temp))]<-p_temp$gain_loss
                r_temp_1<-r_temp
                values(r_temp_1)[values(r_temp_1)>0]<-NA
                r_temp_2<-r_temp
                values(r_temp_2)[values(r_temp_2)<=0]<-NA
                plot(r_temp_1, col=loss(100))
                plot(r_temp_2, col=gain(100), add=T)
                values(r_temp)[!is.na(values(r_temp))]<-p_temp$gain_loss_ratio
                plot(r_temp)
              }
            }
            
          }
          if (y==2014){
            next()
          }
          points_se<-points%>%dplyr::group_by(x, y, v_raw, slope, alt)%>%
            dplyr::summarise(mean_v=mean(v),
                             mean_gain_loss=mean(gain_loss),
                             sd_v=sd(v),
                             sd_gain_loss=sd(gain_loss),
                             CI_v=CI(v)[2]-CI(v)[3],
                             CI_gain_loss=CI(gain_loss)[2]-CI(gain_loss)[3])
          points_se$type<-item$type
          points_se$group<-item$group
          points_se$SSP<-item$SSP
          points_se$M<-item$M
          points_se$N<-item$N
          points_se$year<-item$year
          result<-bind(result, points_se)
        }
      }
    }
saveRDS(result, sprintf("../../Figures/Species_gain_loss/%s_%s.rda", f, g))
  #}
#}
