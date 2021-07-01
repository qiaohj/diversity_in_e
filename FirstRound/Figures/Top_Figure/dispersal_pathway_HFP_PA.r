library(raster)
library(dplyr)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
mask<-raster("../../Raster/mask_index.tif")
source("commonFuns/colors.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
SSP_i<-SSPs[1]
threshold_i=1
HFP<-raster("../../Raster/HFP2009_Low.tif")
PA_International<-raster("../../Raster/PA/International.tif")
if (F){
  f<-"China_PAs"
  for (f in c("China_PAs", "International", "National", "Not Applicable", "Regional")){
    r_pa<-raster(sprintf("../../Raster/PA/moll/%s.tif", f))
    r_re<-projectRaster(r_pa, mask, method="ngb")
    writeRaster(r_re, sprintf("../../Raster/PA/%s.tif", f), overwrite=T)  
  }
  
}
threshold_HFP<-quantile(values(HFP), c(0.1, 0.9), na.rm=T)
mask_p<-data.frame(rasterToPoints(mask))
df_all<-NULL
for (SSP_i in SSPs){
  for (threshold_i in c(5)){
    print(paste(SSP_i, threshold_i))
    r<-raster(sprintf("../../Figures_Full_species/Top_Figure_all/TIF/%s_%d.tif", 
                        SSP_i, threshold_i))
    r_p<-data.frame(rasterToPoints(r))
    colnames(r_p)[3]<-"N_Pathway"
    r_p$HFP<-raster::extract(HFP, r_p[, c("x", "y")])
    for (f in c("China_PAs", "International", "National", "Not Applicable", "Regional")){
      r_pa<-raster(sprintf("../../Raster/PA/%s.tif", f))
      r_p[, f]<-raster::extract(r_pa, r_p[, c("x", "y")])
    }
    r_p[is.na(r_p)]<-0
    r_p$PA<-r_p$China_PAs+r_p$International+r_p$National+r_p$`Not Applicable`+r_p$Regional
    r_p$PA[which(r_p$PA>1)]<-1
   
    threshold<-quantile(r_p$N_Pathway, 0.75)
    r_p<-r_p[which(r_p$N_Pathway>=threshold),]
    r_p$HFP_T<-""
    r_p[which(r_p$HFP<=threshold_HFP[1]), "HFP_T"]<-"Wilderness"
    r_p[which(r_p$HFP>=threshold_HFP[2]), "HFP_T"]<-"High HFP"
    table(r_p$HFP_T)
    #r_p<-r_p[which(r_p$HFP_T!=""),]
    r_p$exposure<-ifelse(threshold_i==1, " no exposure", "5-year exposure")
    r_p$SSP<-SSP_i
    r_p$PA<-ifelse(r_p$PA==0, "Unprotected area", "Protected area")
    df_all<-bind_dplyr(df_all, r_p)
  }
}
    p<-ggplot()+
      geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
      geom_tile(data=df_all, aes(x=x, y=y, fill=PA))+
      geom_point(data=df_all%>%dplyr::filter(HFP_T=="High HFP"), aes(x=x, y=y), 
                 size=0.01, shape=4)+
      facet_wrap(~SSP, ncol=1)+
      scale_fill_manual(values=c(colors_blue[5], colors_red[5]))+
      labs(fill = "")+
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = map_background, color = NA), 
        panel.background = element_blank(), 
        legend.background = element_rect(fill = map_background, color = NA),
        panel.border = element_blank()
      )
    p
ggsave(p, filename="../../Figures_Full_species/Pathway_HFP_PA/Pathway_HFP_PA.pdf",
       width=10, height=20)    
  