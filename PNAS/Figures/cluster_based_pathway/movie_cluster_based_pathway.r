library(rgl)
library(ceramic)
#library(anglr)
library(ggnewscale)
library(ggplot2)
library(raster)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
#alt<-raster("../../Raster/ALT/alt_eck4_high_res.tif")
mask<-raster("../../Raster/mask_10km.tif")
mask_p<-data.frame(rasterToPoints(mask))
source("commonFuns/colors.r")
#p2<-data.frame(rasterToPoints(alt))
#p2[which(p2$x<=-12103059), "alt_eck4_high_res"]<-NA
#p2[which((p2$x>12912000)&(p2$y>5000000)), "alt_eck4_high_res"]<-NA
#values(alt)[!is.na(values(alt))]<-p2$alt_eck4_high_res

p_bak<-ggplot() + 
  geom_tile(data = mask_p, aes(x = x, y = y), fill="grey50", alpha=0.2)+
  map_theme

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=9
df_sp_list<-list()
for (group in c("Birds", "Mammals")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list, fill=T)

#ttt<-2
#df_sp_list<-df_sp_list[area>ttt]
df_sp_list$sp<-gsub(" ", "_", df_sp_list$sp)
j=9
euc.dist <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2) ^ 2+(y1-y2)^2)
}
exposure=5
for (j in c(1:nrow(layer_df))){
  for (exposure in c(0, 5)){
    persent<-0.2
    #for (j in c(3, 9)){
    #  for (exposure in c(5)){
    layer_item<-layer_df[j,]
    for (group in c("Birds", "Mammals")){
      if (group=="Birds"){
        persent<-0.2
      }
      if (group=="Mammals"){
        persent<-0.5
      }
      
      target_rda<-sprintf("../../Objects/cluster_based_pathway/merged/%s_%s_exposure_%d_sub_%d.rda",
                          group, layer_item$LABEL, exposure, persent * 100)
      print(target_rda)
      smooth_path<-readRDS(target_rda)
      
      yy=2050
      for (yy in c(2021:2100)){
        print(paste(exposure, j, nrow(layer_df), yy))
        smooth_path_item<-smooth_path[YEAR<=yy]
        smooth_path_item[smooth_path_item[, .I[YEAR==max(YEAR)], by=list(group, sp, line_group)]$V1]
        setDT(smooth_path_item)[, MAX_YEAR := max(YEAR), by = list(group, sp, line_group)]
        
        smooth_path_item$is_head<-ifelse(smooth_path_item$MAX_YEAR==smooth_path_item$YEAR, T, F)
        
        
        p<-p_bak+geom_path(data=smooth_path_item, aes(x=x, y=y, alpha=alpha, color=group,
                                                      group=line_group))+
          geom_point(data=smooth_path_item[smooth_path_item$is_head,], aes(x=x, y=y, color=group), size=0.05)+
          scale_alpha_continuous()+
          scale_color_manual(values = color_groups)+
          scale_fill_manual(values = color_groups)
        
        width<-13
        height<-6
        folder<-sprintf("../../Figures/cluster_based_pathway_movies/Movies/exposure_%d/%s/%s/per_%d", 
                        exposure, group, layer_item$LABEL, round(persent * 100))
        dir.create(folder, showWarnings = F, recursive = T)
        ggsave(p, filename=sprintf("%s/%d.png", 
                                   folder, yy), width=width, height = height)
        
        
      }
    }
  }
}

library(magick)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
exposure<-5
label<-"UKESM1_SSP585"
g<-"Birds"
year<-2021
persent<-0.2
for (label in c("UKESM1_SSP119", "UKESM1_SSP245", "UKESM1_SSP585")){
  for (year in c(2021:2100)){
    print(year)
    png(sprintf("../../Figures/cluster_based_pathway_movies/Movies/exposure_%d/Combined/%s/%d.png", 
                exposure, label, year), 
        width=1920, height=540, units = "px")
    par(mfrow=c(1,2),
        oma = c(0,0,0,0),
        mar = c(0,0,0,0),
        mgp = c(0, 0, 0),    # axis label at 2 rows distance, tick labels at 1 row
        xpd = NA,
        tcl=-1
    )
    for (g in c("Birds", "Mammals")){
      persent<-ifelse(g=="Birds", 0.2, 0.5)
      folder<-sprintf("../../Figures/cluster_based_pathway_movies/Movies/exposure_%d/%s/%s/per_%d", 
                      exposure, g, label, round(persent * 100))
      
      
      path<-image_read(sprintf("%s/%d.png", folder, year))
      path<-image_crop(path, "3000x1650+600+50")
      
      path<-image_annotate(path, g, gravity = "south", 
                           size = 120,  color = "#000000",
                           strokecolor = NULL, boxcolor = NULL)
      #path<-image_resize(path, "960x580")
      
      plot(path)
      
    }
    text(0,100,year,cex=5,font=2)
    
    dev.off()
  }
}

cd /media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/cluster_based_pathway_movies/Movies/exposure_5/Combined/UKESM1_SSP119
ffmpeg -r 5 -start_number 2021 -i %04d.png -y ../../../UKESM1_SSP119_exposure_5.mp4

cd /media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/cluster_based_pathway_movies/Movies/exposure_5/Combined/UKESM1_SSP245
ffmpeg -r 5 -start_number 2021 -i %04d.png -y ../../../UKESM1_SSP245_exposure_5.mp4

cd /media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Figures/cluster_based_pathway_movies/Movies/exposure_5/Combined/UKESM1_SSP585
ffmpeg -r 5 -start_number 2021 -i %04d.png -y ../../../UKESM1_SSP585_exposure_5.mp4

