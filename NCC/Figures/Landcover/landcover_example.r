library(raster)
library(ggplot2)
library(data.table)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask_10km<-raster("../../Raster/mask_10km.tif")
no_na_mask_10km<-!is.na(values(mask_10km))
mask_points_10km<-data.table(rasterToPoints(mask_10km))
r_sinu<-raster("../../Raster/2020_LC_Type5_sinu_1km.tif")
bi<-"Baryphthengus_ruficapillus"
disp_2100<-readRDS(sprintf("../../Objects/Dispersal/Birds/%s/MRI-ESM2-0_SSP585_0_dispersal_1_10km.rda", bi))
disp_2100<-disp_2100[["2100"]]
disp_2020<-readRDS(sprintf("../../Objects/Dispersal/Birds/%s/initial_disp_10km_exposure_0_dispersal_1.rda", bi))                   

extend<-c(min(disp_2100$x, disp_2020$x)-50000, 
          max(disp_2100$x, disp_2020$x)+50000, 
          min(disp_2100$y, disp_2020$y)-50000, 
          max(disp_2100$y, disp_2020$y)+50000)
mask<-mask_10km
mask_p<-mask_points_10km
mask_p$v<-NA
mask_p[mask_10km %in% disp_2020$mask_10km]$v<-1
values(mask)[no_na_mask_10km]<-mask_p$v
mask<-crop(mask, c(min(disp_2020$x - 5000), max(disp_2020$x + 5000), 
                   min(disp_2020$y - 5000), max(disp_2020$y + 5000)))

target_folder<-sprintf("../../Objects/Dispersal/%s/%s", "Birds", bi)
exposure_threshold<-0
dispersal<-1
mask_2020_file<-sprintf("%s/exposure_%d_dispersal_%d_disp_2020.tif", 
                        target_folder, exposure_threshold, dispersal)
mask_2020_file_1km<-sprintf("%s/exposure_%d_dispersal_%d_disp_2020_1km.tif", 
                            target_folder, exposure_threshold, dispersal)
writeRaster(mask, mask_2020_file, datatype="INT1U",
            overwrite=T)
gdalwarp(mask_2020_file, mask_2020_file_1km, tr=c(1000, 1000), r="near", ot="Byte", 
         srcnodata=255, dstnodata=255, overwrite=T)

mask_1km<-raster(mask_2020_file_1km)
lc_2020<-crop(r_sinu, extent(mask_1km))
extent(lc_2020)<-extent(mask_1km)
lc_2020<-mask(lc_2020, mask_1km)

plot(mask_10km)
plot(mask, add=T, col="black")
plot(lc_2020)
disp_2100<-disp_2100[suitable==1]
mask<-mask_10km
mask_p<-mask_points_10km
mask_p$v<-NA
mask_p[mask_10km %in% disp_2100$mask_10km]$v<-1
values(mask)[no_na_mask_10km]<-mask_p$v
mask<-crop(mask, c(min(disp_2100$x)-5000, max(disp_2100$x)+5000, 
                   min(disp_2100$y)-5000, max(disp_2100$y)+5000))
mask_2100_file<-sprintf("%s/%s_exposure_%d_dispersal_%d_disp_2100.tif", 
                        target_folder, item_str, exposure_threshold, dispersal)
mask_2100_file_1km<-sprintf("%s/%s_exposure_%d_dispersal_%d_disp_2100_1km.tif", 
                            target_folder, item_str, exposure_threshold, dispersal)
writeRaster(mask, mask_2100_file,datatype="INT1U",
            overwrite=T)
gdalwarp(mask_2100_file, mask_2100_file_1km, tr=c(1000, 1000), r= "near", ot="Byte", 
         srcnodata=255, dstnodata=255, overwrite=T)

mask_1km<-raster(mask_2100_file_1km)
lc_2100<-crop(r_sinu, extent(mask_1km))
extent(lc_2100)<-extent(mask_1km)
lc_2100<-mask(lc_2100, mask_1km)
source("commonFuns/colors.r")
mask_map<-crop(mask_10km, extend)
mask_map_p<-data.frame(rasterToPoints(mask_map))
mask_2020<-raster(mask_2020_file)
mask_2020_p<-data.frame(rasterToPoints(mask_2020))
p1<-ggplot()+
  geom_tile(data=mask_map_p, aes(x=x, y=y), fill=mask_color)+
  coord_equal()+
  geom_tile(data=mask_2020_p, aes(x=x, y=y), fill=colors_red[6])+
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
    panel.border = element_blank(),
    legend.position="none"
  )
ggsave(p1, filename="../../Figures/Landcover/Example/disp_2020.png")
lc_map<-crop(r_sinu, extend)
lc_map_p<-data.frame(rasterToPoints(lc_map))
lc_types<-unique(c(values(lc_2020), values(lc_2100)))
lc_types<-lc_types[!is.na(lc_types)]
lc_types<-lc_types[lc_types!=0]
n_colors<-length(lc_types)
colors<-rainbow(n_colors)
names(colors)<-lc_types
lc_map_p<-lc_map_p[which(lc_map_p$X2020_LC_Type5_sinu_1km!=0),]
p2<-ggplot()+
  geom_tile(data=mask_map_p, aes(x=x, y=y), fill=mask_color)+
  coord_equal()+
  geom_tile(data=lc_map_p, aes(x=x, y=y, fill=factor(X2020_LC_Type5_sinu_1km)))+
  scale_fill_manual(breaks=lc_types, values=colors)+
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
    panel.border = element_blank(),
    legend.position="none"
  )
ggsave(p2, filename="../../Figures/Landcover/Example/lc_map.png")
lc_2020_p<-data.frame(rasterToPoints(lc_2020))
lc_2020_p<-lc_2020_p[which(lc_2020_p$X2020_LC_Type5_sinu_1km!=0),]
p3<-ggplot()+
  geom_tile(data=mask_map_p, aes(x=x, y=y), fill=mask_color)+
  coord_equal()+
  geom_tile(data=lc_2020_p, aes(x=x, y=y, fill=factor(X2020_LC_Type5_sinu_1km)))+
  scale_fill_manual(breaks=lc_types, values=colors)+
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
    panel.border = element_blank(),
    legend.position="none"
  )
ggsave(p3, filename="../../Figures/Landcover/Example/lc_2020.png")

#0 At least 60% of area is covered by permanent water bodies.
#Evergreen Needleleaf Trees 1 Dominated by evergreen conifer trees (>2m). Tree cover >10%.
#Evergreen Broadleaf Trees 2 Dominated by evergreen broadleaf and palmate trees (>2m). Tree cover >10%.
#Deciduous Needleleaf Trees 3 Dominated by deciduous needleleaf (larch) trees (>2m). Tree cover >10%.
#Deciduous Broadleaf Trees 4 Dominated by deciduous broadleaf trees (>2m). Tree cover >10%.
#Shrub 5 Shrub (1-2m) cover >10%.
#Grass 6 Dominated by herbaceous annuals (<2m) that are not cultivated.
#Cereal Croplands 7 Dominated by herbaceous annuals (<2m). At least 60% cultivated cereal crops.
#Broadleaf Croplands 8 Dominated by herbaceous annuals (<2m). At least 60% cultivated broadleaf crops.
#Urban and Built-up Lands 9 At least 30% impervious surface area including building materials, asphalt, and vehicles.
#Permanent Snow and Ice 10 At least 60% of area is covered by snow and ice for at least 10 months of the year.
#Barren 11 At least 60% of area is non-vegetated barren (sand, rock, soil) with less than 10% vegetation.

lc_type_label<-c("Water Bodies"=0, "Evergreen Needleleaf Trees"=1,
                 "Evergreen Broadleaf Trees"=2, "Deciduous Needleleaf Trees"=3,
                 "Deciduous Broadleaf Trees"=4, "Shrub"=5,
                 "Grass"=6, "Cereal Croplands"=7, "Broadleaf Croplands"=8,
                 "Urban and Built-up Lands"=9, "Permanent Snow and Ice"=10,
                 "Barren"=11)
lc_type_label<-data.table(lc_v=lc_type_label, lc_label=names(lc_type_label))
lc_2020_p<-data.table(lc_2020_p)
lc_2020_p_e<-lc_2020_p[, .(N=.N), by=list(X2020_LC_Type5_sinu_1km)]
lc_2020_p_e<-merge(lc_2020_p_e, lc_type_label, by.x="X2020_LC_Type5_sinu_1km", by.y="lc_v")
p4<-ggplot(lc_2020_p_e)+geom_bar(aes(x=lc_label, y=N, fill=factor(X2020_LC_Type5_sinu_1km)), stat = "identity")+
  scale_fill_manual(breaks=lc_types, values=colors)+
  coord_flip()+
  theme_bw()+
  labs(y="Area", x="")+
  theme(legend.position = "none")
p4
ggsave(p4, filename="../../Figures/Landcover/Example/lc_2020_bar.png")

mask_2100<-raster(mask_2100_file)
mask_2100_p<-data.frame(rasterToPoints(mask_2100))
p5<-ggplot()+
  geom_tile(data=mask_map_p, aes(x=x, y=y), fill=mask_color)+
  coord_equal()+
  geom_tile(data=mask_2100_p, aes(x=x, y=y), fill=colors_red[6])+
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
    panel.border = element_blank(),
    legend.position="none"
  )
ggsave(p5, filename="../../Figures/Landcover/Example/disp_2100.png")

lc_2100_p<-data.frame(rasterToPoints(lc_2100))
lc_2100_p<-lc_2100_p[which(lc_2100_p$X2020_LC_Type5_sinu_1km!=0),]
p6<-ggplot()+
  geom_tile(data=mask_map_p, aes(x=x, y=y), fill=mask_color)+
  coord_equal()+
  geom_tile(data=lc_2100_p, aes(x=x, y=y, fill=factor(X2020_LC_Type5_sinu_1km)))+
  scale_fill_manual(breaks=lc_types, values=colors)+
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
    panel.border = element_blank(),
    legend.position="none"
  )
ggsave(p6, filename="../../Figures/Landcover/Example/lc_2100.png")

#0 At least 60% of area is covered by permanent water bodies.
#Evergreen Needleleaf Trees 1 Dominated by evergreen conifer trees (>2m). Tree cover >10%.
#Evergreen Broadleaf Trees 2 Dominated by evergreen broadleaf and palmate trees (>2m). Tree cover >10%.
#Deciduous Needleleaf Trees 3 Dominated by deciduous needleleaf (larch) trees (>2m). Tree cover >10%.
#Deciduous Broadleaf Trees 4 Dominated by deciduous broadleaf trees (>2m). Tree cover >10%.
#Shrub 5 Shrub (1-2m) cover >10%.
#Grass 6 Dominated by herbaceous annuals (<2m) that are not cultivated.
#Cereal Croplands 7 Dominated by herbaceous annuals (<2m). At least 60% cultivated cereal crops.
#Broadleaf Croplands 8 Dominated by herbaceous annuals (<2m). At least 60% cultivated broadleaf crops.
#Urban and Built-up Lands 9 At least 30% impervious surface area including building materials, asphalt, and vehicles.
#Permanent Snow and Ice 10 At least 60% of area is covered by snow and ice for at least 10 months of the year.
#Barren 11 At least 60% of area is non-vegetated barren (sand, rock, soil) with less than 10% vegetation.

lc_2100_p<-data.table(lc_2100_p)
lc_2100_p_e<-lc_2100_p[, .(N=.N), by=list(X2020_LC_Type5_sinu_1km)]
lc_2100_p_e<-merge(lc_2100_p_e, lc_type_label, by.x="X2020_LC_Type5_sinu_1km", by.y="lc_v")
p7<-ggplot(lc_2100_p_e)+geom_bar(aes(x=lc_label, y=N, fill=factor(X2020_LC_Type5_sinu_1km)), stat = "identity")+
  scale_fill_manual(breaks=lc_types, values=colors)+
  coord_flip()+
  theme_bw()+
  labs(y="Area", x="")+
  theme(legend.position = "none")
p7
ggsave(p7, filename="../../Figures/Landcover/Example/lc_2100_bar.png")

pp<-list(p1, p2, p3, p4, p5, p2, p6, p7)
ppp<-ggarrange(plotlist=pp, nrow=2, ncol=4, widths=c(1,1,1,2), 
               labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"))
ggsave(ppp, filename="../../Figures/Landcover/Example/lc_all.png", width=12, height=7)
