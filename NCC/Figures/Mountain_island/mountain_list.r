library(sf)
library(raster)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_10km.tif")
if (F){
  mountain<-st_read("../../Shape/GMBA mountain inventory V1.2/GMBA Mountain Inventory_v1.2-World.shp")
  mountain<-st_transform(mountain, crs=st_crs(mask))
  all_info<-readRDS("../../Objects/Island/species_mountain_10km.rda")
  all_info$mountain_p<-all_info$mountain/(all_info$mountain+all_info$plain)
  nrow(all_info[mountain_p>0.95])
  mountain_species<-all_info[mountain_p>0.95]
  mountain_species[(plain+mountain)==0]
  i=1
  mountain_list<-list()
  for (i in c(1:nrow(mountain_species))){
    print(paste(i, nrow(mountain_species)))
    item<-mountain_species[i]
    target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", item$group, item$sp)
    if (!file.exists(target)){
      next()
    }
    disp<-readRDS(target)
    f<-st_intersects(disp, mountain)
    for (j in c(1:length(f))){
      if (length(f[[j]])>0){
        ff<-mountain[j,]
        if (F){
          plot(st_geometry(ff))
          plot(st_geometry(disp), add=T, col="red")
        }
        df_item<-data.frame(Name=ff$Name, Country=ff$Country, sp=item$sp, group=item$group)
        mountain_list[[length(mountain_list)+1]]<-df_item
      }
    }
    
  }
  mountain_list<-rbindlist(mountain_list)
  saveRDS(mountain_list, "../../Figures/Mountain_island/mountain_species.rda")
  write.csv(mountain_list, "../../Figures/Mountain_island/mountain_species.csv", row.names = F)
  
  PRESENCE<-c(1,2,3,4,5)
  ORIGIN<-c(1,2,3,5,6)
  SEASONAL<-c(1,2)
  sp_list<-mountain_list[, .(N=.N), by=list(sp, group)]
  i=1
  species_richness<-mask
  values(species_richness)[!is.na(values(species_richness))]<-0
  for (i in c(1:nrow(sp_list))){
    print(paste(i, nrow(sp_list)))
    item<-sp_list[i]
    target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", item$group, item$sp)
    if (!file.exists(target)){
      next()
    }
    dis<-readRDS(target)
    if (item$group=="Birds"){
      dis<-dis[which((dis$PRESENCE %in% PRESENCE)&
                       (dis$ORIGIN %in% ORIGIN)&
                       (dis$SEASONAL %in% SEASONAL)),]
    }else{
      dis<-dis[which((dis$presence %in% PRESENCE)&
                       (dis$origin %in% ORIGIN)&
                       (dis$seasonal %in% SEASONAL)),]
    }
    if (nrow(dis)==0){
      next()
    }
    x<-mask(mask, dis)
    values(x)[!is.na(values(x))]<-1
    values(x)[is.na(values(x))]<-0
    species_richness<-species_richness+x
  }
  writeRaster(species_richness, "../../Figures/Mountain_island/mountain_species.tif")
}

species_richness<-raster("../../Figures/Mountain_island/mountain_species.tif")
plot(mask)
plot(species_richness)

points<-data.frame(rasterToPoints(species_richness))
library(ggplot2)
source("commonFuns/colors.r")
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
colors<-myPalette(max(points$mountain_species))
p<-ggplot(points)+geom_tile(aes(x=x, y=y, fill=mountain_species))+
  scale_fill_gradientn(colors=colors)+
  coord_fixed()+
  labs(fill="Species richness")+
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
ggsave(p, filename="../../Figures/Mountain_island/mountain_species.png")

mountain_list<-readRDS("../../Figures/Mountain_island/mountain_species.rda")
sp_list<-mountain_list[, .(N=.N), by=list(sp, group)]
GCMs<-c("UKESM1", "EC-Earth3-Veg", "MRI-ESM2-0")
SSPs<-c("SSP245", "SSP585", "SSP119")
coms<-expand.grid(GCM=GCMs, SSP=SSPs, dispersal=c(0, 1), exposure=c(0, 5),
                  stringsAsFactors = F)
i=1
j=1
all_df<-list()
for (i in c(1:nrow(sp_list))){
  print(paste(i, nrow(sp_list)))
  item<-sp_list[i]
  for (j in c(1:nrow(coms))){
    com<-coms[j,]
    source<-sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d_10km.rda",
                    item$group, item$sp, com$GCM, com$SSP, com$exposure, com$dispersal)
    res<-"10km"
    if (!file.exists(source)){
      source<-sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d.rda",
                      item$group, item$sp, com$GCM, com$SSP, com$exposure, com$dispersal)
      res<-"100km"
    }
    if (!file.exists(source)){
      next()
    }
    df<-readRDS(source)
    N_2020<-nrow(df[["2021"]])
    N_2100<-nrow(df[["2100"]])
    if (is.null(N_2100)){
      N_2100<-0
    }
    com$sp<-item$sp
    com$group<-item$group
    com$N2020<-N_2020
    com$N2100<-N_2100
    all_df[[length(all_df)+1]]<-com
  }
}
all_df<-rbindlist(all_df, fill=T)
saveRDS(all_df, "../../Figures/Mountain_island/mountain_species_stat.rda")
all_df_extincted<-all_df[N2100==0]

all_df_extincted
mountain_list<-readRDS("../../Figures/Mountain_island/mountain_species.rda")
mountain_list_SSP119<-mountain_list[sp %in% all_df_extincted[SSP=="SSP119"]$sp]
write.csv(mountain_list_SSP119, "../../Figures/Mountain_island/mountain_species_SSP119.csv", row.names = F)
sp_list<-mountain_list_SSP119[, .(N=.N), by=list(sp, group)]
i=1
species_richness<-mask
values(species_richness)[!is.na(values(species_richness))]<-0
for (i in c(1:nrow(sp_list))){
  print(paste(i, nrow(sp_list)))
  item<-sp_list[i]
  target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", item$group, item$sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  
  if (nrow(dis)==0){
    next()
  }
  x<-mask(mask, dis)
  values(x)[!is.na(values(x))]<-1
  values(x)[is.na(values(x))]<-0
  species_richness<-species_richness+x
}
writeRaster(species_richness, "../../Figures/Mountain_island/mountain_species_SSP119.tif")

mountain_list_SSP245<-mountain_list[sp %in% all_df_extincted[SSP=="SSP245"]$sp]
write.csv(mountain_list_SSP245, "../../Figures/Mountain_island/mountain_species_SSP245.csv", row.names = F)

sp_list<-mountain_list_SSP245[, .(N=.N), by=list(sp, group)]
i=1
species_richness<-mask
values(species_richness)[!is.na(values(species_richness))]<-0
for (i in c(1:nrow(sp_list))){
  print(paste(i, nrow(sp_list)))
  item<-sp_list[i]
  target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", item$group, item$sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  
  if (nrow(dis)==0){
    next()
  }
  x<-mask(mask, dis)
  values(x)[!is.na(values(x))]<-1
  values(x)[is.na(values(x))]<-0
  species_richness<-species_richness+x
}
writeRaster(species_richness, "../../Figures/Mountain_island/mountain_species_SSP245.tif")

mountain_list_SSP585<-mountain_list[sp %in% all_df_extincted[SSP=="SSP585"]$sp]
write.csv(mountain_list_SSP585, "../../Figures/Mountain_island/mountain_species_SSP585.csv", row.names = F)


sp_list<-mountain_list_SSP585[, .(N=.N), by=list(sp, group)]
i=1
species_richness<-mask
values(species_richness)[!is.na(values(species_richness))]<-0
for (i in c(1:nrow(sp_list))){
  print(paste(i, nrow(sp_list)))
  item<-sp_list[i]
  target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", item$group, item$sp)
  if (!file.exists(target)){
    next()
  }
  dis<-readRDS(target)
  
  if (nrow(dis)==0){
    next()
  }
  x<-mask(mask, dis)
  values(x)[!is.na(values(x))]<-1
  values(x)[is.na(values(x))]<-0
  species_richness<-species_richness+x
}
writeRaster(species_richness, "../../Figures/Mountain_island/mountain_species_SSP585.tif")

species_richness<-raster("../../Figures/Mountain_island/mountain_species_SSP119.tif")
plot(mask)
plot(species_richness)

points<-data.frame(rasterToPoints(species_richness))
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
colors<-myPalette(max(points$mountain_species_SSP119))
p<-ggplot(points)+geom_tile(aes(x=x, y=y, fill=mountain_species_SSP119))+
  scale_fill_gradientn(colors=colors)+
  coord_fixed()+
  labs(fill="Species richness")+
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
#p
ggsave(p, filename="../../Figures/Mountain_island/mountain_species_SSP119.png")



species_richness<-raster("../../Figures/Mountain_island/mountain_species_SSP245.tif")
plot(mask)
plot(species_richness)

points<-data.frame(rasterToPoints(species_richness))
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
colors<-myPalette(max(points$mountain_species_SSP245))
p<-ggplot(points)+geom_tile(aes(x=x, y=y, fill=mountain_species_SSP245))+
  scale_fill_gradientn(colors=colors)+
  coord_fixed()+
  labs(fill="Species richness")+
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
#p
ggsave(p, filename="../../Figures/Mountain_island/mountain_species_SSP245.png")


species_richness<-raster("../../Figures/Mountain_island/mountain_species_SSP585.tif")
plot(mask)
plot(species_richness)

points<-data.frame(rasterToPoints(species_richness))
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
colors<-myPalette(max(points$mountain_species_SSP585))
p<-ggplot(points)+geom_tile(aes(x=x, y=y, fill=mountain_species_SSP585))+
  scale_fill_gradientn(colors=colors)+
  coord_fixed()+
  labs(fill="Species richness")+
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
#p
ggsave(p, filename="../../Figures/Mountain_island/mountain_species_SSP585.png")

