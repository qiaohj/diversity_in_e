library(raster)
#library(rgdal)
library(dplyr)
library(vegan)
library(data.table)
library(ggplot2)
library(tidyr)

#rm(list=ls())
setDTthreads(threads=1)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Mammals"
}

exposure<-as.numeric(args[2])
if (is.na(exposure)){
  exposure<-0
}

dispersal<-as.numeric(args[3])
if (is.na(dispersal)){
  dispersal<-0
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

sp_list<-readRDS(sprintf("../../Objects/rerun_10km_%s.rda", group))
sp_list<-sp_list[run_10km==T&run_100km==T]
sp_list<-sp_list[, .(N=.N), by=list(sp, group)]
i=1
j=1
k=2

mask_100km<-raster("../../Raster/mask_100km.tif")
points_100km<-as_tibble(rasterToPoints(mask_100km))
colnames(points_100km)[3]<-"mask"
#mask_10km<-raster("../../Raster/mask_10km.tif")

if (F){
  points_10km<-as_tibble(rasterToPoints(mask_10km))
  cols<-c("x", "y")
  points_10km$mask_100km<-raster::extract(mask_100km, points_10km[, c("x", "y")])
  points_10km<-data.table(points_10km)
  
  saveRDS(points_10km, "../../Raster/points_10km.rda")
}
#points_10km<-readRDS("../../Raster/points_10km.rda")
add_location<-function(indices, location, type){
  location$metric<-indices
  location$type<-type
  location
}


for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  
  target_folder<-sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d_10km_2_100km/%s/%s", 
                         exposure, dispersal, group, layer$LABEL)
  
  target<-sprintf("%s/indices_df.rda", target_folder)
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  #dir.create(sprintf("%s/species.richness", target_folder), showWarnings = F)
  #dir.create(sprintf("%s/species.richness/TIF", target_folder), showWarnings = F)
  #dir.create(sprintf("%s/species.richness/PNG", target_folder), showWarnings = F)
  i=247
  print(paste("READING DATA", target_folder))
  diversity_df_10km<-readRDS(sprintf("%s/diversity_df.rda", target_folder))
  diversity_df_100km<-readRDS(sprintf("%s/diversity_df.rda", 
                                      sprintf("../../Objects/Diversity_exposure_%d_dispersal_%d/%s/%s", 
                                              exposure, dispersal, group, layer$LABEL)))
  
  YYYY<-names(diversity_df_10km)[1]
  indices_df<-list()
  indices_100km_df<-list()
  indices_10km_df<-list()
  sp_dis<-NULL
  sp_dis_100km<-NULL
  sp_dis_10km<-NULL
  for (YYYY in names(diversity_df_10km)){
    print(YYYY)
    indices<-list()
    indices_10km<-list()
    indices_100km<-list()
    
    iii=1
    diversity_100km_only<-list()
    for (iii in c(1:length(diversity_df_100km[[YYYY]]))){
      i_item<-diversity_df_100km[[YYYY]][[iii]]
      if (!is.null(i_item)){
        sp<-i_item[1]$sp
        if (sp %in% sp_list$sp){
          diversity_100km_only[[sp]]<-i_item
        }
      }
    }
    
    diversity_10km_only<-list()
    iii=1
    for (sp in names(diversity_df_10km[[YYYY]])){
      if (sp %in% names(diversity_100km_only)){
        diversity_10km_only[[sp]]<-diversity_df_10km[[YYYY]][[sp]]
      }
    }
    
    
    
    if (F){
      names_10km<-names(diversity_10km_only)
      names_10km[names_10km=="Tragelaphus_eurycerus"]
      names_100km<-names(diversity_100km_only)
      names_10km[!(names_10km %in% names_100km)]
      names_100km[!(names_100km %in% names_10km)]
    }
    
    diversity<-diversity_df_10km[[YYYY]]
    diversity_10km<-diversity_10km_only
    diversity_100km<-diversity_100km_only
    print(paste("BINDING DATA", YYYY, target_folder))
    diversity<-rbindlist(diversity)
    diversity<-unique(diversity)
    diversity_10km<-rbindlist(diversity_10km)
    diversity_10km<-unique(diversity_10km)
    diversity_100km<-rbindlist(diversity_100km)
    diversity_100km<-unique(diversity_100km)
    print(paste("COUNTING DATA", YYYY, target_folder))
    colnames(diversity)[c(1:3)]<-c("x", "y", "mask")
    colnames(diversity_10km)[c(1:3)]<-c("x", "y", "mask")
    colnames(diversity_100km)[c(1:3)]<-c("x", "y", "mask")
    
    n_dis<-diversity%>%dplyr::group_by(YEAR, sp)%>%
      dplyr::summarise(N=n(),
                       xabsmin=min(abs(x)),
                       xabsmax=max(abs(x)),
                       xmin=min(x),
                       xmax=max(x),
                       xrange=range(x)[2]-range(x)[1],
                       yabsmin=min(abs(y)),
                       yabsmax=max(abs(y)),
                       ymin=min(y),
                       ymax=max(y),
                       yrange=range(y)[2]-range(y)[1],
                       xmean=mean(x),
                       ymean=mean(y))
    
    n_dis_10km<-diversity_10km%>%dplyr::group_by(YEAR, sp)%>%
      dplyr::summarise(N=n(),
                       xabsmin=min(abs(x)),
                       xabsmax=max(abs(x)),
                       xmin=min(x),
                       xmax=max(x),
                       xrange=range(x)[2]-range(x)[1],
                       yabsmin=min(abs(y)),
                       yabsmax=max(abs(y)),
                       ymin=min(y),
                       ymax=max(y),
                       yrange=range(y)[2]-range(y)[1],
                       xmean=mean(x),
                       ymean=mean(y))
    
    n_dis_100km<-diversity_100km%>%dplyr::group_by(YEAR, sp)%>%
      dplyr::summarise(N=n(),
                       xabsmin=min(abs(x)),
                       xabsmax=max(abs(x)),
                       xmin=min(x),
                       xmax=max(x),
                       xrange=range(x)[2]-range(x)[1],
                       yabsmin=min(abs(y)),
                       yabsmax=max(abs(y)),
                       ymin=min(y),
                       ymax=max(y),
                       yrange=range(y)[2]-range(y)[1],
                       xmean=mean(x),
                       ymean=mean(y))
    
    
    if (F){
      library(ggplot2)
      ggplot(n_dis)+geom_histogram(aes(x=N, fill=factor(res)))+xlim(0, 10000)
    }
    if (is.null(sp_dis)){
      sp_dis<-n_dis
      sp_dis_10km<-n_dis_10km
      sp_dis_100km<-n_dis_100km
    }else{
      sp_dis<-bind_rows(sp_dis, n_dis)
      sp_dis_10km<-bind_rows(sp_dis_10km, n_dis_10km)
      sp_dis_100km<-bind_rows(sp_dis_100km, n_dis_100km)
    }
    if (T){
      print(paste("EXTRACTING DATA", YYYY, target_folder))
      diversity$index<-raster::extract(mask_100km, diversity[, c("x", "y")])
      diversity<-diversity[complete.cases(diversity),]
      diversity$v<-1
      
      diversity_10km$index<-raster::extract(mask_100km, diversity_10km[, c("x", "y")])
      diversity_10km<-diversity_10km[complete.cases(diversity_10km),]
      diversity_10km$v<-1
      
      diversity_100km$index<-raster::extract(mask_100km, diversity_100km[, c("x", "y")])
      diversity_100km<-diversity_100km[complete.cases(diversity_100km),]
      diversity_100km$v<-1
      
      print(paste("Matrix DATA", YYYY, target_folder))
      t_m<-pivot_wider(diversity, id_cols=index, names_from = sp, values_from=v, values_fill=0)
      d2<-t_m[,-1]
      points_item<-left_join(t_m[,1], points_100km, by=c("index"="mask"))
      
      t_m_10km<-pivot_wider(diversity_10km, id_cols=index, names_from = sp, values_from=v, values_fill=0)
      d2_10km<-t_m_10km[,-1]
      points_item_10km<-left_join(t_m_10km[,1], points_100km, by=c("index"="mask"))
      
      t_m_100km<-pivot_wider(diversity_100km, id_cols=index, names_from = sp, values_from=v, values_fill=0)
      d2_100km<-t_m_100km[,-1]
      points_item_100km<-left_join(t_m_100km[,1], points_100km, by=c("index"="mask"))
      
      
      print(paste("simpson", YYYY, target_folder))
      simpson <- diversity(d2, "simpson")
      indices[["simpson"]]<-add_location(simpson, points_item, "simpson")
      simpson_10km <- diversity(d2_10km, "simpson")
      indices_10km[["simpson"]]<-add_location(simpson_10km, points_item_10km, "simpson")
      simpson_100km <- diversity(d2_100km, "simpson")
      indices_100km[["simpson"]]<-add_location(simpson_100km, points_item_100km, "simpson")
      
      print(paste("shannon", YYYY, target_folder))
      shannon<-diversity(d2, "shannon")
      indices[["shannon"]]<-add_location(shannon, points_item, "shannon")
      shannon_10km <- diversity(d2_10km, "shannon")
      indices_10km[["shannon"]]<-add_location(shannon_10km, points_item_10km, "shannon")
      shannon_100km <- diversity(d2_100km, "shannon")
      indices_100km[["shannon"]]<-add_location(shannon_100km, points_item_100km, "shannon")
      
      print(paste("inv simpson", YYYY, target_folder))
      invsimp <- diversity(d2, "inv")
      indices[["invsimp"]]<-add_location(invsimp, points_item, "invsimp")
      invsimp_10km <- diversity(d2_10km, "invsimp")
      indices_10km[["invsimp"]]<-add_location(invsimp_10km, points_item_10km, "invsimp")
      invsimp_100km <- diversity(d2_100km, "invsimp")
      indices_100km[["invsimp"]]<-add_location(invsimp_100km, points_item_100km, "invsimp")
      
      #print(paste("fisher.alpha", YYYY, target_folder))
      #alpha <- fisher.alpha(d2)
      #indices[["alpha"]]<-add_location(alpha, points_item, "alpha")
      print(paste("species.richness", YYYY, target_folder))
      species.richness <- specnumber(d2) # observed number of species
      indices[["species.richness"]]<-add_location(species.richness, points_item, "species.richness")
      species.richness_10km <- specnumber(d2_10km)
      indices_10km[["species.richness"]]<-add_location(species.richness_10km, points_item_10km, "species.richness")
      species.richness_100km <- specnumber(d2_100km)
      indices_100km[["species.richness"]]<-add_location(species.richness_100km, points_item_100km, "species.richness")
      
      #r_temp<-rasterFromXYZ(indices[["species.richness"]][, c("x", "y", "metric")], res=res(mask), crs=crs(mask))
      #writeRaster(r_temp, sprintf("%s/species.richness/TIF/%s.tif", target_folder, YYYY), overwrite=T)
      #png(filename=sprintf("%s/species.richness/PNG/%s.png", target_folder, YYYY), width=640, height=480)
      #plot(r_temp)
      #dev.off()
      print(paste("raremax", YYYY, target_folder))
      raremax <- min(rowSums(d2))
      print(paste("Srare", YYYY, target_folder))
      Srare <- rarefy(d2, raremax) 
      indices[["Srare"]]<-add_location(Srare, points_item, "Srare")
      
      raremax_10km <- min(rowSums(d2_10km))
      Srare_10km <- rarefy(d2_10km, raremax_10km) 
      indices_10km[["Srare"]]<-add_location(Srare_10km, points_item_10km, "Srare")
      
      raremax_100km <- min(rowSums(d2_100km))
      Srare_100km <- rarefy(d2_100km, raremax_100km) 
      indices_100km[["Srare"]]<-add_location(Srare_100km, points_item_100km, "Srare")
      
      ## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
      print(paste("unbias.simp", YYYY, target_folder))
      unbias.simp <- rarefy(d2, raremax) - 1
      indices[["unbias.simp"]]<-add_location(unbias.simp, points_item, "unbias.simp")
      print(paste("Pielous_evenness", YYYY, target_folder))
      Pielous_evenness <- shannon/log(species.richness)
      indices[["Pielous_evenness"]]<-add_location(Pielous_evenness, points_item, "Pielous_evenness")
      indices_df[[YYYY]]<-indices
      
      unbias.simp_10km <- rarefy(d2_10km, raremax_10km) - 1
      indices_10km[["unbias.simp"]]<-add_location(unbias.simp_10km, points_item_10km, "unbias.simp")
      Pielous_evenness_10km <- shannon_10km/log(species.richness_10km)
      indices_10km[["Pielous_evenness"]]<-add_location(Pielous_evenness_10km, points_item_10km, "Pielous_evenness")
      indices_10km_df[[YYYY]]<-indices_10km
      
      unbias.simp_100km <- rarefy(d2_100km, raremax_100km) - 1
      indices_100km[["unbias.simp"]]<-add_location(unbias.simp_100km, points_item_100km, "unbias.simp")
      Pielous_evenness_100km <- shannon_100km/log(species.richness_100km)
      indices_100km[["Pielous_evenness"]]<-add_location(Pielous_evenness_100km, points_item_100km, "Pielous_evenness")
      indices_100km_df[[YYYY]]<-indices_100km
    }
  }
  saveRDS(indices_df, target)
  saveRDS(indices_10km_df, sprintf("%s/indices_10km_df.rda", target_folder))
  saveRDS(indices_100km_df, sprintf("%s/indices_100km_df.rda", target_folder))
  
  saveRDS(sp_dis, sprintf("%s/sp_dis.rda", target_folder))
  saveRDS(sp_dis_10km, sprintf("%s/sp_dis_10km.rda", target_folder))
  saveRDS(sp_dis_100km, sprintf("%s/sp_dis_100km.rda", target_folder))
  
  
  
}
if (F){
  df1<-readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Diversity_exposure_0_dispersal_1_xxx/Mammals/UKESM1_SSP585/diversity_df.rda")
  df2<-readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Diversity_exposure_0_dispersal_1/Mammals/UKESM1_SSP585/diversity_df.rda")
  
}
