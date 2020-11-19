library(raster)
#library(rgdal)
library(dplyr)
library(vegan)
library(data.table)
library(ggplot2)
library(tidyr)

rm(list=ls())
setDTthreads(threads=1)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (is.na(group)){
  group<-"Amphibians"
}

threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-5
}
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
j=1
k=2
dispersals<-c(0:2)

mask<-raster("../../Raster/mask_index.tif")
points<-as_tibble(rasterToPoints(mask))
add_location<-function(indices, location, type){
  location$metric<-indices
  location$type<-type
  location
}


for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:length(dispersals))){
    layer$M<-dispersals[k]
    
    if (threshold==5){
      target_folder<-sprintf("../../Objects/Diversity_%d/%s/%s_%d", threshold, group, layer$LABEL, layer$M)
    }else{
      target_folder<-sprintf("../../Objects/Diversity/%s/%s_%d", group, layer$LABEL, layer$M)
    }
    
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
    diversity_df<-readRDS(sprintf("%s/diversity_df.rda", target_folder))
    YYYY<-names(diversity_df)[1]
    indices_df<-list()
    sp_dis<-NULL
    for (YYYY in names(diversity_df)){
      indices<-list()
      diversity<-diversity_df[[YYYY]]
      print(paste("BINDING DATA", YYYY, target_folder))
      diversity<-rbindlist(diversity)
      diversity<-unique(diversity)
      print(paste("COUNTING DATA", YYYY, target_folder))
      n_dis<-diversity%>%dplyr::group_by(year, sp)%>%
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
      if (is.null(sp_dis)){
        sp_dis<-n_dis
      }else{
        sp_dis<-bind_rows(sp_dis, n_dis)
      }
      print(paste("EXTRACTING DATA", YYYY, target_folder))
      diversity$index<-raster::extract(mask, diversity[, c("x", "y")])
      diversity<-diversity[complete.cases(diversity),]
      diversity$v<-1
      print(paste("Matrix DATA", YYYY, target_folder))
      t_m<-pivot_wider(diversity, id_cols=index, names_from = sp, values_from=v, values_fill=0)
      d2<-t_m[,-1]
      points_item<-left_join(t_m[,1], points, by=c("index"="mask_index"))
      
      print(paste("simpson", YYYY, target_folder))
      simpson <- diversity(d2, "simpson")
      indices[["simpson"]]<-add_location(simpson, points_item, "simpson")
      print(paste("shannon", YYYY, target_folder))
      shannon<-diversity(d2, "shannon")
      indices[["shannon"]]<-add_location(shannon, points_item, "shannon")
      print(paste("inv simpson", YYYY, target_folder))
      invsimp <- diversity(d2, "inv")
      indices[["invsimp"]]<-add_location(invsimp, points_item, "invsimp")
      #print(paste("fisher.alpha", YYYY, target_folder))
      #alpha <- fisher.alpha(d2)
      #indices[["alpha"]]<-add_location(alpha, points_item, "alpha")
      print(paste("species.richness", YYYY, target_folder))
      species.richness <- specnumber(d2) # observed number of species
      indices[["species.richness"]]<-add_location(species.richness, points_item, "species.richness")
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
      
      ## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
      print(paste("unbias.simp", YYYY, target_folder))
      unbias.simp <- rarefy(d2, raremax) - 1
      indices[["unbias.simp"]]<-add_location(unbias.simp, points_item, "unbias.simp")
      print(paste("Pielous_evenness", YYYY, target_folder))
      Pielous_evenness <- shannon/log(species.richness)
      indices[["Pielous_evenness"]]<-add_location(Pielous_evenness, points_item, "Pielous_evenness")
      indices_df[[YYYY]]<-indices
    }
    saveRDS(indices_df, target)
    saveRDS(sp_dis, sprintf("%s/sp_dis.rda", target_folder))
    
  }
}
