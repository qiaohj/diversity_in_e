library(raster)
library(ggplot2)
library(dplyr)
library(dplyr)
library(Rmisc)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("functions.r")

group<-"Amphibians"


GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2015:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

i=1
j=1
k=2
#dispersals<-data.frame(M=c(1:5, rep(1, 4), 2, 0, -1), N=c(rep(1,5), c(2:5), 2, 1, 1))
dispersals<-data.frame(M=c(0:5), N=1)

mask<-raster("../../Raster/mask_index.tif")
points<-data.frame(rasterToPoints(mask))
alt<-raster("../../Raster/ALT/alt_eck4.tif")
slope<-raster("../../Raster/ALT/slope_eck4.tif")

for (j in c(1:nrow(layer_df))){
  layer<-layer_df[j,]
  for (k in c(1:nrow(dispersals))){
    layer$M<-dispersals[k, "M"]
    layer$N<-dispersals[k, "N"]
    
    target_folder<-sprintf("../../Objects/Diversity/%s/%s_%d_%d", group, layer$LABEL, layer$M, layer$N)
    target<-sprintf("%s/loss_gain.rda", target_folder)
    df<-readRDS(target)
    df_end<-df%>%dplyr::filter(YEAR==2100)
    df_end$alt<-extract(alt, df_end[, c("x", "y")])
    df_end$slope<-extract(slope, df_end[, c("x", "y")])
    plot(df_end$n_loss, df_end$n_gain)
    plot(df_end$alt, df_end$n_loss)
    plot(df_end$alt, df_end$n_gain)
    plot(df_end$slope, df_end$n_loss)
    plot(df_end$slope, df_end$n_gain)
    plot(df_end$y, df_end$n_loss)
    plot(df_end$y, df_end$n_gain)
    df_end$n<-df_end$n_overlap+df_end$n_loss
    df_end$loss_proportion<-df_end$n_loss/df_end$n
    df_end$gain_proportion<-df_end$n_gain/(df_end$n
    plot(df_end$alt, df_end$loss_proportion)
    plot(df_end$alt, df_end$gain_proportion)
    plot(df_end$slope, df_end$loss_proportion)
    plot(df_end$slope, df_end$gain_proportion)
    ggplot(df_end)+geom_tile(aes(x=x, y=y, color=loss_proportion))
    ggplot(df_end)+geom_tile(aes(x=x, y=y, color=n_gain))
    ggplot(df_end)+geom_tile(aes(x=x, y=y, color=n_loss))
  }
}