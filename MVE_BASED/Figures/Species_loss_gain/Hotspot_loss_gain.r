setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sp)
library(Rmisc)
library(ggpubr)
f<-"Diversity"
g<-"Amphibians"

for (f in c("Diversity", "Diversity_with_human")){
  for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
    print(paste(f, g))
    df<-readRDS(sprintf("../../Figures/Species_gain_loss/%s_%s.rda", f, g))
    if (F){
      df_sample<-df[sample(nrow(df), 1000),]
      plot(df_sample$slope, df_sample$mean_gain_loss)
    }
    
    source("functions.r")
    
  
    df[which((df$M==0)&(df$mean_gain_loss>0)), "mean_gain_loss"]<-0
    df$gain_loss<-"GAIN"
    df[which(df$mean_gain_loss<0), ]$gain_loss<-"LOSS"
    df[which(df$mean_gain_loss==0), ]$gain_loss<-"UNCHANGED"
  }
}