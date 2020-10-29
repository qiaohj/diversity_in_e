setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(raster)
library(dplyr)
library(ggplot2)
f<-"Diversity"
g<-"Amphibians"
df<-readRDS(sprintf("../../Figures/Species_gain_loss/%s_%s.rda", f, g))
df_sample<-df[sample(nrow(df), 1000),]
plot(df_sample$slope, df_sample$mean_gain_loss)
