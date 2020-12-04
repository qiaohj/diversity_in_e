library(raster)
#library(rgdal)
#library(rgeos)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(leaps)

rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")


sp_dis_all_sub_N_all<-NULL
sp_dis_extinct<-NULL
for (threshold in c(1, 5)){
  rda<-sprintf("../../Figures/N_Extinction/sp_dis_all_%d.rda", threshold)
  print(paste("Reading", rda))
  sp_dis_all<-readRDS(rda)
  sp_dis_all_sub<-sp_dis_all%>%dplyr::filter(year==2100)
  sp_dis_all_sub$range_ratio<-sp_dis_all_sub$N_CELL/sp_dis_all_sub$st_N_CELL
  plot(sp_dis_all_sub$range_ratio, sp_dis_all_sub$st_N_CELL)
  sp_dis_all_sub$abs_st_ymean<-abs(sp_dis_all_sub$st_ymean)
  models<-regsubsets(range_ratio ~ st_N_CELL + st_xrange + st_yrange + abs_st_ymean, 
                     data = sp_dis_all_sub, nbest=4)
  plot(models, scale = 'adjr2')
  
}
