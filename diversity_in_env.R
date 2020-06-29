library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(heplots)
library(ntbox)

setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
i=1
groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")
group_base<-"../Object/IUCN_Distribution"
raster_base<-"../Raster/IUCN_Distribution"

rasters<-stack(c("../Raster/Bioclim/PCs/Present/pc1.tif", "../Raster/Bioclim/PCs/Present/pc2.tif"))
p<-data.frame(rasterToPoints(rasters))


species<-list.files(sprintf("%s/%s", group_base, groups[i]), pattern = "\\.rda")
species<-gsub("\\.rda", "", species)
species<-species[!grepl("\\.", species)]

sp<-species[3]
mve_list<-list()
NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
source("addEllipse.R")
source("genCircle.R")
in_Ellipsoid <- stats::qchisq(0.95, 2)

plot(x=p$pc1, y=p$pc2, pch=".", xlim=c(-10, 15), ylim=c(-18, 6))
all_p<-NULL
for (j in c(1:length(species))){
  sp<-species[j]
  print(paste(j, length(species), sp))
  mve_file<-sprintf("%s/%s/%s.mve.rda", group_base, groups[i], sp)
  fit_file<-sprintf("%s/%s/%s.fit.rda", group_base, groups[i], sp)
  r_file<-sprintf("%s/%s/%s.tif", raster_base, groups[i], sp)
  p_file<-sprintf("%s/%s/%s.rda", group_base, groups[i], sp)
  best_ellipse_file<-sprintf("%s/%s/%s.best_ellipse.rda", group_base, groups[i], sp)
  p_tt<-readRDS(p_file)
  if (is.na(p_tt)){
    next()
  }
  if (nrow(p_tt)==0){
    next()
  }
  mve<-NULL
  if (file.exists(mve_file)){
    mve<-readRDS(mve_file)
    fit<-readRDS(fit_file)
  }else{
    mve<-raster(r_file)
  }
  #mve_list[[sp]]<-mve
  p_item<-p
  if (class(mve)=="RasterLayer"){
    p_item$v<-extract(mve, p_item[, c("x", "y")])
    p_item<-p_item %>% dplyr::filter(!is.na(v))
    if (nrow(p_item)==0){
      next()
    }
    p_item$v<-1
  }else{
    p_item$v <- stats::mahalanobis(p_item[, c("pc1", "pc2")], center = mve$centroid, 
                                  cov = mve$covariance)
    p_item<-p_item%>%filter(p_item$v<=in_Ellipsoid)
    p_item$v<-1
    addEllipse(mve$centroid, mve$covariance, col="red", p.interval=0.95)
    
    if (F){
      plot(p_item$pc1, p_item$pc2, pch=".")
      addEllipse(mve$centroid, mve$covariance, col="red", p.interval=0.95)
      points(p_tt$pc1, p_tt$pc2, col="blue")
      points(mve$centroid[1], mve$centroid[2], col="red", pch=3)
    }
  }
  p_item$sp<-sp
  if (is.null(all_p)){
    all_p<-p_item
  }else{
    all_p<-bind_rows(all_p, p_item)
  }
}

