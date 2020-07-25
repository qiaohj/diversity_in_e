library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)
library(heplots)
library(ntbox)
library(ggplot2)

setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
source("functions.R")
source("addEllipse.R")
source("genCircle.R")
in_Ellipsoid <- stats::qchisq(0.95, 2)
var_names<-c("pc1", "pc2")
i=1
groups<-c("Amphibians", "Birds", "Mammals", "Reptiles")
group_base<-"../Object/bioclim 2.0/IUCN_Distribution/PCA"
raster_base<-"../Raster/bioclim 2.0/IUCN_Distribution"
mask<-raster("../Raster/mask.tif")

species_file<-list.files(sprintf("%s/%s", group_base, groups[i]), pattern = "\\.rda", full.names=T)
species<-list.files(sprintf("%s/%s", group_base, groups[i]), pattern = "\\.rda")
species<-gsub("\\.rda", "", species)
species<-species[!grepl("\\.", species)]
j=4
sp<-species[j]
ppp<-readRDS(species_file[j])
fit <- cov.rob(ppp[, var_names], quantile.used=NDquntil(nrow(ppp), 0.95),  method = "mve")

rasters<-stack(c(
  sprintf("/Volumes/Disk2/Experiments/Diversity_in_Env/Raster/Bioclim 2.0/PCs/Future/bc%dbi70/pc1.tif", 26),
  sprintf("/Volumes/Disk2/Experiments/Diversity_in_Env/Raster/Bioclim 2.0/PCs/Future/bc%dbi70/pc2.tif", 26)
))
p_all<-data.frame(rasterToPoints(rasters))

p<-p_all
p$v <- stats::mahalanobis(p[, c("pc1", "pc2")], center = fit$center, 
                               cov = fit$cov)
p<-p%>%filter(p$v<=in_Ellipsoid)

plot(p$pc1, p$pc2, pch=".")
points(ppp$pc1, ppp$pc2, col="blue")
addEllipse(fit$center, fit$cov, col="red", p.interval=0.95)
point<-p[1, c("x", "y")]
points<-ppp
min_dist<-function(x, y, points){
  min(sqrt((x-points$x)^2+(y-points$y)^2))
}
max_dist<-5
p<-p%>%dplyr::rowwise()%>%dplyr::mutate(dist=min_dist(x, y, ppp)/100000)
p$probability<-1-p$dist/max_dist
p$random<-0
p<-p%>%dplyr::filter(probability>=random)

r_ppp <- rasterFromXYZ(ppp[, c("x", "y", "layer")])
r_p <- rasterFromXYZ(p[, c("x", "y", "dist")])

plot(r_p)
plot(r_ppp, col="red", add=T, legend=F)


for (k in c(kk)){
  
  if (k==0){
    rasters<-stack(c("../Raster/Bioclim/PCs/Present/pc1.tif", "../Raster/Bioclim/PCs/Present/pc2.tif"))
    p<-data.frame(rasterToPoints(rasters))
    png(filename=sprintf("../Figures/Niche_E/%s.png", groups[i]), width=1000, height=700)
  }else{
    com<-coms[k,]
    rasters<-stack(c("../Raster/Bioclim/PCs/Present/pc1.tif", "../Raster/Bioclim/PCs/Present/pc2.tif"))
    p_present<-data.frame(rasterToPoints(rasters))
    
    rasters<-stack(c(
      sprintf("../Raster/Bioclim/PCs/Future/%s%sbi70/pc1.tif", com$gcm, com$rcp),
      sprintf("../Raster/Bioclim/PCs/Future/%s%sbi70/pc2.tif", com$gcm, com$rcp)
    ))
    p<-data.frame(rasterToPoints(rasters))
    png(filename=sprintf("../Figures/Niche_E/%s_%s%s.png", groups[i], com$gcm, com$rcp), width=1000, height=700)
  }
  
  
  
  if (k==0){
    plot(x=p$pc1, y=p$pc2, pch=".", xlim=c(-10, 15), ylim=c(-18, 6))
  }else{
    plot(x=p_present$pc1, y=p_present$pc2, pch=".", xlim=c(-10, 15), ylim=c(-18, 6), col=alpha("black", alpha=0.2))
    points(x=p$pc1, y=p$pc2, pch=".", xlim=c(-10, 15), ylim=c(-18, 6), col="#005AB5")
  }
  
  all_p<-NULL
  for (j in c(1:length(species))){
    sp<-species[j]
    print(paste(k, nrow(coms), j, length(species), sp))
    mve_file<-sprintf("%s/%s/%s.mve.rda", group_base, groups[i], sp)
    fit_file<-sprintf("%s/%s/%s.fit.rda", group_base, groups[i], sp)
    r_file<-sprintf("%s/%s/%s.tif", raster_base, groups[i], sp)
    p_file<-sprintf("%s/%s/%s.rda", group_base, groups[i], sp)
    best_ellipse_file<-sprintf("%s/%s/%s.best_ellipse.rda", group_base, groups[i], sp)
    p_tt<-readRDS(p_file)
    if (class(p_tt)=="logical"){
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
      if (nrow(p_item)==0){
        next()
      }
      p_item$v<-1
      addEllipse(mve$centroid, mve$covariance, col="red", p.interval=0.95)
      
      if (F){
        plot(p_tt$pc1, p_tt$pc2, pch=".")
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
  
  dev.off()
  all_p_sum<-all_p%>%dplyr::group_by(x, y, pc1, pc2)%>%dplyr::summarise(n_sp=n())
  
  if (k==0){
    saveRDS(all_p, sprintf("../Object/Niche_E/%s.rda", groups[i]))
    g<-ggplot(all_p_sum, aes(x=pc1, y=pc2, color=n_sp))+geom_point()+theme_bw()+
      scale_colour_gradient(
        low = "#005AB5",
        high = "#DC3220",
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    ggsave(g, file=sprintf("../Figures/Niche_E/%s_merged.png", groups[i]))
    
    p_full<-SpatialPointsDataFrame(all_p_sum[, c("x", "y")], 
                                   all_p_sum, proj4string = crs(mask))
    saveRDS(p_full, sprintf("../Object/Niche_E/%s_SpatialPoints.rda", groups[i]))
    r <- rasterFromXYZ(as.data.frame(all_p_sum)[, c("x", "y", "n_sp")], res=res(mask), crs=crs(mask))
    plot(r)
    writeRaster(r, sprintf("../Raster/Niche_E/%s_present.tif", groups[i]), overwrite=TRUE)
    
  }else{
    saveRDS(all_p, sprintf("../Object/Niche_E/%s_%s%s.rda", groups[i], com$gcm, com$rcp))
    g<-ggplot(p_present, aes(x=pc1, y=pc2), color="black", alpha=0.2)+geom_point()+
      geom_point(data=all_p_sum, aes(x=pc1, y=pc2, color=n_sp))+theme_bw()+
      scale_colour_gradient(
        low = "#005AB5",
        high = "#DC3220",
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    ggsave(g, file=sprintf("../Figures/Niche_E/%s_%s%s_merged.png", groups[i], com$gcm, com$rcp))
    
    p_full<-SpatialPointsDataFrame(all_p_sum[, c("x", "y")], 
                                   all_p_sum, proj4string = crs(mask))
    saveRDS(p_full, sprintf("../Object/Niche_E/%s_%s%s_SpatialPoints.rda", groups[i], com$gcm, com$rcp))
    r <- rasterFromXYZ(as.data.frame(all_p_sum)[, c("x", "y", "n_sp")], res=res(mask), crs=crs(mask))
    #plot(r)
    writeRaster(r, sprintf("../Raster/Niche_E/%s_%s%s.tif", groups[i], com$gcm, com$rcp), overwrite=TRUE)
  }
}



