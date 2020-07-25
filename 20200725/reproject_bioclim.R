library(raster)
library(RStoolbox)
library(stringr)
library(factoextra)

setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")
templete<-"/Volumes/Disk2/Experiments/Protected_Diversity/GIS/eck4_10km/bio_10m_tiff/bio%d.tif"
i=1

for (i in c(1:19)){
  print(i)
  r<-raster(sprintf(templete, i))
  s<-r
  res(s)<-c(100000, 100000)
  s <- resample(r, s, method='bilinear')
  writeRaster(s, sprintf("../Raster/Bioclim/Present/bio%d.tif", i), overwrite=TRUE)
}

templete<-"/Volumes/Disk2/Experiments/Protected_Diversity/GIS/eck4_10km/Future/10m/%s%sbi70/%s%sbi70%d.tif"
GCMs<-c("bc", "cc", "gs", "hd", "he", "ip", "mc", "mg", "mi", "mr", "no")
rcps<-c("26", "45", "60", "85")
coms<-expand.grid(gcm=GCMs, rcp=rcps, stringsAsFactors = F)
for (j in c(1:nrow(coms))){
  
  com<-coms[j,]
  for (i in c(1:19)){
    print(paste(j, nrow(coms), i))
    r<-raster(sprintf(templete, com$gcm, com$rcp, com$gcm, com$rcp, i))
    s<-r
    res(s)<-c(100000, 100000)
    s <- resample(r, s, method='bilinear')
    dir.create(sprintf("../Raster/Bioclim/Future/%s%sbi70", com$gcm, com$rcp), showWarnings = F)
    writeRaster(s, sprintf(sprintf("../Raster/Bioclim/Future/%s%sbi70/%s%sbi70%d.tif", com$gcm, com$rcp, com$gcm, com$rcp, i)), 
                overwrite=TRUE)
  }
}

#PCA
rasters<-c()
template<-"../Raster/Bioclim/Present/bio%d.tif"
for (i in c(1:19)){
  rasters<-c(rasters, sprintf(template, i))
}
rasters <- stack(rasters)

pca <- rasterPCA(rasters, spca=T) 

summary(pca$model)
t(pca$model$sdev^2 / sum(pca$model$sdev^2))
fviz_eig(pca$model)
rownames(pca$model$loadings)<-paste("BIO", c(1:19))
fviz_pca_var(pca$model,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
saveRDS(pca, "../Object/pca_bioclim.rda")
for (i in c(1:19)){
  print(i)
  writeRaster(pca$map[[i]], sprintf("../Raster/Bioclim/PCs/Present/pc%d.tif", i), overwrite=T)
}

for (j in c(34:nrow(coms))){
  print(paste(j, nrow(coms)))
  com<-coms[j,]
  rasters<-c()
  for (i in c(1:19)){
    rasters<-c(rasters, sprintf("../Raster/Bioclim/Future/%s%sbi70/%s%sbi70%d.tif", com$gcm, com$rcp, com$gcm, com$rcp, i))
  }
  rasters <- stack(rasters)
  p<-rasterToPoints(rasters)
  p<-data.frame(p)
  colnames(p)[3:21]<-paste("BIO", c(1:19))
  
  p<-p[complete.cases(p),]
  dim(p)
  predicted<-predict(pca$model, p[, c(3:21)])
  no_na<-NULL
  for (i in c(1:19)){
    print(paste(j, nrow(coms), i))
    r<-raster(sprintf("../Raster/Bioclim/Future/%s%sbi70/%s%sbi70%d.tif", com$gcm, com$rcp, com$gcm, com$rcp, i))
    if (is.null(no_na)){
      no_na<-!is.na(values(r))
    }

    values(r)<-NA
    values(r)[which(no_na)]<-predicted[, i]
    dir.create(sprintf("../Raster/Bioclim/PCs/Future/%s%sbi70", com$gcm, com$rcp), showWarnings = F)
    writeRaster(r, sprintf("../Raster/Bioclim/PCs/Future/%s%sbi70/pc%d.tif", com$gcm, com$rcp, i), overwrite=TRUE)
  }
  
}

