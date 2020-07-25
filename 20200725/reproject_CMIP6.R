library(raster)
library(RStoolbox)
library(stringr)
library(factoextra)

setwd("/Volumes/Disk2/Experiments/Diversity_in_Env/Script")


x=1
args<-commandArgs(trailingOnly = T)
x<-as.numeric(args[1])

GCMs<-c("BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1", 
        "MIROC-ES2L", "MIROC6")

GCMs<-GCMs[x]

ssps<-c("ssp126", "ssp245", "ssp370", "ssp585")
years<-c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
vartypes<-c("bioc", "prec", "tmax", "tmin")
coms<-expand.grid(gcm=GCMs, ssp=ssps, year=years, vartype=vartypes, stringsAsFactors = F)
mask<-raster("../Raster/mask.tif")

if (F){
  #future
  templete<-"/Volumes/ST6T/Worldclim2.1/10m/tif/future/cmip6/7_fut/10m/%s/%s/wc2.1_10m_%s_%s_%s_%s.tif"
  for (j in c(1:nrow(coms))){
    
    com<-coms[j,]
    r<-stack(sprintf(templete, com$gcm, com$ssp, com$vartype, com$gcm, com$ssp, com$year))
    if (com$vartype=="bioc"){
      n_var<-19
    }else{
      n_var<-12
    }
    for (i in c(1:n_var)){
      print(paste(j, nrow(coms), i))
      
      s<-r[[i]]
      s<-projectRaster(s, res=res(mask), crs=crs(mask))
      dir.create(sprintf("../Raster/Bioclim 2.1/Future/%s/%s", com$gcm, com$ssp), showWarnings = F, recursive=T)
      writeRaster(s, sprintf(sprintf("../Raster/Bioclim 2.1/Future/%s/%s/%s_%s_%d.tif", 
                                     com$gcm, com$ssp, com$vartype, com$year, i)), 
                  overwrite=TRUE, NAflag=-9999)
    }
  }
  
  #present
  templete<-"/Volumes/ST6T/Worldclim2.1/10m/tif/present/variable/wc2.1_10m_%s_%s.tif"
  vartypes<-c("prec", "srad", "tavg", "tmax", "tmin", "vapr", "wind")
  for (vartype in vartypes){
    for (i in c(1:12)){
      print(paste(vartype, i))
      r<-raster(sprintf(templete, vartype, str_pad(i, 2, pad="0")))
      s<-r
      s<-projectRaster(s, res=res(mask), crs=crs(mask))
      writeRaster(s, sprintf(sprintf("../Raster/Bioclim 2.1/Present/variable/%s_%d.tif", 
                                     vartype, i)), 
                  overwrite=TRUE, NAflag=-9999)
    }
  }
  
  templete<-"/Volumes/ST6T/Worldclim2.1/10m/tif/present/wc2.1_10m_bio/wc2.1_10m_%s_%d.tif"
  vartypes<-c("bio")
  for (vartype in vartypes){
    for (i in c(1:19)){
      print(paste(vartype, i))
      r<-raster(sprintf(templete, vartype, i))
      s<-r
      s<-projectRaster(s, res=res(mask), crs=crs(mask))
      writeRaster(s, sprintf(sprintf("../Raster/Bioclim 2.1/Present/bioclim/%s_%d.tif", 
                                     vartype, i)), 
                  overwrite=TRUE, NAflag=-9999)
    }
  }
  
  r<-raster("/Volumes/ST6T/Worldclim2.1/10m/tif/wc2.1_10m_elev.tif")
  s<-r
  s<-projectRaster(s, res=res(mask), crs=crs(mask))
  writeRaster(s, "../Raster/Bioclim 2.1/elevation.tif", 
              overwrite=TRUE, NAflag=-9999)
  
  
  #PCA
  rasters<-c()
  template<-"../Raster/Bioclim 2.1/Present/bioclim/bio_%d.tif"
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
    writeRaster(pca$map[[i]], sprintf("../Raster/Bioclim 2.1/PCs/Present/pc_%d.tif", i), 
                overwrite=T, NAflag=-9999)
  }
  
  
}

pca<-readRDS("../Object/pca_bioclim.rda")
for (j in c(1:nrow(coms))){
  com<-coms[j,]
  if (com$vartype!="bioc"){
    next()
  }
  
  print(paste(j, nrow(coms), com$gcm, com$ssp, com$year, com$vartype))
  
  rasters<-c()
  for (i in c(1:19)){
    rasters<-c(rasters, sprintf("../Raster/Bioclim 2.1/Future/%s/%s/%s_%s_%d.tif", 
                                com$gcm, com$ssp, com$vartype, com$year, i))
  }
  rasters <- stack(rasters)
  p<-rasterToPoints(rasters)
  p<-data.frame(p)
  colnames(p)[3:21]<-paste("BIO", c(1:19))
  
  p<-p[complete.cases(p),]
  dim(p)
  predicted<-predict(pca$model, p[, c(3:21)])
  no_na<-NULL
  dir.create(sprintf("../Raster/Bioclim 2.1/PCs/Future/%s/%s", com$gcm, com$ssp), 
             showWarnings = F, recursive=T)
  for (i in c(1:19)){
    print(paste(j, nrow(coms), i))
    r<-raster(sprintf("../Raster/Bioclim 2.1/Future/%s/%s/%s_%s_%d.tif", 
                      com$gcm, com$ssp, com$vartype, com$year, i))
    if (is.null(no_na)){
      no_na<-!is.na(values(r))
    }
    
    values(r)<-NA
    values(r)[which(no_na)]<-predicted[, i]
    
    writeRaster(r, sprintf("../Raster/Bioclim 2.1/PCs/Future/%s/%s/pc_%d.tif", 
                           com$gcm, com$ssp, i), overwrite=TRUE, NAflag=-9999)
  }
  
}
