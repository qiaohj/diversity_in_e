library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
VARs<-c("pr", "tasmax", "tasmin")

start_range<-c(1850:2100)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, Y=start_range)

args = commandArgs(trailingOnly=TRUE)
mm<-as.numeric(args[1])
if (F){
  r<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/MRI-ESM2-0/SSP119/1889/bio13_eck4_1km.tif")
  plot(r)
  r<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/MRI-ESM2-0/SSP119/1889/bio13_eck4_10km.tif")
  plot(r)
  
  mask<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/mask_index_1km.tif")
  no_na<-!is.na(values(mask))
  values(mask)[no_na]<-c(1:length(no_na[no_na==TRUE]))
  plot(mask)
  writeRaster(mask, "/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/mask_index_1km.tif", datatype="INT4U", overwrite=T)
}

if (F){
  mask<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/mask_index.tif")
  for (i in c(1:nrow(start_layer_df))){
    print(paste("Init layer list:", i, nrow(start_layer_df)))
    item<-start_layer_df[i,]
    m=1
    var_tamplate<-"/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Monthly/%s/%s/%s/%d/%s_%d.tif"
    var_tamplate_2<-"/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Monthly/%s/%s/%s/%d/%s_%d_eck4.tif"
    
    for (m in c(mm:mm)){
      cmd_1<-raster(sprintf(var_tamplate, item$GCM, item$SSP, "pr", item$Y, "sum", m))
      cmd_2<-sprintf(var_tamplate_2, item$GCM, item$SSP, "pr", item$Y, "sum", m)
      r<-projectRaster(cmd_1, crs=proj4string(mask))
      writeRaster(r, cmd_2, overwrite=T)
      
      cmd_1<-raster(sprintf(var_tamplate, item$GCM, item$SSP, "tasmax", item$Y, "max", m))
      cmd_2<-sprintf(var_tamplate_2, item$GCM, item$SSP, "tasmax", item$Y, "max", m)
      r<-projectRaster(cmd_1, crs=proj4string(mask))
      writeRaster(r, cmd_2, overwrite=T)
      
      cmd_1<-raster(sprintf(var_tamplate, item$GCM, item$SSP, "tasmin", item$Y, "min", m))
      cmd_2<-sprintf(var_tamplate_2, item$GCM, item$SSP, "tasmin", item$Y, "min", m)
      r<-projectRaster(cmd_1, crs=proj4string(mask))
      writeRaster(r, cmd_2, overwrite=T)
    }
  }
  r<-raster("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/EC-Earth3-Veg/SSP119/1850/bio1_eck4_1km.tif")
  plot(r)
  range(values(r))
  v<-values(r)
  head(v[!is.na(v)])
}

if (T){
  i=1
  cmds<-c()
  cmd_template_1<-"gdalwarp -co COMPRESS=LZW -r bilinear -dstnodata -999999 -ot Int32 -tr 1000 1000 %s %s"
  cmd_template_2<-"gdalwarp -co COMPRESS=LZW -r bilinear -dstnodata -999999 -ot Int32 -tr 10000 10000 %s %s"
  
  GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  #VARs<-c("pr", "tasmax", "tasmin")
  
  start_range<-c(1850:2100)
  start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs, Y=start_range)
  i=1
  index<-1
  xxx<-1
  for (i in c(1:nrow(start_layer_df))){
    print(paste("Init layer list:", i, nrow(start_layer_df)))
    item<-start_layer_df[i,]
    target_folder<-sprintf("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d", 
                           item$GCM, item$SSP, item$Y)
    j=1
    
    #for (j in c(1:19)){
    for (j in c(1,5,6,12,13,14)){
      source<-sprintf("%s/bio%d_eck4.tif", target_folder, j)
      target<-sprintf("%s/bio%d_eck4_1km.tif", target_folder, j)
      if (!file.exists(target)){
        cmd<-sprintf(cmd_template_1, source, target)
        cmds<-c(cmds, cmd)
        index<-index+1
      }
      
      source<-sprintf("%s/bio%d_eck4.tif", target_folder, j)
      target<-sprintf("%s/bio%d_eck4_10km.tif", target_folder, j)
      if (!file.exists(target)){
        cmd<-sprintf(cmd_template_2, source, target)
        cmds<-c(cmds, cmd)
        index<-index+1
      }
    }
    if (index==100){
      
      write.table(cmds, file=sprintf("../../temp/gdalwarp_%d.sh", xxx), row.names = F, col.names=F, quote=F)
      index<-1
      xxx<-xxx+1
      cmds<-c()
    }
    
  }
  #write.table(cmds, file="../../temp/gdalwarp.txt", row.names = F, quote=F)
}
