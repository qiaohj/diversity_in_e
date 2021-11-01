library(raster)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
rm(list=ls())
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
years<-c(1850:2100)
bioclims<-c(1, 5, 6, 12, 13, 14)
labels<-c("1km", "10km")
cmds<-expand.grid(GCM=GCMs, SSP=SSPs, y=years, bioclim=bioclims, label=labels)

i=2
p<-data.frame(x=17240000, y=-8410000)
#cmds<-cmds[sample(nrow(cmds), nrow(cmds)),]
rm<-c()
for (i in c(1:nrow(cmds))){
  item<-cmds[i,]
  
  tif<-sprintf("/media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d/bio%d_eck4_%s.tif", 
               item$GCM, item$SSP, item$y, item$bioclim, item$label)
  print(paste(i, nrow(cmds), tif))
  if (file.exists(tif)){
    
    v<-NA
    tryCatch(
      {
        r<-raster(tif)  
        v<<-extract(r, p)
      },
      error=function(cond) {
        v<<-NA
        
      }
    )
    
    if (is.na(v)){
      rm<-c(rm, sprintf("rm -rf %s", tif))
    }
  }
  
}

write.table(rm, "../../temp/rm_error_tif.sh", row.names = F, col.names = F, quote=F)
