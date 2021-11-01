GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
years<-c(1901:2100)

cmds<-expand.grid(GCM=GCMs, SSP=SSPs, y=years)

i=1

ff<-c()
for (i in c(1:nrow(cmds))){
  item<-cmds[i,]
  cmd<-sprintf("mv /media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s/%d /media/huijieqiao/QNAS/Sp_Richness_GCM/Raster/ENV/Bioclim/%s/%s", 
               item$GCM, item$SSP, item$y, item$GCM, item$SSP)
  ff<-c(ff, cmd)
}

write.table(ff, "../../temp/mv2.sh", row.names = F, col.names = F, quote=F)
