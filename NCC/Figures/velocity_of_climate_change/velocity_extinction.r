library(data.table)
library(raster)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
exposure=0
g<-"Mammals"
for (exposure in c(0, 5)){
  for (g in c("Birds", "Mammals")){
    print(paste(g, exposure))
    df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s_10km.rda", exposure, g))
    df<-rbindlist(df)
    df_se<-df[, .(N=.N), by=list(x, y, mask, GCM, SSP, dispersal)]
    if (F){
      ggplot(df_se[GCM=="UKESM1"&SSP=="SSP585"&dispersal==0])+geom_tile(aes(x=x, y=y, color=N))
    }
  }
}
ratio_final<-readRDS("../../Figures/when_where_extinction_all/ratio_final_group_10km.rda")
mask<-raster("../../Raster/mask_100km.tif")
Raster/VoCC/EC-Earth3-Veg/SSP245