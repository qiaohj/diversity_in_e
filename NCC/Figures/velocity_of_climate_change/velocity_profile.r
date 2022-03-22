library(raster)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
if (F){
  gcm<-"UKESM1"
  ssp<-"SSP585"
  var=1
  mask<-raster("../../Raster/mask_100km.tif")
  p<-data.table(rasterToPoints(mask))
  all_df<-list()
  var<-1
  for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
    for (ssp in c("SSP119", "SSP245", "SSP585")){
      for (var in c(1,5,6,12,13,14)){
        bio<-raster(sprintf("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Raster/VoCC_mask/%s/%s/bio%d_voccMag.tif", gcm, ssp, var))
        p_t<-p
        p_t$v<-extract(bio, data.frame(x=p$x, y=p$y))
        p_t$var<-sprintf("bio%d", var)
        p_t$SSP<-ssp
        p_t$ESM<-gcm
        
        for (group in c("Birds", "Mammals")){
          for (dispersal in c(0, 1)){
            for (exposure in c(0, 5)){
              print(paste(gcm, ssp, var, group, dispersal, exposure))
              p_tt<-p_t
              p_tt$group<-group
              p_tt$dispersal<-dispersal
              p_tt$exposure<-exposure
              extinct<-raster(sprintf("../../Figures/when_where_extinction_exposure_%d/where/TIFF/%s_%s_%s_extinct_ratio_dispersal_%d_10km.tif",
                                      exposure, group, gcm, ssp, dispersal))
              p_tt$extinct_ratio<-extract(extinct, data.frame(x=p_tt$x, y=p_tt$y))
              extinct<-raster(sprintf("../../Figures/when_where_extinction_exposure_%d/where/TIFF/%s_%s_%s_extinct_n_dispersal_%d_10km.tif",
                                      exposure, group, gcm, ssp, dispersal))
              p_tt$extinct_n<-extract(extinct, data.frame(x=p_tt$x, y=p_tt$y))
              all_df[[length(all_df)+1]]<-p_tt
            }
          }
        }
      }
    }
  }
  all_df<-rbindlist(all_df)
  saveRDS(all_df, "../../Figures/VoCC/all_df.rda")
}
all_df<-readRDS("../../Figures/VoCC/all_df.rda")
rm("var")
vv="bio1"
#for (vv in c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")){
cols<-colnames(all_df)[c(1:7)]
all_df_vocc<-unique(all_df[, ..cols])
plist<-list()
for (vv in c("bio1", "bio12")){
  p_df<-all_df_vocc[var==vv]
  p<-ggplot(p_df)+geom_histogram(aes(x=v, fill=SSP), bins=50)+
    scale_fill_manual(values=color_ssp)+
    theme_bw()+
    xlab("Velocity of climage change in annual temperature")+
    facet_grid(var~ESM, scale="free")+theme(axis.title.x=element_blank())
    if (vv=="bio12"){
      #p<-p+xlim(-1e3, 1.5e3)
    }
  plist[[vv]]<-p
}
p<-ggarrange(plotlist=plist, nrow=2, ncol=1)
p
ggsave(p, filename="../../Figures/VoCC/VoCC_profile.png", width=10, height=5)    
g<-"Mammals"
dis<-0
ex<-0
vv<-"bio1"
island<-raster("../../Objects/Island/islands.tif")
gcm<-"UKESM1"
ssp<-"SSP245"
for (gcm in c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")){
  for (ssp in c("SSP119", "SSP245", "SSP585")){
    for (vv in c("bio1", "bio12")){
      for (g in c("Birds", "Mammals")){
        for (dis in c(0, 1)){
          for (ex in c(0, 5)){
            print(paste(gcm, ssp, vv, group, dispersal, exposure))
            item<-all_df[ESM==gcm&SSP==ssp&var==vv&group==g&dispersal==dis&exposure==ex]
            item$is_tropic<-between(item$y, -2e6, 2e6)
            ggplot(item)+geom_tile(aes(x=x, y=y, fill=item$v))+
              scale_fill_gradient(low=color_two_map[1], high=color_two_map[2])+
              theme(
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_rect(fill = map_background, color = NA), 
                panel.background = element_blank(), 
                legend.background = element_rect(fill = map_background, color = NA),
                panel.border = element_blank())
            
            item<-item[!is.na(extinct_ratio)]
            item<-item[extinct_ratio<1]
            
            item$island<-extract(island, data.frame(x=item$x, y=item$y))
            model<-glm(data=item, extinct_ratio~v)
            item$v_int<-round(item$v)
            item_se<-item[, .(mean_ratio=mean(extinct_ratio)), by=list(island, v_int)]
            ggplot(item, aes(x=v, y=extinct_ratio))+geom_point(aes(color=factor(island)))+
              facet_grid(is_tropic~island, scale="free")
            ggplot(item_se, aes(x=abs(v_int), y=mean_ratio))+geom_point(aes(color=factor(island)))+
              geom_smooth(method="glm")+facet_grid(is_tropic~island, scale="free")
          }
        }
      }
    }
  }
}
