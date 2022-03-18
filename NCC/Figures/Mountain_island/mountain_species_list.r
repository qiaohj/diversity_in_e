library(sf)
library(raster)
library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_10km.tif")
if (F){
  mountain<-raster("../../Objects/Island/mountain_index_10km.tif")
  extant_mask<-raster::extent(mountain)
  all_info<-readRDS("../../Objects/Island/species_mountain_10km.rda")
  mountain_species<-all_info[mountain>0]
  i=1
  all_mountain_species<-list()
  for (i in c(1:nrow(mountain_species))){
    print(paste(i, nrow(mountain_species)))
    sp_item<-mountain_species[i]
    sp<-sp_item$sp
    target<-sprintf("../../Objects/IUCN_Distribution/%s/RAW/%s.rda", sp_item$group, sp)
    if (!file.exists(target)){
      next()
    }
    dis<-readRDS(target)
    if (nrow(dis)==0){
      next()
    }
    extend<-st_bbox(dis)
    extend<-c(extend[1], extend[3], extend[2], extend[4])
    
    if (between(extend[1], extant_mask[1], extant_mask[2])&
        between(extend[2], extant_mask[1], extant_mask[2])&
        between(extend[3], extant_mask[3], extant_mask[4])&
        between(extend[4], extant_mask[3], extant_mask[4])){
      dis_c<-crop(mountain, extend)
      dis_c<-mask(dis_c, dis)
      v<-values(dis_c)
      v_index<-data.frame(table(v))
      v_index$sp<-sp
      v_index$group<-sp_item$group
      all_mountain_species[[sp]]<-v_index
    }else{
      print("no overlap, skip")
    }
  }
  all_mountain_species<-rbindlist(all_mountain_species)
  colnames(all_mountain_species)<-c("index", "area", "sp", "group")
  mountain_index<-readRDS("../../Objects/Island/mountain_id_list.rda")
  all_mountain_species$index<-as.numeric(as.character(all_mountain_species$index))
  all_mountain_species<-merge(all_mountain_species, mountain_index, by="index")
  saveRDS(all_mountain_species, "../../Objects/Island/all_mountain_species.rda")
}


if (F){
  all_mountain_species<-readRDS("../../Objects/Island/all_mountain_species.rda")
  extinct_sp_list<-list()
  for (exposure in c(0, 5)){
    for (dispersal in c(0, 1)){
      
      extinct_sp<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/extinct_sp_%d_10km_2_100km.rda", 
                                  exposure, dispersal))
      cols<-c("sp", "group", "GCM", "SSP")
      extinct_sp<-extinct_sp[, ..cols]
      extinct_sp$exposure<-exposure
      extinct_sp$dispersal<-dispersal
      extinct_sp_list[[length(extinct_sp_list)+1]]<-extinct_sp
    }
  }
  extinct_sp_list<-rbindlist(extinct_sp_list)
  all_mountain_species_label<-merge(extinct_sp_list, all_mountain_species,
                                    by=c("sp", "group"))
  all_mountain_species_label<-
    all_mountain_species_label[, .(N=.N), by=list(GCM, SSP, exposure, dispersal, Name, Country)]
  coms<-expand.grid(GCM=unique(all_mountain_species_label$GCM),
                    SSP=unique(all_mountain_species_label$SSP),
                    exposure=unique(all_mountain_species_label$exposure),
                    dispersal=unique(all_mountain_species_label$dispersal),
                    Name=unique(all_mountain_species_label$Name),
                    stringsAsFactors = F)
  all_mountain_species_label<-merge(coms, all_mountain_species_label,
                                    by=c("GCM", "SSP", "exposure", "dispersal", "Name"),
                                    all.x=T)
  all_mountain_species_label<-data.table(all_mountain_species_label)
  all_mountain_species_label[is.na(N)]$N<-0
  
  all_mountain_species_label<-all_mountain_species_label[,.(N=mean(N)),
                                                         by=c("SSP", "exposure", "dispersal", "Name")]
  
  all_mountain_species_number<-all_mountain_species[, .(N_sp=.N),
                                                    by=list(Name, index)]
  
  
  all_mountain_species_info<-merge(all_mountain_species_label, all_mountain_species_number, 
                                   by=c("Name"), all=T)
  
  all_mountain_species_info$extinct_proportion<-
    all_mountain_species_info$N/all_mountain_species_info$N_sp
  saveRDS(all_mountain_species_info, "../../Objects/Island/all_mountain_species_info.rda")
  all_mountain_f<-all_mountain_species_info[extinct_proportion>0]
  
  all_mountain_f <- all_mountain_f[order(all_mountain_f$extinct_proportion, decreasing = TRUE), ]  
  keycols = c("SSP", "exposure", "dispersal")
  setkeyv(all_mountain_f, keycols)
  top5<-all_mountain_f[, head(.SD, 5), 
                            by=list(SSP, exposure, dispersal)]
  top5<-merge(top5, all_mountain_species[, .(Nspp=.N), by=list(Name, Country)], by="Name")
  write.csv(top5, "../../Figures/Mountain_island/extinction_proportion_top5.csv", row.names = F)
  mountain<-raster("../../Objects/Island/mountain_index_10km.tif")
  mountain_p<-data.table(rasterToPoints(mountain))
  SSPi<-"SSP119"
  cccc<-expand.grid(SSP= unique(all_mountain_species_info$SSP),
                    exposure=c(0, 5), dispersal=c(0, 1), stringsAsFactors = F)
  i=1
  all_p<-list()
  for (i in 1:nrow(cccc)){
    com<-cccc[i,]
    item<-all_mountain_species_info[SSP==com$SSP&dispersal==com$dispersal&exposure==com$exposure]
    p<-merge(mountain_p, item, by.x="mountain_index_10km", by.y="index")
    all_p[[i]]<-p
  }
  all_p<-rbindlist(all_p)
  saveRDS(all_p, "../../Figures/Mountain_island/extinct_proportion.rda")
}
library(ggplot2)
source("commonFuns/colors.r")
all_p<-readRDS("../../Figures/Mountain_island/extinct_proportion.rda")

sspi<-"SSP245"
myPalette <- colorRampPalette(c(mask_color, color_two_map[2]))
colors<-myPalette(max(all_p$extinct_proportion*100))
all_p$exposure<-ifelse(all_p$exposure==0, " no climate resilience", "climate resilience")
all_p$dispersal<-ifelse(all_p$dispersal==0, "no dispersal", "with dispersal")
for (sspi in c("SSP119", "SSP245", "SSP585")){
  p<-ggplot(all_p[SSP==sspi])+
    geom_tile(aes(x=x, y=y, fill=extinct_proportion*100))+
    facet_grid(exposure~dispersal)+
    scale_fill_gradientn(colors=colors)+
    coord_fixed()+
    labs(fill="Extinction proportion")+
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
      panel.border = element_blank()
    )
  ggsave(p, filename=sprintf("../../Figures/Mountain_island/extinct_proportion/extinct_proportion_%s.png", 
                             sspi), width=10, height=5)
}
  
