library(raster)
library(dplyr)
library(concaveman)
library(sf)
args = commandArgs(trailingOnly=TRUE)
group<-args[1]
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mask<-raster("../../Raster/mask_index.tif")
mask_p<-data.frame(rasterToPoints(mask))
if (is.na(group)){
  group<-"Amphibians"
}
threshold<-as.numeric(args[2])
if (is.na(threshold)){
  threshold<-1
}
if (F){
  r<-raster("../../Raster/Continent.tif")
  plot(r)
  p<-data.frame(rasterToPoints(r))
  p[which(p$x<=-167), "Continent"]<-NA
  values(r)[!is.na(values(r))]<-p$Continent
  
  
  r2<-projectRaster(r, mask, crs(mask), method="ngb")
  p2<-data.frame(rasterToPoints(r2))
  p2[which(p2$x<=-12103059), "Continent"]<-NA
  
  p2[which((p2$x>12912000)&(p2$y>5000000)), "Continent"]<-NA
  p2[which(p2$Continent==4), "Continent"]<-2
  values(r2)[!is.na(values(r2))]<-p2$Continent
  plot(r2)
  writeRaster(r2, "../../Raster/Continent_ect4.tif", overwrite=T)
  plot(r2)
}
#threshold<-5
continent<-raster("../../Raster/Continent_ect4.tif")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")

df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))
i=1
#dispersals<-data.frame(M=c(0:5, rep(1, 4), 2), N=c(rep(1,6), c(2:5), 2))
dispersals<-c(1)
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
final_df<-NULL
colors<-rainbow(length(2021:2100))
for (i in c(1:nrow(df_list))){
  print(paste(i, nrow(df_list)))
  item<-df_list[i,]
  item$sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, item$sp)
  
  target<-sprintf("%s/dispersal_%d", target_folder, threshold)
  
  j=1
  no_na<-!is.na(values(mask))
  for (j in c(1:nrow(layer_df))){
    layer_item<-layer_df[j,]
    k=1
    for (k in c(1:length(dispersals))){
      dispersal<-dispersals[k]
      ttt<-sprintf("../../Objects/dispersal_path_%d/%s/%s_%s_%d.rda", 
                   threshold, group, item$sp, layer_item$LABEL, dispersal)
      if (file.exists(ttt)){
        #ddddd<-readRDS(ttt)
        #if (!is.null(ddddd)){
          next()
        #}
        
      }
      dir.create(sprintf("../../Objects/dispersal_path_%d/%s", 
                         threshold, group), showWarnings = F, recursive = T)
      saveRDS(NULL, ttt)
      dispersal_log<-readRDS(sprintf("%s/%s_%d.rda", target, layer_item$LABEL, dispersal))
      if (is.null(dispersal_log)){
        next()
      }
      
      start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
      #start_dis<-start_dis%>%dplyr::filter(in_out==1)
      start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y, mask_index)
      start_dis$YEAR<-2020
      dispersal_log<-bind_rows(dispersal_log[, c("x", "y", "mask_index", "YEAR")], start_dis)
      dispersal_log$continent_i<-raster::extract(continent, dispersal_log[, c("x", "y")])
      dispersal_log<-dispersal_log%>%dplyr::filter(!is.na(continent_i))
      if (nrow(dispersal_log)==0){
        next()
      }
      year=2020
      if (F){
        for (year in c(2020:2100)){
          item_y<-dispersal_log%>%dplyr::filter(YEAR==year)
          item_y$N<-1
          p<-left_join(mask_p, item_y, by="mask_index")
          #plot(item_y$x, item_y$y, xlim=range(dispersal_log$x), ylim=range(dispersal_log$y))
          if (F){
            r<-mask
            values(r)[no_na]<-p$N
            plot(r, main=year, col=colors[year-2019])
          }
          v<-readline(prompt="X=exit ")
          if (toupper(v)=="X"){
            break()
          }
        }
      }
      dispersal_log_se<-dispersal_log%>%dplyr::group_by(YEAR, continent_i)%>%
        dplyr::summarise(gravity_x=mean(x), gravity_y=mean(y))
      saveRDS(dispersal_log_se, ttt)
      if (F){
        keyyears<-c(2020, 2040, 2060, 2080, 2100)
        ggplot(dispersal_log_se)+geom_path(aes(x=gravity_x, y=gravity_y, color=YEAR))
        ggplot(dispersal_log_se%>%dplyr::filter(YEAR %in% keyyears))+geom_path(aes(x=gravity_x, y=gravity_y, color=YEAR))
        ggplot(start_dis)+geom_point(aes(x=x, y=y))+
          geom_point(data=dispersal_log_se%>%filter(YEAR==2014), aes(x=gravity_x, y=gravity_y, color="red"))
        
        lsdf <- list()
        plot.new()
        sm_path<-data.frame(xspline(dispersal_log_se$gravity_x, dispersal_log_se$gravity_y, shape=1, draw=F))
        dispersal_log_se_few<-dispersal_log_se%>%dplyr::filter(YEAR %in% keyyears)
        sm_path<-data.frame(xspline(dispersal_log_se_few$gravity_x, dispersal_log_se_few$gravity_y, shape=1, draw=F, repEnds=T))
        sm_path$YEAR<-seq(keyyears[1], keyyears[length(keyyears)], by=(keyyears[length(keyyears)]-keyyears[1])/(nrow(sm_path)-1))
        ggplot(dispersal_log_se_few)+
          #geom_path(aes(x=gravity_x, y=gravity_y, color=YEAR))+
          geom_path(data=sm_path, aes(x=x, y=y, color=YEAR))
        
        
        sdf <- do.call(rbind,lsdf)   
        orbit.plot <- ggplot(orbit.data, aes(x=OpM, y=INVT, colour=Subj, label=Year)) +
          geom_point(size=5, shape=20) + 
          geom_point(data=orbit.data,size=7, shape=20,color="black") + 
          geom_path(size=1) +
          geom_path(data=sdf,aes(x=x,y=y,label="",color=Subj),size=1) + 
          ggtitle("Title Orbits") +
          geom_text(data=subset(orbit.data,Year==2006 | Year==2014), 
                    aes(label=Year, vjust=1, hjust=1)) +
          theme(panel.background = element_rect(fill = 'white', colour = 'red'),  
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          geom_vline(xintercept=0, size=1) +
          geom_hline(yintercept=7, size=1) +
          scale_y_continuous(limits = c(7, 15), breaks=seq(7,15,1/2))
        print(orbit.plot)
      }
    }
  }
}
