library(raster)
library(dplyr)
#library(alphahull)
library(concaveman)
library(sf)
setwd("Y:/Script/diversity_in_e")
source("colors.r")
group<-"Amphibians"
mask<-raster("../../Raster/mask_index.tif")
target_sp<-c()

df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
for (i in c(1:nrow(df_list))){
  item<-df_list[i,]
  example_sp<-gsub(" ", "_", item$sp)
  if (item$area<=0){
    next()
  }
  
  
  print(paste(i, nrow(df_list), item$sp))
  #example_sp<-"Adenomera_bokermanni"
  target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, example_sp)
  target<-sprintf("%s/dispersal", target_folder)
  dispersal_log<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1_1"))
  if (is.null(dispersal_log)){
    next()
  }
  start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
  if (is.null(start_dis)){
    next()
  }
  #start_dis<-start_dis%>%dplyr::filter(in_out==1)
  start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
  colnames(start_dis)<-c("x", "y")
  start_dis$mask_index<-extract(mask, start_dis)
  
  dispersal_log_end<-dispersal_log%>%dplyr::filter(YEAR==2100)
  dispersal_log_others<-dispersal_log%>%dplyr::filter((YEAR!=2100)&!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))%>%
    dplyr::distinct(x, y, mask_index, YEAR)
  
  #points<-st_sfc(st_point(as.matrix(start_dis[, c("x", "y")])))
  points<-as.matrix(bind_rows(start_dis[, c("x", "y")], dispersal_log_end[, c("x", "y")]))
  alpha_hull<-concaveman(points)
  dispersal_log_others$inout<-point.in.polygon(dispersal_log_others$x, dispersal_log_others$y,
                                               alpha_hull[,1], alpha_hull[,2])
  dispersal_log_others<-dispersal_log_others%>%
    dplyr::filter(!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))
  
  
  p<-ggplot()+
    geom_raster(data=start_dis, aes(x=x, y=y), fill=colors_blue[8])+
    geom_raster(data=dispersal_log_end, aes(x=x, y=y), fill=colors_red[8])+
    geom_raster(data=dispersal_log_others%>%dplyr::filter(inout==1), aes(x=x, y=y), fill=colors_green[8])+
    geom_path(data=as.data.frame(alpha_hull), aes(x=V1, y=V2), color="grey")
  print(p)
  xx<-readline(prompt = "X=EXIT:")
  if (toupper(xx)=="X"){
    break()
  }
  if (xx=="1"){
    target_sp<-c(target_sp, example_sp)    
  }
}

target2<-c()
#for (i in c(1:length(target_sp))){
target2<-c("Acris_gryllus",                "Adelphobates_castaneoticus",   "Adelphobates_quinquevittatus", "Adenomera_diptyx"  ,          
           "Afrixalus_weidholzi",          "Allobates_fuscellus",          "Ambystoma_macrodactylum",      "Plethodon_dorsalis",          
           "Phrynomantis_affinis" )
ps<-list()
for (i in c(1:length(target2))){
  example_sp<-target2[i]
  if (item$area<=0){
    next()
  }
  
  
  print(paste(i, length(target_sp), example_sp))
  #example_sp<-"Adenomera_bokermanni"
  target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, example_sp)
  target<-sprintf("%s/dispersal", target_folder)
  dispersal_log<-readRDS(sprintf("%s/%s.rda", target, "UKESM1_SSP585_1_1"))
  if (is.null(dispersal_log)){
    next()
  }
  start_dis<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
  if (is.null(start_dis)){
    next()
  }
  #start_dis<-start_dis%>%dplyr::filter(in_out==1)
  start_dis<-start_dis%>%ungroup()%>%dplyr::distinct(x, y)
  colnames(start_dis)<-c("x", "y")
  start_dis$mask_index<-extract(mask, start_dis)
  
  dispersal_log_end<-dispersal_log%>%dplyr::filter(YEAR==2100)
  dispersal_log_others<-dispersal_log%>%dplyr::filter((YEAR!=2100)&!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))%>%
    dplyr::distinct(x, y, mask_index, YEAR)
  
  #points<-st_sfc(st_point(as.matrix(start_dis[, c("x", "y")])))
  points<-as.matrix(bind_rows(start_dis[, c("x", "y")], dispersal_log_end[, c("x", "y")]))
  alpha_hull<-concaveman(points)
  dispersal_log_others$inout<-point.in.polygon(dispersal_log_others$x, dispersal_log_others$y,
                                               alpha_hull[,1], alpha_hull[,2])
  dispersal_log_others<-dispersal_log_others%>%
    dplyr::filter(!(mask_index %in% c(start_dis$mask_index, dispersal_log_end$mask_index)))
  
  
  
  p<-ggplot()+
    geom_tile(data=start_dis, aes(x=x, y=y), fill=colors_blue[8])+
    geom_tile(data=dispersal_log_end, aes(x=x, y=y), fill=colors_red[8])+
    geom_tile(data=dispersal_log_others%>%dplyr::filter(inout==1), aes(x=x, y=y), fill=colors_green[8])+
    geom_tile(data=dispersal_log_others%>%dplyr::filter(inout==0), aes(x=x, y=y), fill=colors_black[3])+
    geom_path(data=as.data.frame(alpha_hull), aes(x=V1, y=V2), color=colors_black[9], size=1)+theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  ps[[length(ps)+1]]<-p

}

p<-ggarrange(ps[[9]], ggarrange(ps[[2]],ps[[3]],ps[[4]],ps[[5]],ps[[6]],ps[[7]],ps[[8]],ps[[1]],
                             labels=paste("(", c(2:9), ")", sep=""), ncol=4, nrow=2),
          labels=c("(1)", ""), nrow=2, heights=c(12,8))
ggsave(p, filename="../../Figures/Methods/dispersal_usage.png", width=10, height=10)
