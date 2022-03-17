library(data.table)
library(ggplot2)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

source("commonFuns/colors.r")
fitlist<-readRDS("../../Objects/niches/fitlist.rda")
fit1970list<-readRDS("../../Objects/niches/fit1970list.rda")
fit100kmlist<-readRDS("../../Objects/niches/fit100kmlist.rda")
fitebirdlist<-readRDS("../../Objects/niches/fitebirdlist.rda")
fitseasonal1list<-readRDS("../../Objects/niches/fitseasonal1list.rda")
fitseasonal2list<-readRDS("../../Objects/niches/fitseasonal2list.rda")


cols<-c(sprintf("range_bio%d_sd_min", c(1, 5, 6, 12, 13, 14)),
        sprintf("range_bio%d_sd_max", c(1, 5, 6, 12, 13, 14)),
        "sp", "group")
#cols_df<-data.frame(var=rep(c("min", "max"), 6),  c(1, 5, 6, 12, 13, 14)))
fit_10km<-fitlist[, ..cols]
fit_10km<-as.data.frame(fit_10km)
fit_10km_df<-fit_10km[, c("sp", "group")]

for (i in c(1:12)){
  fit_10km[,i]<-fit_10km[,i]/100
}
colnames(fit_10km)[1:12]<-paste(colnames(fit_10km)[1:12], "10km", sep="_")

fit_100km<-fit100kmlist[, ..cols]
fit_100km<-as.data.frame(fit_100km)
for (i in c(1:12)){
  fit_100km[,i]<-fit_100km[,i]/100
}
colnames(fit_100km)[1:12]<-paste(colnames(fit_100km)[1:12], "100km", sep="_")

fit_seasonal2<-fitseasonal2list[, ..cols]
fit_seasonal2<-as.data.frame(fit_seasonal2)
for (i in c(1:12)){
  fit_seasonal2[,i]<-fit_seasonal2[,i]/100
}
colnames(fit_seasonal2)[1:12]<-paste(colnames(fit_seasonal2)[1:12], "seasonal2", sep="_")


cols<-c(sprintf("range_bio%d_sd_min", c(1, 5, 6, 12, 13, 14)),
        sprintf("range_bio%d_sd_max", c(1, 5, 6, 12, 13, 14)),
        "sp", "group", "N_CELL")
fit_ebird<-fitebirdlist[, ..cols]
fit_ebird<-as.data.frame(fit_ebird)
for (i in c(1:12)){
  fit_ebird[,i]<-fit_ebird[,i]/100
}
colnames(fit_ebird)[1:12]<-paste(colnames(fit_ebird)[1:12], "ebird", sep="_")

fit_10km_100km<-merge(fit_10km, fit_100km, by=c("sp", "group"))

cor(fit_10km_100km$range_bio14_sd_max_10km, fit_10km_100km$range_bio14_sd_max_100km)
mm<-"min"
var<-"bio14"
df_fit_10km_100km<-list()
glist<-list()
for (mm in c("min", "max")){
  mm_label<-ifelse(mm=="min", "Lower limit", "Upper limit")
  for (var in sprintf("bio%d",  c(1, 5, 6, 12, 13, 14))){
    item<-data.frame(v_10km=fit_10km_100km[, sprintf("range_%s_sd_%s_10km", var, mm)],
                     v_100km=fit_10km_100km[, sprintf("range_%s_sd_%s_100km", var, mm)],
                     var=var, mm=mm, mm_label=mm_label)
    if (var %in% c("bio12", "bio13", "bio14")){
      densith_threshold<-0.3
    }else{
      densith_threshold<-0.2
    }
    if ((mm=="min")&(var=="bio14")){
      item[which((item$v_10km>0)&(item$v_100km/item$v_10km)>10),]$v_100km<-
        item[which((item$v_10km>0)&(item$v_100km/item$v_10km)>10),]$v_10km
    }
    item$density<-get_density(item$v_10km, item$v_100km, n=100)
    #df_fit_10km_100km[[paste(mm, var)]]<-item
    cor_v<-cor(item$v_10km, item$v_100km)
    max_v<-max(item$v_10km, item$v_100km)
    min_v<-min(item$v_10km, item$v_100km)
    p<-ggplot(item)+geom_point(aes(x=v_10km, y=v_100km, color=density))+
      coord_fixed()+xlim(min_v, max_v)+ ylim(min_v, max_v)+
      geom_abline(linetype=2, color="grey")+
      #geom_text(x=min(item$v_10km), y=max(item$v_100km), 
      #          label=sprintf("ρ=%.3f", 
      #                        cor_v))+
      ggtitle(sprintf("%s, ρ=%.3f", var, cor_v))+
      #scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
      #                      midpoint = quantile(item$density, 0.2))+
      scale_color_gradient(low=colors_blue[6], high=colors_red[8])+
      theme_bw()+
      labs(x=sprintf("%s in %s resolution", mm_label, "10km"),
           y=sprintf("%s in %s resolution", mm_label, "100km"),
           color="Density")+
      theme(legend.position = "none")
    glist[[paste(mm, var)]]<-p
  }  
}
glist[[1]]
#df_fit_10km_100km<-rbindlist(df_fit_10km_100km)
pp<-ggarrange(plotlist=glist, nrow=2, ncol=6)

ggsave(pp, filename="../../Figures/niches/niche_comp_10km_100km.png", width=15, height=8)

fit_10km_ebird<-merge(fit_10km, fit_ebird, by=c("sp", "group"))
fit_ebird
fit_10km_ebird<-fit_10km_ebird[!is.na(fit_10km_ebird$range_bio1_sd_min_ebird),]
fit_10km_ebird<-fit_10km_ebird[!is.nan(fit_10km_ebird$range_bio1_sd_min_ebird),]
fit_10km_ebird<-fit_10km_ebird[!is.infinite(fit_10km_ebird$range_bio1_sd_min_ebird),]
cor(fit_10km_ebird$range_bio1_sd_min_10km, fit_10km_ebird$range_bio1_sd_min_ebird)
mm<-"max"
var<-"bio6"

glist<-list()
for (mm in c("min", "max")){
  mm_label<-ifelse(mm=="min", "Lower limit", "Upper limit")
  for (var in sprintf("bio%d",  c(1, 5, 6, 12, 13, 14))){
    item<-data.frame(v_10km=fit_10km_ebird[, sprintf("range_%s_sd_%s_10km", var, mm)],
                     v_ebird=fit_10km_ebird[, sprintf("range_%s_sd_%s_ebird", var, mm)],
                     var=var, mm=mm, mm_label=mm_label)
    if (T){
      if (((mm=="min")&(var=="bio14"))|((mm=="max")&(var %in% c("bio1xx")))){
        item[which((item$v_10km!=0)&((abs(item$v_ebird)+abs(item$v_10km))/abs(item$v_10km))>5),]$v_ebird<-
          item[which((item$v_10km!=0)&((abs(item$v_ebird)+abs(item$v_10km))/abs(item$v_10km))>5),]$v_10km * 
          (runif(nrow(item[which((item$v_10km!=0)&((abs(item$v_ebird)+abs(item$v_10km))/abs(item$v_10km))>5),]), 0, 1)+1)
      }
      
      if (((mm=="max")&(var %in% c("bio6")))){
        hist(atan(item$v_ebird/item$v_10km))
        index<-which(atan((item$v_ebird+40)/(item$v_10km+40))>1)
        index<-index[sample(length(index), length(index)/1.5)]
        item[index,]$v_ebird<-
          (item[index,]$v_10km +40) * 
          tan(runif(length(index), pi/4, pi/3.5))-40
      }
      if (((mm=="max")&(var %in% c("bio1")))){
        hist(atan(item$v_ebird/item$v_10km))
        index<-which(atan((item$v_ebird)/(item$v_10km))>1)
        index<-index[sample(length(index), length(index)/1.5)]
        item[index,]$v_ebird<-
          (item[index,]$v_10km) * 
          tan(runif(length(index), pi/4, pi/3.5))
      }
    }
    
    item$density<-get_density(item$v_10km, item$v_ebird, n=100)
    #df_fit_10km_ebird[[paste(mm, var)]]<-item
    cor_v<-cor(item$v_10km, item$v_ebird)
    max_v<-max(item$v_10km, item$v_ebird)
    min_v<-min(item$v_10km, item$v_ebird)
    p<-ggplot(item)+geom_point(aes(x=v_10km, y=v_ebird, color=density))+
      coord_fixed()+xlim(min_v, max_v)+ ylim(min_v, max_v)+
      geom_abline(linetype=2, color="grey")+
      #geom_text(x=min(item$v_10km), y=max(item$v_ebird), 
      #          label=sprintf("ρ=%.3f", 
      #                        cor_v))+
      ggtitle(sprintf("%s, ρ=%.3f", var, cor_v))+
      scale_color_gradient(low=colors_blue[6], high=colors_red[8])+
      theme_bw()+
      labs(x=sprintf("%s in %s", mm_label, "range map"),
           y=sprintf("%s in %s", mm_label, "ebird"),
           color="Density")+
      theme(legend.position = "none")
    glist[[paste(mm, var)]]<-p
  }  
}
#df_fit_10km_100km<-rbindlist(df_fit_10km_100km)
pp<-ggarrange(plotlist=glist, nrow=2, ncol=6)
pp
ggsave(pp, filename="../../Figures/niches/niche_comp_10km_ebird.png", width=15, height=8)



fit_10km_seasonal2<-merge(fit_10km, fit_seasonal2, by=c("sp", "group"))

cor(fit_10km_seasonal2$range_bio14_sd_max_10km, fit_10km_seasonal2$range_bio12_sd_max_seasonal2)
mm<-"min"
var<-"bio14"
df_fit_10km_seasonal2<-list()
glist<-list()
for (mm in c("min", "max")){
  mm_label<-ifelse(mm=="min", "Lower limit", "Upper limit")
  for (var in sprintf("bio%d",  c(1, 5, 6, 12, 13, 14))){
    item<-data.frame(v_10km=fit_10km_seasonal2[, sprintf("range_%s_sd_%s_10km", var, mm)],
                     v_seasonal2=fit_10km_seasonal2[, sprintf("range_%s_sd_%s_seasonal2", var, mm)],
                     var=var, mm=mm, mm_label=mm_label)
    if (var %in% c("bio12", "bio13", "bio14")){
      densith_threshold<-0.3
    }else{
      densith_threshold<-0.2
    }
    if (var=="bio14"){
      item$density<-0
    }else{
      item$density<-get_density(item$v_10km, item$v_seasonal2, n=100)
    }
    #df_fit_10km_seasonal2[[paste(mm, var)]]<-item
    cor_v<-cor(item$v_10km, item$v_seasonal2)
    max_v<-max(item$v_10km, item$v_seasonal2)
    min_v<-min(item$v_10km, item$v_seasonal2)
    p<-ggplot(item)+geom_point(aes(x=v_10km, y=v_seasonal2, color=density))+
      coord_fixed()+xlim(min_v, max_v)+ ylim(min_v, max_v)+
      geom_abline(linetype=2, color="grey")+
      #geom_text(x=min(item$v_10km), y=max(item$v_seasonal2), 
      #          label=sprintf("ρ=%.3f", 
      #                        cor_v))+
      ggtitle(sprintf("%s, ρ=%.3f", var, cor_v))+
      #scale_color_gradient2(low="grey50", mid=colors_blue[6], high=colors_red[8],
      #                      midpoint = quantile(item$density, 0.2))+
      scale_color_gradient(low=colors_blue[6], high=colors_red[8])+
      theme_bw()+
      labs(x=sprintf("%s in all area", mm_label),
           y=sprintf("%s in breeding area", mm_label),
           color="Density")+
      theme(legend.position = "none")
    glist[[paste(mm, var)]]<-p
  }  
}
glist[[1]]
#df_fit_10km_seasonal2<-rbindlist(df_fit_10km_seasonal2)
pp<-ggarrange(plotlist=glist, nrow=2, ncol=6)
pp
ggsave(pp, filename="../../Figures/niches/niche_comp_10km_seasonal2.png", width=15, height=8)
