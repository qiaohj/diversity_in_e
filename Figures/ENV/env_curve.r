library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

if (F){
  GCMs<-c("UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  VARs<-c("pr", "tasmax", "tasmin")
  SUBs<-c("mean", "sum", "min", "max")
  years<-c(1850:2100)
  months<-c(1:12)
  coms<-expand.grid(GCM=GCMs, SSP=SSPs, VAR=VARs, SUB=SUBs, Y=years, M=months)
  base<-"../../Raster"
  i=1
  df<-NULL
  for (i in c(1:nrow(coms))){
    s<-coms[i,]
    r_file<-sprintf("%s/%s/%s/%s/%d/%s_%d.tif", base, s$GCM, s$SSP, s$VAR, s$Y, s$SUB, s$M)
    
    if (!file.exists(r_file)){
      next()
    }
    print(r_file)
    r<-raster(r_file)
    s$v<-mean(values(r), na.rm=T)
    if (is.null(df)){
      df<-s
    }else{
      df<-bind_rows(df, s)
    }
  }
  df$DATE<-as.Date(paste(df$Y, df$M, 1, sep="/"), format="%Y/%m/%d")
  df$Label<-paste(df$GCM, df$SSP)
  
  saveRDS(df, "../../Objects/mean_env_year.rda")
}
source("commonFuns/colors.r")
source("commonFuns/functions.r")
df<-readRDS("../../Objects/mean_env_year.rda")
df_se_2020<-df%>%dplyr::filter(Y<=2020)%>%
  dplyr::group_by(Y, GCM, VAR)%>%
  dplyr::summarise(annul_prec=mean(V),
                   annul_max_temp=max(V),
                   annual_min_temp=min(V))
df_se_2020$SSP<-"SSP_ALL"
df_se_2100<-df%>%dplyr::filter(Y>2020)%>%
  dplyr::group_by(Y, GCM, SSP, VAR)%>%
  dplyr::summarise(annul_prec=sum(V),
                   annul_max_temp=max(V),
                   annual_min_temp=min(V))
df_se<-bind_rows(df_se_2020, df_se_2100)
p1<-ggplot(df_se%>%filter((VAR=="pr")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_prec, color=SSP, linetype=GCM))+geom_line()+
  geom_line(data=df_se%>%filter((VAR=="pr")&(between(Y, 1850, 2020))), aes(x=Y, y=annul_prec, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  scale_x_continuous(breaks = c(seq(1850, 2000, by=50), 2020, seq(2050, 2100, by=50)))+
  xlab("Year")+
  ylab("Annual Precipitation")+
  theme_bw()
legend_g<-get_legend(p1)
p1<-p1+theme(legend.position="none")

p2<-ggplot(df_se%>%filter((VAR=="tasmax")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_max_temp, color=SSP, linetype=GCM))+geom_line()+
  geom_line(data=df_se%>%filter((VAR=="tasmax")&(between(Y, 1850, 2020))), aes(x=Y, y=annul_max_temp, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  scale_x_continuous(breaks = c(seq(1850, 2000, by=50), 2020, seq(2050, 2100, by=50)))+
  xlab("Year")+
  ylab("Annual Maximum Temperature")+
  theme_bw()+theme(legend.position="none")

p2
p<-ggarrange(p1, p2, nrow=2, common.legend = T, legend.grob=legend_g, legend="right")
p
ggsave(p, file="../../Figures/Env/GCM_Curves.png", width=15, height=8)
