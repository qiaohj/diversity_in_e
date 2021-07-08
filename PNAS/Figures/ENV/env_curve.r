library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

if (F){
  GCMs<-c("UKESM1")
  SSPs<-c("SSP119", "SSP245", "SSP585")
  VARs<-c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")
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
df<-readRDS("../../Objects_Full_species/mean_env_year.rda")
df[which(df$VAR!="pr"),]$V<-df[which(df$VAR!="pr"),]$V*10-273.16
df_se_2020<-df%>%dplyr::filter(Y<=2020)%>%
  dplyr::group_by(Y, GCM, VAR)%>%
  dplyr::summarise(annul_prec=mean(V),
                   annul_max_temp=mean(V),
                   annual_min_temp=mean(V))
df_se_2020$SSP<-"SSP_ALL"
df_se_2100<-df%>%dplyr::filter(Y>2020)%>%
  dplyr::group_by(Y, GCM, SSP, VAR)%>%
  dplyr::summarise(annul_prec=sum(V),
                   annul_max_temp=mean(V),
                   annual_min_temp=mean(V))
df_se<-bind_rows(df_se_2020, df_se_2100)
p1<-ggplot(df_se%>%filter((VAR=="pr")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_prec, color=SSP, linetype=GCM))+
  geom_line()+
  geom_line(data=df_se%>%filter((VAR=="pr")&(between(Y, 1850, 2020))), aes(x=Y, y=annul_prec, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  scale_x_continuous(breaks = c(seq(1850, 2000, by=50), 2020, seq(2050, 2100, by=50)))+
  xlab("Year")+
  ylab("Annual precipitation")+
  labs(linetype="ESM")+
  theme_bw()
legend_g<-get_legend(p1)
p1<-p1+
  theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank())

p2<-ggplot(df_se%>%filter((VAR=="tasmax")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_max_temp, color=SSP, linetype=GCM))+geom_line()+
  geom_line(data=df_se%>%filter((VAR=="tasmax")&(between(Y, 1850, 2020))), aes(x=Y, y=annul_max_temp, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  scale_x_continuous(breaks = c(seq(1850, 2000, by=50), 2020, seq(2050, 2100, by=50)))+
  scale_y_continuous(labels = paste("   ", seq(18, 24, by=2)), breaks=seq(18, 24, by=2))+
  xlab("Year")+
  ylab("Annual maximum temperature")+
  theme_bw()+
  theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank())

p2

p3<-ggplot(df_se%>%filter((VAR=="tasmin")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_max_temp, color=SSP, linetype=GCM))+
  geom_line()+
  geom_line(data=df_se%>%filter((VAR=="tasmin")&(between(Y, 1850, 2020))), aes(x=Y, y=annul_max_temp, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  scale_x_continuous(breaks = c(seq(1850, 2000, by=50), 2020, seq(2050, 2100, by=50)))+
  scale_y_continuous(labels=paste("  ", seq(-12, -4, by=4)), breaks=seq(-12, -4, by=4))+
  xlab("Year")+
  ylab("Annual minimum temperature")+
  theme_bw()+
  theme(legend.position="none")

p3

p<-ggarrange(p1, p2, p3, nrow=3, 
             common.legend = T, 
             legend.grob=legend_g, 
             legend="right", labels = c("(a)", "(b)", "(c)"),
             label.x=0.065, label.y=0.95)
p
ggsave(p, file="../../Figures_Full_species/Env/GCM_Curves.png", width=10, height=8)
ggsave(p, file="../../Figures_Full_species/Env/GCM_Curves.pdf", width=10, height=8)
