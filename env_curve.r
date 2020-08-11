library(raster)
library(ggplot2)
library(dplyr)
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
df_se<-df%>%dplyr::group_by(Y, GCM,  SSP, VAR, SUB, Label)%>%dplyr::summarise(annul_prec=sum(v),
                                                annul_max_temp=max(v),
                                                annual_min_temp=min(v))
p<-ggplot(df_se%>%filter((VAR=="pr")&(SUB=="sum")), aes(x=Y, y=annul_prec, color=factor(Label)))+geom_line()+
  theme_bw()
ggsave(p, file="../../Figures/Env/GCM_Curves_prec.png", width=6, height=4)
p<-ggplot(df_se%>%filter((VAR=="tasmax")&(SUB=="max")), aes(x=Y, y=annul_max_temp, color=factor(Label)))+geom_line()+
  geom_line(data=df_se%>%filter((VAR=="tasmin")&(SUB=="min")), aes(x=Y, y=annual_min_temp, color=factor(Label)))+
  theme_bw()
ggsave(p, file="../../Figures/Env/GCM_Curves_temp.png", width=6, height=4)

p<-ggplot(df%>%filter(between(Y, 2010, 2020)), aes(x=DATE, y=v, color=factor(Label)))+geom_line()+
  facet_wrap( ~ VAR+SUB, ncol=1, scales = 'free')+
  theme_bw()
ggsave(p, file="../../Figures/Env/GCM_Curves_2010_2020.png", width=6, height=6)
