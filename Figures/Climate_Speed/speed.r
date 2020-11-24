setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
library(dplyr)
library(ggplot2)
if (F){
  library("DBI")
  library("mgcv")
  library("pastecs")
  mydb <- dbConnect(RSQLite::SQLite(), sprintf("%s/ISEA3H8/SQLITE/env_Hadley3D.sqlite", base))  
  max_prec<-dbReadTable(mydb, "Debiased_Maximum_Monthly_Precipitation")
  max_temp<-dbReadTable(mydb, "Debiased_Maximum_Monthly_Temperature")
  min_temp<-dbReadTable(mydb, "Debiased_Minimum_Monthly_Temperature")
  dbDisconnect(mydb) 
  head(max_prec)
  max_prec_se<-max_prec%>%dplyr::group_by(year)%>%
    dplyr::summarise(mean=mean(v),
                     sd=sd(v),
                     ci=CI(v)[2]-CI(v)[3])
  max_prec_se$year<-max_prec_se$year*-1
  max_prec_se$type<-"Maximum Monthly Precipitation"
  
  max_temp_se<-max_temp%>%dplyr::group_by(year)%>%
    dplyr::summarise(mean=mean(v),
                     sd=sd(v),
                     ci=CI(v)[2]-CI(v)[3])
  max_temp_se$type<-"Maximum Monthly Temperature"
  
  max_temp_se$year<-max_temp_se$year*-1
  plot(max_temp_se$year, max_temp_se$mean, type="l", col="red")
  
  min_temp_se<-min_temp%>%dplyr::group_by(year)%>%
    dplyr::summarise(mean=mean(v),
                     sd=sd(v),
                     ci=CI(v)[2]-CI(v)[3])
  min_temp_se$type<-"Minimum Monthly Temperature"
  min_temp_se$year<-min_temp_se$year*-1
  lines(min_temp_se$year, min_temp_se$mean, col="blue")
  
  env_se<-bind_rows(bind_rows(max_prec_se, max_temp_se), min_temp_se)
  
  env_df_new<-NULL
  for (var in c("Maximum Monthly Precipitation",
                "Maximum Monthly Temperature",
                "Minimum Monthly Temperature")){
    print(var)
    if (var=="Maximum Monthly Precipitation"){
      k<-40
    }else{
      k<-20
    }
    if (var=="Maximum Monthly Temperature"){
      type_first<-1
    }else{
      type_first<--1
    } 
    env_df<-env_se%>%dplyr::filter(type==var)
    
    max_temp_gam<-gam(mean~s(year, k=k),data=env_df)
    #plot(max_temp_gam)
    
    env_df$predicted<-predict(max_temp_gam, env_df)
    
    v<-as.vector(pull(env_df, predicted))
    tp<-turnpoints(v)
    env_df$peak<-tp$peaks
    env_df$pit<-tp$pits
    env_df$direction<-0
    env_df$group<-0
    env_df$length<-0
    direction<-0
    length<-0
    group<-0
    for (year in c(-1200:0)){
      if (pull(env_df[which(env_df$year==year), "peak"])){
        env_df[which(env_df$group==group), "length"]<-length
        direction<- -1
        length<-0
        group<-group+1
      }
      if (pull(env_df[which(env_df$year==year), "pit"])){
        env_df[which(env_df$group==group), "length"]<-length
        direction<- 1
        length<-0
        group<-group+1
      }
      length<-length+1
      env_df[which(env_df$year==year), "direction"]<- direction
      env_df[which(env_df$year==year), "group"]<- group
    }
    env_df[which(env_df$group==group), "length"]<-length
    unique(env_df[, c("group", "length")])
    if (type_first==1){
      env_df[which((env_df$direction==0)&(env_df$year==-1200)), "peak"]<-T
    }else{
      env_df[which((env_df$direction==0)&(env_df$year==-1200)), "pit"]<-T
    }
    env_df[which(env_df$direction==0), "direction"]<-type_first
    if (is.null(env_df_new)){
      env_df_new<-env_df
    }else{
      env_df_new<-bind_rows(env_df_new, env_df)
    }
    ggplot(env_df, aes(x=year, y=predicted, color=factor(direction)))+
      geom_point(aes(y=mean), color="black")+
      geom_point()
  }
  saveRDS(env_df_new, sprintf("%s/Data/env_se.rda", base))
}

library(dplyr)
library(ggplot2)
library(hrbrthemes)
source("commonFuns/colors.r")
env_se<-readRDS("../../Figures/Climate_Speed/env_se.rda")
env_se_item<-env_se%>%dplyr::group_by(group, type, direction)%>%dplyr::summarise(
  mean_year=mean(year),
  max=max(mean),
  min=min(mean),
  length=mean(length)
)
env_se_item$speed<-(env_se_item$max-env_se_item$min)/env_se_item$length
env_se_item_se<-env_se_item%>%ungroup()%>%dplyr::group_by(type, direction, mean_year)%>%
  dplyr::summarise(mean_speed=mean(speed),
                   sd_speed=sd(speed))
write.csv(env_se_item, "../../Figures/Climate_Speed/Env_speed.csv")


min_prec<-6
scale<-23
offset<-35


colors_up_down<-c("-1"=colors_blue[7], "1"=colors_red[7])

env_df<-env_se%>%dplyr::filter(type=="Maximum Monthly Temperature")
env_df[which(env_df$year==0), "peak"]<-T
env_df[which(env_df$year==-1200), "peak"]<-F
env_df[which(env_df$year==-1200), "pit"]<-T
env_point<-env_df%>%dplyr::filter((peak)|(pit))
env_point$hjust_v<--0.1
env_point[which(env_point$peak), "hjust_v"]<--0.1
env_point$vjust_v<--1
env_point[which(env_point$pit), "vjust_v"]<-1


env_se_item_se_item<-env_se_item_se%>%dplyr::filter(type=="Maximum Monthly Temperature")
p<-ggplot() +
  geom_path(data=env_df,
            aes(x=year, y=predicted), size=2, color=colors_black[6]) + 
  geom_point(data=env_df,
             aes(x=year, y=predicted, group=type, color=as.character(direction)), size=0.3) + 
  geom_vline(data=env_point,
             aes(xintercept=year), linetype=2)+
  geom_text(data=env_se_item_se_item, 
            aes(x=mean_year, y=20, label=round(mean_speed, digits=3), color=as.character(direction)))+
  geom_text(data=env_point, aes(x=year, y=predicted, label=round(predicted, digits=2), 
                                hjust=hjust_v, vjust=vjust_v))+
  ylim(c(19.5, 31.5))+
  scale_x_continuous(breaks=env_point$year, labels=env_point$year/10)+
  xlab("Year (kbp)")+
  ylab("Annual Maximum Temperature")+
  scale_color_manual(values=colors_up_down)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Temperature (120kbp ~ present)")
ggsave(p, filename="../../Figures/Climate_Speed/temp_past.png", width=10, height=5)
ggsave(p, filename="../../Figures/Climate_Speed/temp_past.pdf", width=10, height=5)


env_df<-env_se%>%dplyr::filter(type=="Maximum Monthly Precipitation")
env_df[which(env_df$year==0), "peak"]<-F
env_df[which(env_df$year==0), "pit"]<-T
env_df[which(env_df$year==-1200), "peak"]<-F
env_df[which(env_df$year==-1200), "pit"]<-T
env_point<-env_df%>%dplyr::filter((peak)|(pit))
env_point$hjust_v<--0.1
env_point[which(env_point$peak), "hjust_v"]<--0.1
env_point$vjust_v<--1
env_point[which(env_point$pit), "vjust_v"]<-1


env_se_item_se_item<-env_se_item_se%>%dplyr::filter(type=="Maximum Monthly Precipitation")
p<-ggplot() +
  geom_path(data=env_df,
            aes(x=year, y=predicted), size=2, color=colors_black[6]) + 
  geom_point(data=env_df,
             aes(x=year, y=predicted, group=type, color=as.character(direction)), size=0.3) + 
  geom_vline(data=env_point,
             aes(xintercept=year), linetype=2)+
  geom_text(data=env_se_item_se_item, 
            aes(x=mean_year, y=6.9, label=round(mean_speed, digits=3), color=as.character(direction)))+
  geom_text(data=env_point, aes(x=year, y=predicted, label=round(predicted, digits=2), 
                                hjust=hjust_v, vjust=vjust_v))+
  ylim(c(6.8, 9))+
  scale_x_continuous(breaks=env_point$year, labels=abs(env_point$year)/10)+
  xlab("Year (kbp)")+
  ylab("Maximum Monthly Precipitation")+
  scale_color_manual(values=colors_up_down)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Precipitation (120kbp ~ present)")
ggsave(p, filename="../../Figures/Climate_Speed/prec_past.png", width=10, height=5)

df<-readRDS("../../Objects/mean_env_year.rda")
df_se<-df%>%dplyr::group_by(Y, GCM,  SSP, VAR)%>%dplyr::summarise(annul_prec=sum(V),
                                                                  annul_max_temp=max(V),
                                                                  annual_min_temp=min(V))
df_se_2020<-df_se%>%dplyr::filter(Y==2020)
colnames(df_se_2020)<-c("Y_2020", "GCM", "SSP", "VAR", "prec_2020", "max_temp_2020", "min_temp_2020")
df_se_1850<-df_se%>%dplyr::filter(Y==1850)
colnames(df_se_1850)<-c("Y_1850", "GCM", "SSP", "VAR", "prec_1850", "max_temp_1850", "min_temp_1850")

df_se_speed<-left_join(df_se, df_se_2020, by=c("GCM", "SSP", "VAR"))
df_se_speed<-left_join(df_se_speed, df_se_1850, by=c("GCM", "SSP", "VAR"))

df_se_speed$max_temp_speed_2020<-(df_se_speed$annul_max_temp-df_se_speed$max_temp_2020)*100/(df_se_speed$Y-df_se_speed$Y_2020)
df_se_speed$max_temp_speed_1850<-(df_se_speed$annul_max_temp-df_se_speed$max_temp_1850)*100/(df_se_speed$Y-df_se_speed$Y_1850)

p2<-ggplot(df_se%>%filter((VAR=="tasmax")&(between(Y, 2021, 2100))), aes(x=Y, y=annul_max_temp, color=SSP, linetype=GCM))+
  geom_line()+
  geom_line(data=df_se%>%filter((VAR=="tasmax")&(between(Y, 1850, 2020))), 
            aes(x=Y, y=annul_max_temp, linetype=GCM),  color=colors_black[5])+
  geom_vline(xintercept = 2020, color=colors_black[6], linetype=3)+
  geom_text(data=df_se_speed%>%dplyr::filter((VAR=="tasmax")&(Y==2020)&(SSP=="SSP119")), 
            aes(x=Y, y=round(annul_max_temp, 2), 
                label=paste(round(annul_max_temp, 2), " Speed=", round(max_temp_speed_1850, 3), sep="")),
                color=colors_black[9])+
  geom_text(data=df_se_speed%>%dplyr::filter((VAR=="tasmax")&(Y==1850)&(SSP=="SSP119")), 
            aes(x=Y, y=round(annul_max_temp, 2), 
                label=round(annul_max_temp, 2)),
            color=colors_black[9])+
  geom_text(data=df_se_speed%>%dplyr::filter((VAR=="tasmax")&(Y==2100)), 
            aes(x=Y, y=round(annul_max_temp, 2), 
                label=paste(round(annul_max_temp, 2), " Speed=", round(max_temp_speed_2020, 3), sep=""),
                color=SSP, hjust=-0.1))+
  xlim(c(1850, 2160))+
  scale_color_manual(values=color_ssp)+
  scale_linetype_manual(values=linetype_gcm)+
  xlab("Year")+
  ylab("Annual Maximum Temperature")+
  theme_bw()

ggsave(p2, filename="../../Figures/Climate_Speed/temp_future.png", width=10, height=5)
ggsave(p2, filename="../../Figures/Climate_Speed/temp_future.pdf", width=10, height=5)
write.csv(df_se_speed, "../../Figures/Climate_Speed/speed_future.csv")
