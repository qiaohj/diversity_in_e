library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

source("commonFuns/colors.r")

args = commandArgs(trailingOnly=TRUE)
group<-args[1]

if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(1850:2020)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_1850_2020_df.rda")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))


var_pair<-start_layer_df[which(start_layer_df$VAR=="pr"),]
var_pair1<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmax"),], by=c("GCM", "SSP", "Y"))
#var_pair2<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmin"),], by=c("GCM", "SSP", "Y"))
var_pair1<-var_pair1%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair1)<-c("GCM", "PR", "Y", "TEMP")
#var_pair2<-var_pair2%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
#colnames(var_pair2)<-c("GCM", "PR", "Y", "TEMP")
#var_pair<-bind_rows(var_pair1, var_pair2)
var_pair<-var_pair1
var_pair$PR_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$PR, var_pair$Y, sep="_")
var_pair$TEMP_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$TEMP, var_pair$Y, sep="_")
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
i=1

if (F){
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area<=0){
      next()
    }
    
    
    print(paste(i, nrow(df_list), item$sp))
    occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
    if (nrow(occ)>1000){
      asdf
    }else{
      next()
    }
  }
}
example_sp<-"Scinax_x-signatus"
target_folder<-sprintf("../../Objects/Niche_Models/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
all_v<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))
all_v$in_out<-F
all_v$in_out<-(between(all_v$PR, fit$range_PR_sd_min, fit$range_PR_sd_max))&
  (between(all_v$TEMP_MAX, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max))&
  (between(all_v$TEMP_MIN, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max))
  


range_PR<-c(fit$range_PR_sd_min, fit$range_PR_sd_max)
range_TEMP<-c(fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)
start_env_layers_sample<-start_env_layers[sample(nrow(start_env_layers), 1000),]
all_vs_sample<-all_v[sample(nrow(all_v), 10000),]
table(all_vs_sample$in_out)
all_vs_sample1<-all_v[which(in_out)]
all_vs_sample2<-all_v[which(!in_out)]

TEMP_MAX<-seq(from=25, 45, by=0.1)
curves<-dnorm(TEMP_MAX, mean=fit$mean_TEMP_MAX, sd=fit$sd_TEMP_MAX)
curves<-curves*2e3/max(curves)
TEMP_MAX_Curve<-data.frame(x=TEMP_MAX, y=curves)

TEMP_MIN<-seq(from=-5, 30, by=0.1)
curves<-dnorm(TEMP_MIN, mean=fit$mean_TEMP_MIN, sd=fit$sd_TEMP_MIN)
curves<-curves*2e3/max(curves)
TEMP_MIN_Curve<-data.frame(x=TEMP_MIN, y=curves)

PR<-seq(from=0, 5000, by=10)
curves<-dnorm(PR, mean=fit$mean_PR, sd=fit$sd_PR)
curves<-curves*20/max(curves)-10
PR_Curve<-data.frame(x=PR, y=curves)


p<-ggplot()+
  geom_point(data=start_env_layers_sample, aes(x=TEMP_MAX, y=PR), 
             color=colors_black[4], alpha=0.4, size=0.5)+
  geom_point(data=start_env_layers_sample, aes(x=TEMP_MIN, y=PR), 
             color=colors_black[4], alpha=0.4, size=0.5)+
  
  geom_point(data=all_vs_sample2, aes(x=TEMP_MAX, y=PR, color=in_out), size=0.5)+
  geom_point(data=all_vs_sample2, aes(x=TEMP_MIN, y=PR, color=in_out), size=0.5)+
  geom_point(data=all_vs_sample1, aes(x=TEMP_MAX, y=PR, color=in_out), size=0.5)+
  geom_point(data=all_vs_sample1, aes(x=TEMP_MIN, y=PR, color=in_out), size=0.5)+
  scale_color_manual(breaks = c(T, F), 
                    values=color_dispersal)+
  geom_hline(yintercept=range_PR, linetype="dashed", color = colors_black[7], size=1.5)+
  geom_vline(xintercept=range_TEMP, linetype="dashed", color = colors_black[7], size=1.5)+
  geom_vline(xintercept=fit$mean_TEMP_MAX, linetype="dashed", color = colors_purple[5], size=1)+
  geom_vline(xintercept=fit$mean_TEMP_MAX-fit$sd_TEMP_MAX*3, linetype="dashed", color = colors_black[5], size=1)+
  geom_vline(xintercept=fit$mean_TEMP_MIN, linetype="dashed", color = colors_purple[5], size=1)+
  geom_vline(xintercept=fit$mean_TEMP_MIN+fit$sd_TEMP_MIN*3, linetype="dashed", color = colors_black[5], size=1)+
  geom_hline(yintercept=fit$mean_PR, linetype="dashed", color = colors_purple[5], size=1)+
  geom_path(data=TEMP_MAX_Curve, aes(x=x, y=y), color = colors_purple[5], size=1)+
  geom_path(data=TEMP_MIN_Curve, aes(x=x, y=y), color = colors_purple[5], size=1)+
  geom_path(data=PR_Curve, aes(x=y, y=x), color = colors_purple[5], size=1)+
  xlim(-10, 50)+
  ylim(0, 6000)+
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20))+
  xlab("Annual Maximum/Minimum Temperature")+ylab("Annual Precipitation")
p

ggsave(p, filename="../../Figures/Methods/niche.png", width=10, height=8)
ggsave(p, filename="../../Figures/Methods/niche.pdf", width=10, height=8)
  

