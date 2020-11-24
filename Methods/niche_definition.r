library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)

setwd("Y:/Script/diversity_in_e")

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
  (between(all_v$TEMP, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max))
range_PR<-c(fit$range_PR_sd_min, fit$range_PR_sd_max)
range_TEMP<-c(fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)
start_env_layers_sample<-start_env_layers[sample(nrow(start_env_layers), 5000),]
all_vs_sample<-all_v[sample(nrow(all_v), 1000),]
p<-ggplot()+
  geom_point(data=start_env_layers_sample, aes(x=TEMP, y=PR), 
             color=colors_black[4], alpha=0.4, size=0.5)+
  geom_point(data=all_v, aes(x=TEMP, y=PR, color=in_out), size=0.5)+
  scale_color_manual(breaks = c(T, F), 
                    values=color_two)+
  geom_hline(yintercept=range_PR, linetype="dashed", color = colors_black[7], size=1.5)+
  geom_vline(xintercept=range_TEMP, linetype="dashed", color = colors_black[7], size=1.5)+
  xlim(c(30, 32.2))+
  ylim(c(0, 5000))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20))+
  xlab("Annual Maximum Temperature")+ylab("Annual Precipitation")
p
ggsave(p, filename="../../Figures/Methods/niche.png", width=10, height=8)
ggsave(p, filename="../../Figures/Methods/niche.pdf", width=10, height=8)
  

