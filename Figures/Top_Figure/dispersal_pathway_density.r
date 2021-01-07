

library(rgl)
library(ceramic)
library(anglr)
library(ggnewscale)
library(ggplot2)
library(raster)
library(data.table)
library(dplyr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#alt<-raster("../../Raster/ALT/alt_eck4.tif")
mask<-raster("../../Raster/mask_index.tif")
source("commonFuns/colors.r")
GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")


predict_range<-c(2021:2100)
layer_df<-expand.grid(GCM=GCMs, SSP=SSPs)
layer_df$LABEL<-paste(layer_df$GCM, layer_df$SSP, sep="_")
j=1
if (T){
  #add mask_index
  for (j in c(1:nrow(layer_df))){
    for (threshold in c(1, 5)){
      layer_item<-layer_df[j,]
      target<-sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s_with_index.rda", 
                      threshold, layer_item$LABEL)
      print(target)
      if (file.exists(target)){
        print("SKIP")
        next()
      }
      saveRDS(NULL, target)
      smooth_path<-readRDS(sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s.rda", 
                                   threshold, layer_item$LABEL))
      #smooth_path$YEAR<-round(smooth_path$YEAR)
      smooth_path$index<-raster::extract(mask, data.frame(x=smooth_path$x, y=smooth_path$y))
      vars<-c("group", "sp", "survive", "continent_i", "index")
      smooth_path<-smooth_path[, ..vars]
      smooth_path<-unique(smooth_path)
      saveRDS(smooth_path, target)
    }
  }
  
  df_all<-NULL
  for (j in c(1:nrow(layer_df))){
    for (threshold in c(1, 5)){
      layer_item<-layer_df[j,]
      target<-sprintf("../../Figures_Full_species/Top_Figure_%d/smooth_path_%s_with_index.rda", 
                      threshold, layer_item$LABEL)
      print(target)
      df<-readRDS(target)
      df$SSP<-layer_item$SSP
      df$GCM<-layer_item$GCM
      df$threshold<-threshold
      df_all<-bind(df_all, df)
    }
  }
  saveRDS(df_all, "../../Figures_Full_species/Top_Figure_all/Data/smooth_path_all.rda")
}

df_all<-readRDS("../../Figures_Full_species/Top_Figure_all/Data/smooth_path_all.rda")

df_all<-as_tibble(df_all)
df_se_survive<-df_all%>%dplyr::group_by(survive, index, SSP, threshold, GCM)%>%
  dplyr::summarise(N=n())
df_se_survive<-df_se_survive%>%dplyr::group_by(survive, index, SSP, threshold)%>%
  dplyr::summarise(mean_N=mean(N, na.rm=T))


df_se_all<-df_all%>%dplyr::group_by(index, SSP, threshold, GCM)%>%
  dplyr::summarise(N=n())
df_se_all<-df_se_all%>%dplyr::group_by(index, SSP, threshold)%>%
  dplyr::summarise(mean_N=mean(N, na.rm=T))
df_se_all<-df_se_all[!is.na(df_se_all$index),]

df_se_group<-df_all%>%dplyr::group_by(group, index, SSP, threshold, GCM)%>%
  dplyr::summarise(N=n())
df_se_group<-df_se_group%>%dplyr::group_by(group, index, SSP, threshold)%>%
  dplyr::summarise(mean_N=mean(N, na.rm=T))
df_se_group<-df_se_group[!is.na(df_se_group$index),]

df_se_group_survive<-df_all%>%dplyr::group_by(survive, group, index, SSP, threshold, GCM)%>%
  dplyr::summarise(N=n())
df_se_group_survive<-df_se_group_survive%>%dplyr::group_by(survive, group, index, SSP, threshold)%>%
  dplyr::summarise(mean_N=mean(N, na.rm=T))
df_se_group_survive<-df_se_group_survive[!is.na(df_se_group_survive$index),]


mask_p<-data.frame(rasterToPoints(mask))
df_se_survive<-df_se_survive[!is.na(df_se_survive$index),]
colnames(mask_p)[3]<-"index"
yyy<-2021
SSP_i<-SSPs[1]
no_na<-!is.na(values(mask))
threshold_i<-1
for (SSP_i in SSPs){
  for (threshold_i in c(1,5)){
    print(paste(SSP_i, threshold_i))
    df_se_all_item<-df_se_all%>%dplyr::filter((SSP==SSP_i)&(threshold==threshold_i))
    df_se_all_item_p<-left_join(mask_p, df_se_all_item, by="index")
    r_mask<-mask
    values(r_mask)[no_na]<-df_se_all_item_p$mean_N
    writeRaster(r_mask, 
                sprintf("../../Figures_Full_species/Top_Figure_all/TIF/%s_%d.tif", 
                        SSP_i, threshold), overwrite=T)
    for (su in unique(df_se_survive$survive)){
      df_se_survive_item<-df_se_survive%>%
        dplyr::filter((SSP==SSP_i)&(threshold==threshold_i)&(survive==su))
      df_se_survive_item_p<-left_join(mask_p, df_se_survive_item, by="index")
      r_mask<-mask
      values(r_mask)[no_na]<-df_se_survive_item_p$mean_N
      writeRaster(r_mask, 
                  sprintf("../../Figures_Full_species/Top_Figure_all/TIF/%s_%d_%s.tif", 
                          SSP_i, threshold, su), overwrite=T)
      for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
        df_se_group_survive_item<-df_se_group_survive%>%
          dplyr::filter((SSP==SSP_i)&(threshold==threshold_i)&(survive==su)&(group==g))
        df_se_group_survive_item_p<-left_join(mask_p, df_se_group_survive_item, by="index")
        r_mask<-mask
        values(r_mask)[no_na]<-df_se_group_survive_item_p$mean_N
        writeRaster(r_mask, 
                    sprintf("../../Figures_Full_species/Top_Figure_all/TIF/%s_%d_%s_%s.tif", 
                            SSP_i, threshold, su, g), overwrite=T)
      }
    }
    for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
      df_se_group_item<-df_se_group%>%dplyr::filter((SSP==SSP_i)&(threshold==threshold_i)&(group==g))
      df_se_group_item_p<-left_join(mask_p, df_se_group_item, by="index")
      r_mask<-mask
      values(r_mask)[no_na]<-df_se_group_item_p$mean_N
      writeRaster(r_mask, 
                  sprintf("../../Figures_Full_species/Top_Figure_all/TIF/%s_%d_%s.tif", 
                          SSP_i, threshold, g), overwrite=T)
    }
  }
}
tail(table(df_se_all$mean_N))
df_se_all$exposure<-ifelse(df_se_all$threshold==1, " no exposure", "5-year exposure")
hist(df_se_all$mean_N)
df_se_all_p<-inner_join(df_se_all, mask_p,  by="index")
my_breaks = c(1, 10, 100, 900)
my_breaks_label<-c("1", "10", "100", round(max(df_se_all$mean_N)))
p<-ggplot()+
  geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  geom_tile(data=df_se_all_p, aes(x=x, y=y, fill=mean_N))+
  facet_grid(exposure~SSP)+
  scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                      trans = "log", breaks=my_breaks, labels=my_breaks_label)+
  labs(fill = "Number of pathways")+
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
p

ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways.png", 
       width=9, height=4)
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways.pdf", 
       width=9, height=4)

df_se_survive$exposure<-ifelse(df_se_survive$threshold==1, " no exposure", "5-year exposure")
hist(df_se_survive$mean_N)
df_se_survive_p<-inner_join(df_se_survive, mask_p,  by="index")

xxx<-df_se_survive_p%>%dplyr::filter(survive=="EXTINCT")
range(xxx$mean_N)
my_breaks = c(1, 10, 100, 270)
my_breaks_label<-c("1", "10", "100", round(max(xxx$mean_N)))

p<-ggplot()+
  geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  geom_tile(data=xxx, aes(x=x, y=y, fill=mean_N))+
  facet_grid(exposure~SSP)+
  scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                      trans = "log", breaks=my_breaks, labels=my_breaks_label)+
  labs(fill = "Number of pathways")+
  ggtitle("Pathways of extinct species")+
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
p

ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways_Extinct.png", 
       width=9, height=4)
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways_Extinct.pdf", 
       width=9, height=4)

xxx<-df_se_survive_p%>%dplyr::filter(survive=="SURVIVE")
range(xxx$mean_N)
my_breaks = c(1, 10, 100, 900)
my_breaks_label<-c("1", "10", "100", round(max(xxx$mean_N)))

p<-ggplot()+
  geom_tile(data=mask_p, aes(x=x, y=y), fill=mask_color)+
  geom_tile(data=xxx, aes(x=x, y=y, fill=mean_N))+
  facet_grid(exposure~SSP)+
  scale_fill_gradient(low=color_two_map[1], high=color_two_map[2],
                      trans = "log", breaks=my_breaks, labels=my_breaks_label)+
  labs(fill = "Number of pathways")+
  ggtitle("Pathways of survive species")+
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
p

ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways_Survive.png", 
       width=9, height=4)
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/N_Pathways_Survive.pdf", 
       width=9, height=4)

library(Rmisc)

HFP<-raster("../../Raster/HFP2009_Low.tif")
HFP_p<-data.frame(rasterToPoints(HFP))
HFP_p$round_HFP<-round(HFP_p$HFP2009_Low/5)*5
p<-ggplot(HFP_p)+geom_histogram(aes(x=HFP2009_Low))+
  labs(x="Humen foot print", y="Count")+
  theme_bw()
ggsave(p, filename = "../../Figures_Full_species/Top_Figure_all/HFP_Hist.png", width=6, height=5)
ggsave(p, filename = "../../Figures_Full_species/Top_Figure_all/HFP_Hist.pdf", width=6, height=5)

HFP_t<-data.frame(table(HFP_p$round_HFP), stringsAsFactors = F)

colnames(HFP_t)<-c("round_HFP", "N_round_HFP")
HFP_t$round_HFP<-as.numeric(as.character(HFP_t$round_HFP))

df_se_survive_p$HFP<-raster::extract(HFP, df_se_survive_p[, c("x", "y")])
df_se_survive_p$round_HFP<-round(df_se_survive_p$HFP/5)*5
df_se_survive_p<-df_se_survive_p[which(!is.na(df_se_survive_p$round_HFP)),]
df_se_survive_p$survive_factor<-factor(df_se_survive_p$survive, levels=c("EXTINCT", "SURVIVE"))

df_se_survive_p_se<-df_se_survive_p%>%dplyr::group_by(exposure, SSP, round_HFP, survive_factor, survive)%>%
  dplyr::summarise(N=sum(mean_N),
                   sd_N=sd(mean_N),
                   CI_N=CI(mean_N)[1]-CI(mean_N)[2])
df_se_survive_p_se<-inner_join(df_se_survive_p_se, HFP_t, by="round_HFP")
df_se_survive_p_se[is.na(df_se_survive_p_se)]<-0
df_se_survive_p_se$usage_proportion<-df_se_survive_p_se$N/df_se_survive_p_se$N_round_HFP
p<-ggplot(df_se_survive_p_se)+geom_line(aes(x=round_HFP, y=usage_proportion, color=survive_factor))+
  #geom_ribbon(aes(x=round_HFP, y=N, ymin=N-CI_N,ymax=N+CI_N, group=survive_factor), fill=colors_black[2], alpha=0.5)+
  scale_y_log10()+
  facet_grid(exposure~SSP)+
  labs(x="Human foot print", y="Mean number of pathways per pixel", color="")+
  theme_bw()
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_proportion.pdf", 
       width=16, height=6)
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_proportion.png", 
       width=16, height=6)

p<-ggplot(df_se_survive_p_se)+geom_line(aes(x=round_HFP, y=N, color=survive_factor))+
  #geom_ribbon(aes(x=round_HFP, y=N, ymin=N-CI_N,ymax=N+CI_N, group=survive_factor), fill=colors_black[2], alpha=0.5)+
  scale_y_log10()+
  facet_grid(exposure~SSP)+
  labs(x="Human foot print", y="Number of pathways", color="")+
  theme_bw()
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve.pdf", 
       width=16, height=6)
ggsave(p, filename="../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve.png", 
       width=16, height=6)

for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  df_se_group_survive_item<-df_se_group_survive%>%dplyr::filter(group==g)
  df_se_group_survive_item<-inner_join(df_se_group_survive_item, mask_p, by="index")
  df_se_group_survive_item$HFP<-raster::extract(HFP, df_se_group_survive_item[, c("x", "y")])
  df_se_group_survive_item$round_HFP<-round(df_se_group_survive_item$HFP/5)*5
  df_se_group_survive_item<-df_se_group_survive_item[which(!is.na(df_se_group_survive_item$round_HFP)),]
  df_se_group_survive_item$survive_factor<-factor(df_se_group_survive_item$survive, levels=c("EXTINCT", "SURVIVE"))
  df_se_group_survive_item$exposure<-ifelse(df_se_group_survive_item$threshold==1, " no exposure", "5-year exposure")
  df_se_group_survive_item_se<-df_se_group_survive_item%>%
    dplyr::group_by(exposure, SSP, round_HFP, survive_factor, survive)%>%
    dplyr::summarise(N=sum(mean_N),
                     sd_N=sd(mean_N),
                     CI_N=CI(mean_N)[1]-CI(mean_N)[2])
  df_se_group_survive_item_se<-inner_join(df_se_group_survive_item_se, HFP_t, by="round_HFP")
  df_se_group_survive_item_se[is.na(df_se_group_survive_item_se)]<-0
  df_se_group_survive_item_se$usage_proportion<-df_se_group_survive_item_se$N/df_se_group_survive_item_se$N_round_HFP
  
  p<-ggplot(df_se_group_survive_item_se)+geom_line(aes(x=round_HFP, y=usage_proportion, color=survive_factor))+
    #geom_ribbon(aes(x=round_HFP, y=N, ymin=N-CI_N,ymax=N+CI_N, group=survive_factor), fill=colors_black[2], alpha=0.5)+
    scale_y_log10()+
    ggtitle(g)+
    facet_grid(exposure~SSP)+
    labs(x="Human foot print", y="Mean number of pathways per pixel", color="")+
    theme_bw()
  ggsave(p, filename=sprintf(
    "../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_proportion_%s.pdf",g), 
         width=16, height=6)
  ggsave(p, filename=sprintf(
    "../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_proportion_%s.png",g), 
    width=16, height=6)
  
  p<-ggplot(df_se_group_survive_item_se)+geom_line(aes(x=round_HFP, y=N, color=survive_factor))+
    #geom_ribbon(aes(x=round_HFP, y=N, ymin=N-CI_N,ymax=N+CI_N, group=survive_factor), fill=colors_black[2], alpha=0.5)+
    scale_y_log10()+
    facet_grid(exposure~SSP)+
    ggtitle(g)+
    labs(x="Human foot print", y="Number of pathways", color="")+
    theme_bw()
  ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_%s.pdf", g),
         width=16, height=6)
  ggsave(p, filename=sprintf("../../Figures_Full_species/Top_Figure_all/Pathways_HFP_Curve_%s.png", g),
         width=16, height=6)
}
