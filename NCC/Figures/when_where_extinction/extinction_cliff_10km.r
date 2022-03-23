library(dplyr)
library(raster)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
source("commonFuns/functions.r")
source("commonFuns/colors.r")
exposure<-0

if (F){
  points_10km<-readRDS("../../Raster/points_10km.rda")
  colnames(points_10km)[c(1,2)]<-c("x_10km", "y_10km")
  mask_100km<-raster("../../Raster/mask_100km.tif")
  points_100km<-data.table(rasterToPoints(mask_100km))
  colnames(points_100km)[c(1,2)]<-c("x_100km", "y_100km")
  points_10km<-merge(points_10km, points_100km, by="mask_100km")
  
  result<-list()
  for (exposure in c(0,5)){
    g<-"Mammals"
    for (g in c("Birds", "Mammals")){
      df<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s_10km.rda", exposure, g))
      df<-rbindlist(df)
      df<-df%>%dplyr::filter(!is.infinite(extinct_year))
      
      sp_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", g))
      sp_list$sp<-sp_list$SP
      sp_list$sp2<-gsub(" ", "_", sp_list$sp)
      df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
      extinct_sp_list<-df%>%dplyr::distinct(group, sp, GCM, SSP, dispersal, extinct_year)
      i=10
      for (i in c(1:nrow(extinct_sp_list))){
        print(paste(exposure, g, i, nrow(extinct_sp_list)))
        item<-extinct_sp_list[i,]
        df_item<-df%>%dplyr::filter((sp==item$sp)&(GCM==item$GCM)&
                                      (SSP==item$SSP)&(dispersal==item$dispersal))
        N_st_Cell<-nrow(df_item)
        if (item$extinct_year==2021){
          N_extinct_Cell<-N_st_Cell
          Max_N_Cell<-N_st_Cell
        }else{
          f<-sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d_10km.rda", 
                     item$group, item$sp, item$GCM, item$SSP, exposure, item$dispersal)
          res<-"10km"
          if (!file.exists(f)){
            f<-sprintf("../../Objects/Dispersal/%s/%s/%s_%s_%d_dispersal_%d.rda", 
                       item$group, item$sp, item$GCM, item$SSP, exposure, item$dispersal)
            res<-"100km"
          }
          if (!file.exists(f)){
            next()
          }
          
          future_dis<-readRDS(f)
          future_dis<-rbindlist(future_dis)
          if (res=="10km"){
            
            future_dis<-merge(future_dis, points_10km, by.x="mask_10km", by.y="mask_10km")
            future_dis<-future_dis[, .(N=.N), by=c("x_100km", "y_100km", "mask_100km", "YEAR")]
            colnames(future_dis)[1:3]<-c("x", "y", "mask")
          }
          
          future_dis_sum<-future_dis%>%dplyr::group_by(YEAR)%>%
            dplyr::summarise(N_CELL=n())
          Max_N_Cell<-max(future_dis_sum$N_CELL)
          
          N_extinct_Cell<-pull(future_dis_sum[which(future_dis_sum$YEAR==(item$extinct_year-1)), "N_CELL"])
        }
        item$N_st_Cell<-N_st_Cell
        item$N_extinct_Cell<-N_extinct_Cell
        item$Max_N_Cell<-Max_N_Cell
        item$res<-res
        item$exposure<-ifelse(exposure==0, " no climate resilience", "climate resilience")
        result[[length(result)+1]]<-item
      }
    }
  }
  saveRDS(result, "../../Figures/Extinction_cliff/Extinction_cliff_10km.rda")
}

cliff<-readRDS("../../Figures/Extinction_cliff/Extinction_cliff_10km.rda")
cliff<-rbindlist(cliff)
cliff[exposure==" no exposure"]$exposure<-" no climate resilience"
cliff[exposure=="5-year exposure"]$exposure<-"climate resilience"
cliff<-cliff[!is.na(N_extinct_Cell)]
range(cliff$N_extinct_Cell)
hist(cliff$N_extinct_Cell)




cliff_se<-cliff%>%dplyr::group_by(SSP, dispersal, exposure, group)%>%
  dplyr::summarise(mean_extinct_Cell=mean(N_extinct_Cell),
                   sd_extinct_Cell=sd(N_extinct_Cell),
                   CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])

cliff_se_all<-cliff%>%dplyr::group_by(SSP, dispersal, exposure)%>%
  dplyr::summarise(group="ALL",
                   mean_extinct_Cell=mean(N_extinct_Cell),
                   sd_extinct_Cell=sd(N_extinct_Cell),
                   CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])
cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")
cliff$da=ifelse(cliff$dispersal==0, "no dispersal", "with dispersal")
cliff[which(cliff$exposure=="no exposure"), "exposure"]<-" no climate resilience"
cliff[which(cliff$exposure==" 5-year exposure"), "exposure"]<-"climate resilience"

p<-ggplot(cliff)+geom_density(aes(x=N_extinct_Cell, color=da))+
  facet_grid(exposure~SSP)+
  xlab("Number of cells")+
  ylab("Density")+
  labs(color="")+
  scale_x_log10()+
  scale_color_manual(values=color_da)+
  theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_10km.pdf"), 
       width=12, height=6)
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_10km.png"), 
       width=12, height=6)


cliff_se_all<-cliff%>%dplyr::group_by(SSP, dispersal, exposure, extinct_year)%>%
  dplyr::summarise(group="ALL",
                   mean_extinct_Cell=mean(N_extinct_Cell),
                   sd_extinct_Cell=sd(N_extinct_Cell),
                   CI_extinct_Cell=CI(N_extinct_Cell)[1]-CI(N_extinct_Cell)[2])
cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")

p<-ggplot(cliff_se_all)+
  geom_ribbon(aes(x=extinct_year,
                  ymin=mean_extinct_Cell-CI_extinct_Cell, 
                  ymax=mean_extinct_Cell+CI_extinct_Cell, 
                  fill=da), alpha=0.2)+
  geom_line(aes(x=extinct_year, y=mean_extinct_Cell, color=da))+
  facet_grid(exposure~SSP, scale="free")+
  xlab("Year")+
  ylab("Number of cells")+
  labs(color="", fill="")+
  scale_color_manual(values=color_da)+
  scale_fill_manual(values=color_da)+
  theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_10km.pdf"), 
       width=12, height=6)
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_10km.png"), 
       width=12, height=6)

cliff$extinct_proportion_st<-cliff$N_extinct_Cell/cliff$N_st_Cell
cliff$Max_N_Cell<-ifelse(cliff$N_extinct_Cell>cliff$Max_N_Cell, cliff$N_extinct_Cell, cliff$Max_N_Cell)
cliff$extinct_proportion_max<-cliff$N_extinct_Cell/cliff$Max_N_Cell

cliff_se_all<-cliff%>%dplyr::group_by(SSP, dispersal, exposure, extinct_year)%>%
  dplyr::summarise(group="ALL",
                   N=n(),
                   mean_extinct_st=mean(extinct_proportion_st),
                   sd_extinct_st=sd(extinct_proportion_st),
                   CI_extinct_st=CI(extinct_proportion_st)[1]-CI(extinct_proportion_st)[2],
                   mean_extinct_max=mean(extinct_proportion_max),
                   sd_extinct_max=sd(extinct_proportion_max),
                   CI_extinct_max=CI(extinct_proportion_max)[1]-CI(extinct_proportion_max)[2])
cliff_se_all$da=ifelse(cliff_se_all$dispersal==0, "no dispersal", "with dispersal")

cliff_se_all%>%dplyr::filter(N==1)
cliff%>%dplyr::filter((SSP=="SSP119")&(da=="no dispersal")&(exposure==" no exposure")&(extinct_year==2098))

p<-ggplot(cliff_se_all)+
  geom_ribbon(aes(x=extinct_year,
                  ymin=mean_extinct_st-sd_extinct_st, 
                  ymax=mean_extinct_st+sd_extinct_st, 
                  fill=da), alpha=0.2)+
  geom_line(aes(x=extinct_year, y=mean_extinct_st, color=da))+
  facet_grid(exposure~SSP, scale="free")+
  xlab("Year")+
  ylab("Number of cells/Initial number of cells")+
  labs(color="", fill="")+
  scale_color_manual(values=color_da)+
  scale_fill_manual(values=color_da)+
  theme_bw()
p
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_inital_cell_10km.pdf"), 
       width=12, height=6)
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_inital_cell_10km.png"), 
       width=12, height=6)

p<-ggplot(cliff_se_all)+
  geom_ribbon(aes(x=extinct_year,
                  ymin=mean_extinct_max-sd_extinct_max, 
                  ymax=mean_extinct_max+sd_extinct_max, 
                  fill=da), alpha=0.2)+
  geom_line(aes(x=extinct_year, y=mean_extinct_max, color=da))+
  facet_grid(exposure~SSP, scale="free")+
  xlab("Year")+
  ylab("Number of cells/Max number of cells")+
  labs(color="", fill="")+
  scale_color_manual(values=color_da)+
  scale_fill_manual(values=color_da)+
  scale_y_continuous(breaks=seq(0, 1, 0.2), labels = seq(0, 1, 0.2))+
  theme_bw()
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_max_cell_10km.pdf"), 
       width=12, height=6)
ggsave(p, filename=sprintf("../../Figures/Extinction_cliff/Extinction_cliff_by_year_max_cell_10km.png"), 
       width=12, height=6)
p
 