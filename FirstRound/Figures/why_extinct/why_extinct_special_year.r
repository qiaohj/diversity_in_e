library(dplyr)
library(ggplot2)
library(raster)
source("commonFuns/colors.r")
source("commonFuns/functions.r")

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
df<-readRDS("../../Objects_Full_species/why_extinct/why_extinct.rda")
head(df)
df$exposure<-ifelse(df$threshold==1, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==0, "no dispersal", "with dispersal")

#factor by year
df_se_temp<-df%>%dplyr::filter(!is_temp_in)%>%dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_temp$causation<-"Temperature"
df_se_prec<-df%>%dplyr::filter(!is_prec_in)%>%dplyr::group_by(SSP, extinct_year, da, exposure)%>%
  dplyr::summarise(N=n())
df_se_prec$causation<-"Precipitation"

df_se<-bind_rows(df_se_prec, df_se_temp)

df_se_prec%>%ungroup()%>%dplyr::filter((N>1000))

threshold<-1
ttt<-2
g<-"Amphibians"
env_layers<-readRDS("../../Objects_Full_species/stacked_layers_2021_2100.rda")
all_result<-NULL
da<-0
target_year<-2032
target_SSP<-"SSP119"
df_all<-NULL
all_result<-NULL
for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
  sp_list<-readRDS(sprintf("../../Objects_Full_species/IUCN_List/%s.rda", g))
  sp_list<-sp_list[which(sp_list$area>ttt),]
  
  df<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
  sp_list$sp2<-gsub(" ", "_", sp_list$sp)
  df<-df%>%dplyr::filter(sp%in%sp_list$sp2)
  df<-df%>%dplyr::filter((extinct_year==target_year)&(dispersal==da)&(SSP==target_SSP))
  df_all<-bind_dplyr(df_all, df)
  when_extinct<-df%>%dplyr::distinct(group, sp, GCM, SSP, extinct_year, dispersal)
  for (i in c(1:nrow(when_extinct))){
    print(paste(threshold, g, i, nrow(when_extinct)))
    item<-when_extinct[i,]
    
    target_folder<-sprintf("../../Objects_Full_species/Niche_Models/%s/%s", g, item$sp)
    fit<-readRDS(sprintf("%s/fit.rda", target_folder))
    dis<-readRDS(sprintf("%s/dispersal_%d/%s_%s_%d.rda",
                         target_folder, threshold, item$GCM, item$SSP, item$dispersal))
    dis<-dis%>%dplyr::filter(YEAR==(item$extinct_year-1))
    env_layer<-env_layers[[sprintf("%s_%s_%d", item$GCM, item$SSP, item$extinct_year)]]
    env_layer<-env_layer%>%dplyr::filter(mask_index %in% dis$mask_index)
    env_layer$is_temp_in<-between(env_layer$TEMP_MAX, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)&
      between(env_layer$TEMP_MIN, fit$range_TEMP_sd_min, fit$range_TEMP_sd_max)
    env_layer$is_prec_in<-between(env_layer$PR, fit$range_PR_sd_min, fit$range_PR_sd_max)
    for (nn in names(item)){
      env_layer[, nn]<-item[, nn]
    }
    env_layer$threshold<-threshold
    all_result<-bind_dplyr(all_result, env_layer)
  }
}
df_pos<-df_all%>%dplyr::group_by(x, y, SSP)%>%
  dplyr::summarise(N=n())

range(df_pos$N)

all_result_se<-all_result%>%dplyr::filter(!is_prec_in)%>%dplyr::group_by(x, y, GCM)%>%
  dplyr::summarise(N=n())


ggplot(all_result_se)+geom_tile(aes(x=x, y=y, fill=N))+
  facet_wrap(~GCM)

target_years<-c(2025:2040)
target_env<-"UKESM1_SSP119_%d"
r_continent<-raster("../../Raster/Continent_ect4.tif")
yy=2032
env_se_all<-NULL
for (yy in target_years){
  env<-env_layers[[sprintf(target_env, yy)]]
  env$continent<-extract(r_continent, env[, c("x", "y")])
  env<-env%>%filter(!is.na(continent))
  env_se<-env%>%dplyr::group_by(continent)%>%
    dplyr::summarise(mean_PR=mean(PR),
                     sd_PR=sd(PR),
                     mean_TEMP_MAX=mean(TEMP_MAX),
                     mean_TEMP_MIN=mean(TEMP_MIN),
                     year=yy)
  env_se_all<-bind_dplyr(env_se_all, env_se)
  #plot(env$x, env$y)
}

ggplot(env_se_all) + 
  geom_ribbon(aes(x=year, ymin=mean_PR-sd_PR, ymax=mean_PR+sd_PR, fill=factor(continent)), alpha=0.2)+
  geom_line(aes(x=year, y=mean_PR, color=factor(continent)))
