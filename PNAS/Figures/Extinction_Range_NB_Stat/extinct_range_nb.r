library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(data.table)
source("commonFuns/functions.r")
source("commonFuns/colors.r")
df<-NULL
g<-"Birds"
for (g in c("Birds", "Mammals")){
  print(g)
  nb<-readRDS(sprintf("../../Objects/Species_property/%s_property.rda", g))
  
  nb<-nb%>%dplyr::select(nb_bio1_sd, nb_bio5_sd, nb_bio6_sd, nb_bio12_sd, nb_bio13_sd, nb_bio14_sd, N_CELL, sp)
  nb<-as.data.frame(nb)
  
  for (exposure in c(0, 5)){
    extinct<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
    extinct<-extinct%>%dplyr::distinct(sp, GCM, SSP, extinct_year, dispersal)
    extinct[is.infinite(extinct$extinct_year), "extinct_year"]<-2020
    df1<-full_join(nb, extinct, by="sp")
    df1<-df1%>%filter(!is.na(extinct_year))
    df1$exposure<-exposure
    df1$group<-g
    df<-bind(df, df1)
  }
}
plot(df1$extinct_year, df1$N_CELL)
source("commonFuns/colors.r")

lm_eqn <- function(df, lm_object) {
  v<-coef(lm_object)
  names(v)<-NULL
  eq <-
    substitute(
      italic(y) == a + b %.% italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2 ~ "," ~  ~ italic(p-value) ~ "=" ~ pp,
      list(
        a = format(v[1], digits = 2),
        b = format(v[2], digits = 2),
        r2 = format(summary(lm_object)$r.squared, digits = 3),
        pp = format(lmp(lm_object), digits = 2)
      )
    )
  as.character(as.expression(eq))
}
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


group<-"Mammals"
df_sp_list<-list()
for (group in c("Mammals", "Birds")){
  if (F){
    group<-"Birds"
    df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df.rda", group))
    iucn_bird<-read.csv("../../Data/Dispersal_distance/Birds/HBW-BirdLife_List_of_Birds_v5.csv", head=T, stringsAsFactors = F)
    iucn_bird<-iucn_bird[, c("Order", "Family.name", "Scientific.name")]
    iucn_bird<-data.table(iucn_bird)
    colnames(iucn_bird)<-c("Order", "family", "SP")
    df_list<-merge(df_list, iucn_bird, by="SP", all.x=T, all.y=F)
    saveRDS(df_list, sprintf("../../Objects/IUCN_List/%s_df_with_family.rda", group))
  }
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df_with_family.rda", group))
  df_list$sp<-gsub(" ", "_", df_list$SP)
  #df_list$estimated_disp
  cols<-c("sp", "estimated_disp", "family")
  df_list<-df_list[, ..cols]
  df_sp_list[[group]]<-df_list
}
df_sp_list<-rbindlist(df_sp_list)
df<-merge(df, df_sp_list, by="sp")

df_se<-df%>%dplyr::group_by(nb_bio1_sd, nb_bio5_sd, nb_bio6_sd,
                            nb_bio12_sd, nb_bio13_sd, nb_bio14_sd,
                            N_CELL, sp, dispersal, exposure, group)%>%
  dplyr::summarise(mean_extinct_year=mean(extinct_year))


  

df_se$exposure<-ifelse(df_se$exposure==0, " no exposure", "5-year exposure")
df_se$da<-ifelse(df_se$dispersal==1, "with dispersal", "no dispersal")
labels<-c("Distribution range", "Range of temperature", "Range of precipitation")
vars<-c("N_CELL", "nb_bio1_sd", "nb_bio12_sd")
i=1

             
for (i in c(1:3)){
  var<-vars[i]
  label<-labels[i]
  pp<-ggplot(df_se, aes_string(x="mean_extinct_year", y=var))+
    geom_point(size=0.1)+
    stat_smooth(method = "lm", formula = y ~ x, size = 1)+
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
                 parse = TRUE, p.digits=3)+
    theme_bw()+
    scale_x_continuous(breaks=seq(2020, 2100, 10), labels=seq(2020, 2100, 10))+
    facet_grid(da~exposure, scale="free")+
    xlab("Extinct year")+
    ylab(label)
  
  if (i==1){
    pp<-pp+scale_y_sqrt(breaks=seq(5, 30, 5)^2)
  }
  
  if (i==1){
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, N_CELL)%>%dplyr::select(sp, group, N_CELL)
    colnames(df_hist)[3]<-"V"
  }
  if (i==2){
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, nb_bio1_sd)%>%dplyr::select(sp, group, nb_bio1_sd)
    colnames(df_hist)[3]<-"V"
  }
  if (i==3){
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, nb_bio12_sd)%>%dplyr::select(sp, group, nb_bio12_sd)
    colnames(df_hist)[3]<-"V"
  }
  
  p<-ggplot(df_hist)+
    geom_density(aes(x = V, y=..count.., color=group, fill = group), position = "identity", alpha=0.3)+
    scale_fill_manual(values=color_groups)+
    scale_color_manual(values=color_groups)+
    theme_bw()+
    labs(fill="Group", color="Group", y="Number of species", x=label)
  if (i==1){
    p<-p+scale_x_log10()
  }
  p
  ppp<-ggarrange(pp, p, widths=c(3, 2))
  ggsave(ppp, filename=sprintf("../../Figures/Extinction_Range_NB_Stat/%s.png", var),
         width=15, height=4)
  ggsave(ppp, filename=sprintf("../../Figures/Extinction_Range_NB_Stat/%s.pdf", var),
         width=15, height=4)
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
dispersals<-c(0, 1)

df<-NULL
for (g in c("Birds", "Mammals")){
  print(g)
  nb<-readRDS(sprintf("../../Objects/Species_property/%s_property.rda", g))
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s_df_with_family.rda", g))
  df_list$sp<-gsub(" ", "_", df_list$SP)
  #df_list$estimated_disp
  if (g=="Mammals"){
    colnames(df_list)[8]<-"Order"
  }
  cols<-c("sp", "estimated_disp", "family", "Order")
  df_list<-df_list[, ..cols]
  nb<-inner_join(nb, df_list, by="sp")
  nb<-nb%>%dplyr::select(nb_bio1_sd, nb_bio5_sd, nb_bio6_sd, 
                         nb_bio12_sd, nb_bio13_sd, nb_bio14_sd,
                         N_CELL, estimated_disp, sp, family, Order, 
                         range_bio1_sd_min, range_bio1_sd_max,
                         range_bio5_sd_min, range_bio5_sd_max,
                         range_bio6_sd_min, range_bio6_sd_max,
                         range_bio12_sd_min, range_bio12_sd_max,
                         range_bio13_sd_min, range_bio13_sd_max,
                         range_bio14_sd_min, range_bio14_sd_max)
  
  nb<-as.data.frame(nb)
  full_sp<-expand.grid(GCM=GCMs, SSP=SSPs, sp=unique(nb$sp), dispersal=dispersals, stringsAsFactors = F)
  full_sp_with_nb<-full_join(full_sp, nb, by="sp")
  for (exposure in c(0, 5)){
    extinct<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/%s.rda", exposure, g))
    extinct<-extinct%>%dplyr::distinct(sp, GCM, SSP, extinct_year, dispersal)
    extinct[is.infinite(extinct$extinct_year), "extinct_year"]<-2020
    df1<-full_join(full_sp_with_nb, extinct, by=c("sp", "GCM", "SSP", "dispersal"))
    df1$exposure<-exposure
    df1$group<-g
    df<-bind(df, df1)
  }
}





df$is_extinct<-ifelse(is.na(df$extinct_year), "NO", "YES")
df$exposure<-ifelse(df$exposure==0, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==1, "with dispersal", "no dispersal")
df$Label<-paste(df$SSP, df$da)
p1<-ggplot(df)+
  geom_histogram(aes(x=nb_bio1_sd/100), fill=colors_black[4], bins=50)+
  geom_histogram(data=df%>%dplyr::filter(is_extinct=="YES"), 
                 aes(x=nb_bio1_sd/100), fill=colors_red[9], bins=50)+
  theme_bw()+
  xlab("Niche breadth in annual mean temperature")+
  ylab("Number of species")+
  facet_grid(exposure~Label, scale="free")
p1
saveRDS(p1, "../../Figures/NB_hist_combined/nb_temp.rda")
ggsave(p1, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_temp.pdf"), width=12, height=6)
ggsave(p1, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_temp.png"), width=12, height=6)


p2<-ggplot(df)+
  geom_histogram(aes(x=nb_bio12_sd/100), fill=colors_black[4], bins=50)+
  geom_histogram(data=df%>%dplyr::filter(is_extinct=="YES"), 
                 aes(x=nb_bio12_sd/100), fill=colors_red[9], bins=50)+
  theme_bw()+
  xlab("Niche breadth in annual precipitation")+
  ylab("Number of species")+
  facet_grid(exposure~Label)
p2
saveRDS(p2, "../../Figures/NB_hist_combined/nb_prec.rda")
#hist((df%>%dplyr::filter(is_extinct=="YES"))$estimated_disp)
p3<-ggplot(df)+
  geom_histogram(aes(x=estimated_disp), fill=colors_black[4], bins=50)+
  geom_histogram(data=df%>%dplyr::filter(is_extinct=="YES"), 
                 aes(x=estimated_disp), fill=colors_red[9], bins=50)+
  scale_x_log10()+
  scale_y_sqrt(breaks=seq(0, 80, by=10)^2)+
  theme_bw()+
  xlab("Estimated natal dispersal distance (km, log transformed)")+
  ylab("Number of species (sq root transformed)")+
  facet_grid(exposure~Label)
p3
saveRDS(p3, "../../Figures/NB_hist_combined/disp_dist.rda")


p2<-p2+
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank())
ggsave(p2, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_prec.pdf"), width=12, height=6)
ggsave(p2, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_prec.png"), width=12, height=6)


pp<-ggarrange(p1, p2, ncol=1, nrow=2)
ggsave(pp, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_combined.pdf"), width=12, height=6)
ggsave(pp, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_combined.png"), width=12, height=6)


df$range_bio1_sd_min<-df$range_bio1_sd_min/100
df$range_bio1_sd_max<-df$range_bio1_sd_max/100
df$range_bio5_sd_min<-df$range_bio5_sd_min/100
df$range_bio5_sd_max<-df$range_bio5_sd_max/100
df$range_bio6_sd_min<-df$range_bio6_sd_min/100
df$range_bio6_sd_max<-df$range_bio6_sd_max/100
df$range_bio12_sd_min<-df$range_bio12_sd_min/100
df$range_bio12_sd_max<-df$range_bio12_sd_max/100
df$range_bio13_sd_min<-df$range_bio13_sd_min/100
df$range_bio13_sd_max<-df$range_bio13_sd_max/100
df$range_bio14_sd_min<-df$range_bio14_sd_min/100
df$range_bio14_sd_max<-df$range_bio14_sd_max/100

df$scaled_range_bio1_sd_min<-scale(c(df$range_bio1_sd_min, df$range_bio1_sd_max))[1:nrow(df)]
df$scaled_range_bio1_sd_max<-scale(c(df$range_bio1_sd_max, df$range_bio1_sd_min))[1:nrow(df)]
df$scaled_range_bio5_sd_min<-scale(c(df$range_bio5_sd_min, df$range_bio5_sd_max))[1:nrow(df)]
df$scaled_range_bio5_sd_max<-scale(c(df$range_bio5_sd_max, df$range_bio5_sd_min))[1:nrow(df)]
df$scaled_range_bio6_sd_min<-scale(c(df$range_bio6_sd_min, df$range_bio6_sd_max))[1:nrow(df)]
df$scaled_range_bio6_sd_max<-scale(c(df$range_bio6_sd_max, df$range_bio6_sd_min))[1:nrow(df)]
df$scaled_range_bio12_sd_min<-scale(c(df$range_bio12_sd_min, df$range_bio12_sd_max))[1:nrow(df)]
df$scaled_range_bio12_sd_max<-scale(c(df$range_bio12_sd_max, df$range_bio12_sd_min))[1:nrow(df)]
df$scaled_range_bio13_sd_min<-scale(c(df$range_bio13_sd_min, df$range_bio13_sd_max))[1:nrow(df)]
df$scaled_range_bio13_sd_max<-scale(c(df$range_bio13_sd_max, df$range_bio13_sd_min))[1:nrow(df)]
df$scaled_range_bio14_sd_min<-scale(c(df$range_bio14_sd_min, df$range_bio14_sd_max))[1:nrow(df)]
df$scaled_range_bio14_sd_max<-scale(c(df$range_bio14_sd_max, df$range_bio14_sd_min))[1:nrow(df)]

df$scaled_nb_bio1_sd<-df$scaled_range_bio1_sd_max - df$scaled_range_bio1_sd_min
df$scaled_nb_bio5_sd<-df$scaled_range_bio5_sd_max - df$scaled_range_bio5_sd_min
df$scaled_nb_bio6_sd<-df$scaled_range_bio6_sd_max - df$scaled_range_bio6_sd_min
df$scaled_nb_bio12_sd<-df$scaled_range_bio12_sd_max - df$scaled_range_bio12_sd_min
df$scaled_nb_bio13_sd<-df$scaled_range_bio13_sd_max - df$scaled_range_bio13_sd_min
df$scaled_nb_bio14_sd<-df$scaled_range_bio14_sd_max - df$scaled_range_bio14_sd_min


df$nb_bio1_sd<-df$nb_bio1_sd/100
df$nb_bio5_sd<-df$nb_bio5_sd/100
df$nb_bio6_sd<-df$nb_bio6_sd/100
df$nb_bio12_sd<-df$nb_bio12_sd/100
df$nb_bio13_sd<-df$nb_bio13_sd/100
df$nb_bio14_sd<-df$nb_bio14_sd/100

df$N_CELL<-ceiling(df$N_CELL/100)
df$nb_volume<-df$scaled_nb_bio1_sd*df$scaled_nb_bio5_sd*df$scaled_nb_bio6_sd*
  df$scaled_nb_bio12_sd*df$scaled_nb_bio13_sd*df$scaled_nb_bio14_sd
write.csv(df, "../../Figures/raw_result.csv", row.names=F)
df$is_extinct<-as.factor(df$is_extinct)
table(df$is_extinct)
df[is.nan(is_extinct)]
item<-NULL

library(lme4)
library(RLRsim)
ggplot(df_bio12)+
  geom_histogram(aes(x=range_bio12_sd_min), fill="red")+
  geom_histogram(aes(x=range_bio12_sd_max), fill="blue")
cols<-c("sp", "range_bio12_sd_min", "range_bio12_sd_max", "group")
df_bio12<-unique(df[, ..cols])
table(df_bio12[range_bio12_sd_min<0]$group)


mean_disp_dist<-readRDS("../../Objects/Dispersal_distances/all_mean_disp_dist.rda")
mean_disp_dist$exposure<-ifelse(mean_disp_dist$exposure==0, " no exposure", "5-year exposure")
colnames(mean_disp_dist)[1]<-"sp"
mean_disp_dist$sp<-gsub(" ", "_", mean_disp_dist$sp)
mean_disp_dist[sp=="Abditomys_latidens"]
df_with_mean<-merge(df, mean_disp_dist, by=c("sp", "SSP", "GCM", "exposure", "group"))
df_with_mean<-df_with_mean[mean_dist!=estimated_disp]
df_with_mean$mean_dist<-df_with_mean$mean_dist/1000
df_with_mean$scaled_nb_volume<-scale(df_with_mean$nb_volume)
df_with_mean$scaled_mean_dist<-scale(df_with_mean$mean_dist)*-1
df_with_mean$scaled_estimated_disp<-scale(df_with_mean$estimated_disp)
df_with_mean$scaled_N_CELL<-scale(df_with_mean$N_CELL)
all_result<-NULL
for (SSPx in SSPs){
  for (GCMx in GCMs){
    for (dda in unique(df$da)){
      for (exp in unique(df$exposure)){
        for (g in unique(df$group)){
          item<-df_with_mean[(SSP==SSPx)&(GCM==GCMx)&(da==dda)&(exposure==exp)&(group==g)]
          item<-item[!is.na(family)]

          #m_glmer <- glmer(is_extinct ~ nb_volume+N_CELL+estimated_disp +
          #             (1 | family), data = item, family = binomial)
          #m_glmer <- glmer(is_extinct ~ nb_volume+N_CELL+estimated_disp+mean_dist+
          #                                (1 | family), data = item, family = binomial)
          m_glmer <- glmer(is_extinct ~ scaled_nb_volume+scaled_N_CELL+scaled_estimated_disp+scaled_mean_dist+
                             (1 | family), data = item, family = binomial)
          
          #m_glm <- glm(is_extinct ~ nb_volume+N_CELL+estimated_disp, data = item, family = binomial)
          #m_glm <- glm(is_extinct ~ nb_volume+N_CELL+estimated_disp+mean_dist, data = item, family = binomial)
          m_glm <- glm(is_extinct ~ scaled_nb_volume+scaled_N_CELL+scaled_estimated_disp+scaled_mean_dist, 
                       data = item, family = binomial)
          #logistic.display(m_glm)
          
          ano<-anova(m_glmer, m_glm)
          
          summ<-summary(m_glm)
          #m_glm<-glm(is_extinct~nb_volume+N_CELL+estimated_disp, data=item, family="binomial")
          result_item<-data.frame(SSP=SSPx, GCM=GCMx, da=dda, exposure=exp, group=g)
          
          result_item$aic<-summ$aic
          result_item$anova_pr<-ano$`Pr(>Chisq)`[2]
          result_item$anova_Chisq<-ano$Chisq[2]
          ccc3<-summ$coefficients
          result_item$Intercept_Estimate<-ccc3[1,1]
          result_item$Intercept_Std_Error<-ccc3[1,2]
          result_item$Intercept_z_value<-ccc3[1,3]
          result_item$Intercept_Pr<-ccc3[1,4]
          
          result_item$nb_volume_Estimate<-ccc3[2,1]
          result_item$nb_volume_Std_Error<-ccc3[2,2]
          result_item$nb_volume_z_value<-ccc3[2,3]
          result_item$nb_volume_Pr<-ccc3[2,4]
          
          result_item$N_CELL_Estimate<-ccc3[3,1]
          result_item$N_CELL_Std_Error<-ccc3[3,2]
          result_item$N_CELL_z_value<-ccc3[3,3]
          result_item$N_CELL_Pr<-ccc3[3,4]
          
          result_item$estimated_disp_Estimate<-ccc3[4,1]
          result_item$estimated_disp_Std_Error<-ccc3[4,2]
          result_item$estimated_disp_z_value<-ccc3[4,3]
          result_item$estimated_disp_Pr<-ccc3[4,4]
          
          result_item$mean_disp_Estimate<-ccc3[5,1]
          result_item$mean_disp_Std_Error<-ccc3[5,2]
          result_item$mean_disp_z_value<-ccc3[5,3]
          result_item$mean_disp_Pr<-ccc3[5,4]
          
          
          ccc<-summ$coefficients[,4]
          result_item$p_value_nb_volume<-ccc[2]
          result_item$p_value_nb_volume_label<-""
          if (result_item$p_value_nb_volume<0.05){
            result_item$p_value_nb_volume_label<-"p<0.05*"
          }
          if (result_item$p_value_nb_volume<0.01){
            result_item$p_value_nb_volume_label<-"p<0.01**"
          }
          if (result_item$p_value_nb_volume<0.001){
            result_item$p_value_nb_volume_label<-"p<0.001***"
          }
          result_item$p_value_N_CELL<-ccc[3]
          result_item$p_value_N_CELL_label<-""
          if (result_item$p_value_N_CELL<0.05){
            result_item$p_value_N_CELL_label<-"p<0.05*"
          }
          if (result_item$p_value_N_CELL<0.01){
            result_item$p_value_N_CELL_label<-"p<0.01**"
          }
          if (result_item$p_value_N_CELL<0.001){
            result_item$p_value_N_CELL_label<-"p<0.001***"
          }
          result_item$p_value_estimated_disp<-ccc[4]
          result_item$p_value_estimated_disp_label<-""
          if (result_item$p_value_estimated_disp<0.05){
            result_item$p_value_estimated_disp_label<-"p<0.05*"
          }
          if (result_item$p_value_estimated_disp<0.01){
            result_item$p_value_estimated_disp_label<-"p<0.01**"
          }
          if (result_item$p_value_estimated_disp<0.001){
            result_item$p_value_estimated_disp_label<-"p<0.001***"
          }
          
          result_item$p_value_mean_disp<-ccc[5]
          result_item$p_value_mean_disp_label<-""
          if (result_item$p_value_mean_disp<0.05){
            result_item$p_value_mean_disp_label<-"p<0.05*"
          }
          if (result_item$p_value_mean_disp<0.01){
            result_item$p_value_mean_disp_label<-"p<0.01**"
          }
          if (result_item$p_value_mean_disp<0.001){
            result_item$p_value_mean_disp_label<-"p<0.001***"
          }
          
          all_result<-bind_dplyr(all_result, result_item)
          #for m_glmer
          if (F){
            table(item$family)
            hist(item$nb_volume)
            hist(item$estimated_disp)
            
            summ<-summary(m_glmer)
            #m_glm<-glm(is_extinct~nb_volume+N_CELL+estimated_disp, data=item, family="binomial")
            result_item<-data.frame(SSP=SSPx, GCM=GCMx, da=dda, exposure=exp, group=g)
            
            result_item$aic<-summ$AICtab[1]
            result_item$bic<-summ$AICtab[2]
            
            ccc<-summ$coefficients[,4]
            result_item$p_value_nb_volume<-ccc[2]
            result_item$p_value_nb_volume_label<-""
            if (result_item$p_value_nb_volume<0.05){
              result_item$p_value_nb_volume_label<-"p<0.05*"
            }
            if (result_item$p_value_nb_volume<0.01){
              result_item$p_value_nb_volume_label<-"p<0.01**"
            }
            if (result_item$p_value_nb_volume<0.001){
              result_item$p_value_nb_volume_label<-"p<0.001***"
            }
            result_item$p_value_N_CELL<-ccc[3]
            result_item$p_value_N_CELL_label<-""
            if (result_item$p_value_N_CELL<0.05){
              result_item$p_value_N_CELL_label<-"p<0.05*"
            }
            if (result_item$p_value_N_CELL<0.01){
              result_item$p_value_N_CELL_label<-"p<0.01**"
            }
            if (result_item$p_value_N_CELL<0.001){
              result_item$p_value_N_CELL_label<-"p<0.001***"
            }
            result_item$p_value_estimated_disp<-ccc[4]
            result_item$p_value_estimated_disp_label<-""
            if (result_item$p_value_estimated_disp<0.05){
              result_item$p_value_estimated_disp_label<-"p<0.05*"
            }
            if (result_item$p_value_estimated_disp<0.01){
              result_item$p_value_estimated_disp_label<-"p<0.01**"
            }
            if (result_item$p_value_estimated_disp<0.001){
              result_item$p_value_estimated_disp_label<-"p<0.001***"
            }
            
            all_result<-bind_dplyr(all_result, result_item)
          }
        }
      }
    }
  }
}

#write.csv(all_result, "../../Figures/NB_hist_combined/p_values.csv", row.names=F)
#write.csv(all_result, "../../Figures/NB_hist_combined/p_values_with_mean_disp_dist.csv", row.names=F)
write.csv(all_result, "../../Figures/NB_hist_combined/p_values_with_scaled_mean_disp_dist.csv", row.names=F)

all_result_raw<-read.csv("../../Figures/NB_hist_combined/p_values_with_mean_disp_dist.csv")
all_result_raw<-data.table(all_result_raw)
all_result_raw<-all_result_raw[da=="with dispersal"]
all_result_raw<-all_result_raw[exposure==" no exposure"]


all_result<-read.csv("../../Figures/NB_hist_combined/p_values_with_scaled_mean_disp_dist.csv")
all_result<-data.table(all_result)
all_result<-all_result[da=="with dispersal"]
all_result<-all_result[exposure==" no exposure"]

vars<-c("nb_volume", "N_CELL", "estimated_disp", "mean_disp")
df_g<-list()
v<-vars[1]
for (v in vars){
  cols<-c("SSP", "GCM", "da", "exposure", "group", "aic", "anova_pr", "anova_Chisq", 
          "Intercept_Estimate", "Intercept_Std_Error", "Intercept_z_value" , "Intercept_Pr",                
          sprintf("%s_Estimate", v), sprintf("%s_Std_Error", v), sprintf("%s_z_value", v),
          sprintf("%s_Pr", v), sprintf("p_value_%s", v), sprintf("p_value_%s_label", v))
  item<-all_result[, ..cols]
  colnames(item)<-c("SSP", "GCM", "da", "exposure", "group", "aic", "anova_pr", "anova_Chisq", 
                    "Intercept_Estimate", "Intercept_Std_Error", "Intercept_z_value" , "Intercept_Pr",                
                    sprintf("%s_Estimate", "v"), sprintf("%s_Std_Error", "v"), sprintf("%s_z_value", "v"),
                    sprintf("%s_Pr", "v"), sprintf("p_value_%s", "v"), sprintf("p_value_%s_label", "v"))
  item$predictor<-v
  df_g[[v]]<-item
}
df_g<-rbindlist(df_g)

p<-ggplot(df_g)+
  geom_errorbar(data=df_g, aes(x=GCM, ymin=v_Estimate-v_Std_Error, 
                               ymax=v_Estimate+v_Std_Error, color=predictor), width=0.1)+
  geom_point(data=df_g, aes(x=GCM, y=v_Estimate, color=predictor, 
                            shape=factor(p_value_v_label)), size=3)+
  
  facet_grid(SSP~group, scale="free")+theme_bw()
p
ggsave(p, filename="../../Figures/NB_hist_combined/Estimate.pdf")
