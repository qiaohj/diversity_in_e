library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(data.table)
df<-NULL
g<-"Birds"
for (g in c("Birds", "Mammals")){
  print(g)
  nb<-readRDS(sprintf("../../Objects/Species_property/%s_property_compared.rda", g))
  
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

df_se<-df%>%dplyr::group_by(nb_bio1_sd, nb_bio5_sd, nb_bio6_sd,
                            nb_bio12_sd, nb_bio13_sd, nb_bio14_sd,
                            N_CELL, sp, dispersal, exposure, group)%>%
  dplyr::summarise(mean_extinct_year=mean(extinct_year))


fit<-lm(extinct_year~N_CELL, data=df1%>%filter(dispersal==0)) 
eqn <- lm_eqn(item_g, fit)
result<-summary(fit)
r<-sqrt(result$r.squared)
f <- summary(fit)$fstatistic
p <- pf(f[1],f[2],f[3],lower.tail=F)
b<-result$coefficients[1,1]
a<-result$coefficients[2,1]


#p<-ggscatter(item_g, x = "temperature", y = "mean_v", add = "reg.line") +
#  stat_cor(label.x = 3, label.y = 34) +
#  stat_regline_equation(label.x = 3, label.y = 32)+
#  labs(x="Temperature", y="N", title=w)+
#  theme(text = element_text(family = 'STSong'))
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
  nb<-readRDS(sprintf("../../Objects/Species_property/%s_property_compared.rda", g))
  
  nb<-nb%>%dplyr::select(nb_bio1_sd, nb_bio5_sd, nb_bio6_sd, nb_bio12_sd, nb_bio13_sd, nb_bio14_sd, N_CELL, sp)
  
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
  xlab("Niche breadth in temperature")+
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
  xlab("Niche breadth in precipitation")+
  ylab("Number of species")+
  facet_grid(exposure~Label)
p2
saveRDS(p2, "../../Figures/NB_hist_combined/nb_prec.rda")

p2<-p2+
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank())
ggsave(p2, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_prec.pdf"), width=12, height=6)
ggsave(p2, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_prec.png"), width=12, height=6)


pp<-ggarrange(p1, p2, ncol=1, nrow=2)
ggsave(pp, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_combined.pdf"), width=12, height=6)
ggsave(pp, filename=sprintf("../../Figures/N_Extinction/NB_Extinct/Extinction_hist_nb_combined.png"), width=12, height=6)
