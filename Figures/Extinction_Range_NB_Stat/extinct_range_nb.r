library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
df<-NULL
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  nb<-readRDS(sprintf("../../Objects_Full_species/Species_property/%s_property.rda", g))
  
  nb<-nb%>%dplyr::select(nb_TEMP_sd, nb_PR_sd, N_CELL, sp)
  nb<-as.data.frame(nb)
  ttt<-2
  for (threshold in c(1, 5)){
    extinct<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
    extinct<-extinct%>%dplyr::distinct(sp, GCM, SSP, extinct_year, dispersal)
    extinct[is.infinite(extinct$extinct_year), "extinct_year"]<-2020
    df1<-full_join(nb, extinct, by="sp")
    df1<-df1%>%filter(!is.na(extinct_year))
    df1<-df1%>%dplyr::filter(N_CELL>ttt)
    df1$threshold<-threshold
    df1$group<-g
    df<-bind(df, df1)
  }
}
plot(df1$extinct_year, df1$N_CELL)


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

df_se<-df%>%dplyr::group_by(nb_TEMP_sd, nb_PR_sd, N_CELL, sp, dispersal, threshold, group)%>%
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
df_se$exposure<-ifelse(df_se$threshold==1, " no exposure", "5-year exposure")
df_se$da<-ifelse(df_se$dispersal==1, "with dispersal", "no dispersal")
labels<-c("Distribution range", "Range of temperature", "Range of precipitation")
vars<-c("N_CELL", "nb_TEMP_sd", "nb_PR_sd")
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
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, N_CELL)%>%select(sp, group, N_CELL)
    colnames(df_hist)[3]<-"V"
  }
  if (i==2){
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, nb_TEMP_sd)%>%select(sp, group, nb_TEMP_sd)
    colnames(df_hist)[3]<-"V"
  }
  if (i==3){
    df_hist<-df_se%>%ungroup()%>%distinct(sp, group, nb_PR_sd)%>%select(sp, group, nb_PR_sd)
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
  ggsave(ppp, filename=sprintf("../../Figures_Full_species/Extinction_Range_NB_Stat/%s.png", var),
         width=15, height=4)
  ggsave(ppp, filename=sprintf("../../Figures_Full_species/Extinction_Range_NB_Stat/%s.pdf", var),
         width=15, height=4)
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
dispersals<-c(0, 1)

df<-NULL
for (g in c("Amphibians", "Birds", "Mammals", "Reptiles")){
  print(g)
  nb<-readRDS(sprintf("../../Objects_Full_species/Species_property/%s_property.rda", g))
  
  nb<-nb%>%dplyr::select(nb_TEMP_sd, nb_PR_sd, N_CELL, sp)
  nb<-as.data.frame(nb)
  full_sp<-expand.grid(GCM=GCMs, SSP=SSPs, sp=unique(nb$sp), dispersal=dispersals, stringsAsFactors = F)
  full_sp_with_nb<-full_join(full_sp, nb, by="sp")
  ttt<-2
  for (threshold in c(1, 5)){
    extinct<-readRDS(sprintf("../../Objects_Full_species/when_where_extinction_%d/%s.rda", threshold, g))
    extinct<-extinct%>%dplyr::distinct(sp, GCM, SSP, extinct_year, dispersal)
    extinct[is.infinite(extinct$extinct_year), "extinct_year"]<-2020
    df1<-full_join(full_sp_with_nb, extinct, by=c("sp", "GCM", "SSP", "dispersal"))
    df1<-df1%>%dplyr::filter(N_CELL>ttt)
    df1$threshold<-threshold
    df1$group<-g
    df<-bind(df, df1)
  }
}





df$is_extinct<-ifelse(is.na(df$extinct_year), "NO", "YES")
df$exposure<-ifelse(df$threshold==1, " no exposure", "5-year exposure")
df$da<-ifelse(df$dispersal==1, "with dispersal", "no dispersal")
df$Label<-paste(df$SSP, df$da)
p<-ggplot(df)+
  geom_histogram(aes(x=nb_TEMP_sd), fill=colors_black[4], bins=50)+
  geom_histogram(data=df%>%dplyr::filter(is_extinct=="YES"), 
                 aes(x=nb_TEMP_sd), fill=colors_red[9], bins=50)+
  theme_bw()+
  xlab("Niche breadth in temperature)")+
  ylab("Number of species")+
  facet_grid(exposure~Label, scale="free")
p
ggsave(p, filename=sprintf("../../Figures_Full_species/N_Extinction/NB_Extinct/Extinction_hist_nb_temp_%d.pdf", ttt), width=12, height=6)
ggsave(p, filename=sprintf("../../Figures_Full_species/N_Extinction/NB_Extinct/Extinction_hist_nb_temp_%d.png", ttt), width=12, height=6)


p<-ggplot(df)+
  geom_histogram(aes(x=nb_PR_sd), fill=colors_black[4], bins=50)+
  geom_histogram(data=df%>%dplyr::filter(is_extinct=="YES"), 
                 aes(x=nb_PR_sd), fill=colors_red[9], bins=50)+
  theme_bw()+
  xlab("Niche breadth in precipitation)")+
  ylab("Number of species")+
  facet_grid(exposure~Label)
p
ggsave(p, filename=sprintf("../../Figures_Full_species/N_Extinction/NB_Extinct/Extinction_hist_nb_prec_%d.pdf", ttt), width=12, height=6)
ggsave(p, filename=sprintf("../../Figures_Full_species/N_Extinction/NB_Extinct/Extinction_hist_nb_prec_%d.png", ttt), width=12, height=6)
