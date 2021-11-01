library(data.table)
library(ggplot2)
library(randomForest)
library(ggpubr)
rm(list=ls())
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#pnas
pnas<-read.table("../../Data/Dispersal distance/pnas/pnas-2012.csv", sep=",", head=T, stringsAsFactors = F)
pnas$sp<-sprintf("%s %s", pnas$Species, pnas$S2)
colnames(pnas)[c(5,7,8)]<-c("Bodybass", "generationlength", "estimated_disp")
pnas<-data.table(pnas)
family<-unique(pnas$Family)[1]
pnas$Bodybass
print(family)
item<-pnas[Family==family]
rfm_with_family<-randomForest(generationlength~Bodybass+Family+DietType, data=pnas)
rfm_with_family$rsq
actual<-pnas$generationlength
predicted<-predict(rfm_with_family, pnas)
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
pnas$pred_rf_generation_length_with_family<-predicted
p1<-ggplot(pnas)+geom_point(aes(x=generationlength, y=pred_rf_generation_length_with_family, color=factor(Family)))+
  facet_wrap(~DietType, scale="free", nrow=1)+theme_bw()+theme(legend.position = "none")+
  ggtitle("Estimate the generation length with random forest with family")

lmm_with_family<-lm(generationlength~Bodybass+Family+DietType, data=pnas)
summary(lmm_with_family)
pred<-predict(lmm_with_family, pnas)
pnas$pred_lm_generation_length_with_family<-pred
p2<-ggplot(pnas)+geom_point(aes(x=generationlength, y=pred_lm_generation_length_with_family, color=factor(Family)))+
  facet_wrap(~DietType, scale="free", nrow=1)+theme_bw()+theme(legend.position = "none")+
  ggtitle("Estimate the generation length with linear model with family")


p<-ggpubr::ggarrange(p1, p2, nrow=2)
p
ggsave(p, filename = "../../Figures/Estimate_Disp/estimate_generation_length_with_family.png", width=10, height=8)


rfm_without_family<-randomForest(generationlength~Bodybass+DietType, data=pnas)
rfm_without_family$rsq
actual<-pnas$generationlength
predicted<-predict(rfm_without_family, pnas)
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
pnas$pred_rf_generation_length_without_family<-predicted
p1<-ggplot(pnas)+geom_point(aes(x=generationlength, y=pred_rf_generation_length_without_family, color=factor(Family)))+
  facet_wrap(~DietType, scale="free", nrow=1)+theme_bw()+theme(legend.position = "none")+
  ggtitle("Estimate the generation length with random forest without family")

lmm_without_family<-lm(generationlength~Bodybass+DietType, data=pnas)
summary(lmm_without_family)
pred<-predict(lmm_without_family, pnas)
pnas$pred_lm_generation_length_without_family<-pred
p2<-ggplot(pnas)+geom_point(aes(x=generationlength, y=pred_lm_generation_length_without_family, color=factor(Family)))+
  facet_wrap(~DietType, scale="free", nrow=1)+theme_bw()+theme(legend.position = "none")+
  ggtitle("Estimate the generation length with linear model without family")


p<-ggpubr::ggarrange(p1, p2, nrow=2)
p
ggsave(p, filename = "../../Figures/Estimate_Disp/estimate_generation_length_without_family.png", width=10, height=8)


pnas$estimated_disp_lm_without_family<-(ifelse(pnas$DietType=="carnivore", 
                                     3.45*pnas$Bodybass^0.89,
                                     1.45*pnas$Bodybass^0.54)/pnas$pred_lm_generation_length_without_family)*1.5

pnas$estimated_disp_lm_with_family<-(ifelse(pnas$DietType=="carnivore", 
                                               3.45*pnas$Bodybass^0.89,
                                               1.45*pnas$Bodybass^0.54)/pnas$pred_lm_generation_length_with_family)*1.5

p1<-ggplot(pnas)+geom_point(aes(x=estimated_disp, y=estimated_disp_lm_without_family, color=factor(Family)))+
  ggtitle("Estimated dispersal distance with linear model and without family")+theme_bw()+theme(legend.position = "none")
p2<-ggplot(pnas)+geom_point(aes(x=estimated_disp, y=estimated_disp_lm_with_family, color=factor(Family)))+
  ggtitle("Estimated dispersal distance with linear model and family")+theme_bw()+theme(legend.position = "none")

p<-ggpubr::ggarrange(p1, p2, nrow=1)
p
ggsave(p, filename = "../../Figures/Estimate_Disp/estimate_disp_distance_without_family.png", width=10, height=4)



#Mammals Trait
mammals_trait<-read.table("../../Data/Dispersal_distance/EltonTraits/MamFuncDat.txt", sep="\t", head=T, stringsAsFactors = F)
mammals_trait_bak<-mammals_trait
DietType<-c("omnivore", "herbivore", "carnivore")
herbivore<-c("Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO")
carnivore<-c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk")
for (col in c(herbivore, carnivore)){
  mammals_trait[, col]<-as.numeric(mammals_trait[, col])
}
#mammals_trait<-data.table(mammals_trait)
mammals_trait$herbivore_score<-rowSums(mammals_trait[, herbivore], na.rm = T)
mammals_trait$carnivore_score<-rowSums(mammals_trait[, carnivore], na.rm = T)
mammals_trait$diet_type<-"omnivore"
mammals_trait[which(mammals_trait$herbivore_score==100), "diet_type"]<-"herbivore"
mammals_trait[which(mammals_trait$carnivore_score==100), "diet_type"]<-"carnivore"

mammals_trait$BodyMass.Value<-as.numeric(mammals_trait$BodyMass.Value)
mammals_trait$BodyMass.Value<-mammals_trait$BodyMass.Value/1000

unique(mammals_trait$ForStrat.Value)
#M - marine: 126
#G - ground level, including aquatic foraging (see ForStrat-Comment): 3098
#S - scansorial: 370
#Ar- arboreal: 1069
#A - aerial: 737
mammals_trait_g<-data.table(mammals_trait)[ForStrat.Value=="G"]
mammals_trait_g_with_family<-mammals_trait_g[MSWFamilyLatin %in% pnas$Family]
mammals_trait_g_without_family<-mammals_trait_g[!(MSWFamilyLatin %in% pnas$Family)]
#Bodybass+Family+DietType
cols = c("BodyMass.Value", "MSWFamilyLatin", "diet_type")
item<-mammals_trait_g_with_family[, mget(cols)]
colnames(item)<-c("Bodybass", "Family", "DietType")
mammals_trait_g_with_family$generation_length<-predict(lmm_with_family, item)

mammals_trait_g_with_family$estimated_disp<-(ifelse(mammals_trait_g_with_family$diet_type=="carnivore", 
                                           3.45*mammals_trait_g_with_family$BodyMass.Value^0.89,
                                           1.45*mammals_trait_g_with_family$BodyMass.Value^0.54)/
                                             mammals_trait_g_with_family$generation_length)*1.5

#Bodybass+DietType
cols = c("BodyMass.Value", "diet_type")
item<-mammals_trait_g_without_family[, mget(cols)]
colnames(item)<-c("Bodybass", "DietType")
mammals_trait_g_without_family$generation_length<-predict(lmm_without_family, item)

mammals_trait_g_without_family$estimated_disp<-(ifelse(mammals_trait_g_without_family$diet_type=="carnivore", 
                                                    3.45*mammals_trait_g_without_family$BodyMass.Value^0.89,
                                                    1.45*mammals_trait_g_without_family$BodyMass.Value^0.54)/
                                               mammals_trait_g_without_family$generation_length)*1.5


mammals_trait_g_full<-rbind(mammals_trait_g_with_family, mammals_trait_g_without_family)


#Empirical distance
#Sutherland
Sutherland<-read.table("../../Data/Dispersal distance/Bird_Mammals_Dispersal_Distance/Mammals.csv", head=T, sep=",")
colnames(Sutherland)<-c("sp", "body_mass", "obs_type", "dispersal_distance_median", 
                        "dispersal_distance_max", "x1", "x2", "n", "source")

Sutherland$max_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_max))))
Sutherland$max_float<-as.numeric(Sutherland$max_str)

Sutherland$median_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_median))))
Sutherland$median_float<-as.numeric(Sutherland$median_str)

Sutherland<-data.table(Sutherland)
Sutherland_se<-Sutherland[, .(max_dis=max(c(max_float,median_float), na.rm=T)), by=list(sp)]
Sutherland_se$source<-"Sutherland, G. D., et. al., 2000"

#Santini
Santini<-read.table("../../Data/Dispersal distance/Mammals_Distance/data.csv", head=T, sep=",", stringsAsFactors = F)
Santini_max_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Maximum)), c(" "))
Santini_mean_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Mean)), c(" "))
Santini_median_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Median)), c(" "))
Santini$max_dis<-NA
i=1
for (i in c(1:nrow(Santini))){
  print(paste(i, nrow(Santini)))
  dist<-c(Santini_max_str[[i]], Santini_mean_str[[i]], Santini_median_str[[i]])
  Santini[i, "max_dis"]<-max(as.numeric(dist), na.rm=T)
}
Santini<-data.table(Santini)
Santini_se<-Santini[, .(max_dis=max(max_dis, na.rm=T)), by=list(Species)]
Santini_se<-Santini_se[Species!=""]
colnames(Santini_se)[1]<-"sp"
Santini_se$source<-"Santini, L., et. al., 2013"

empirical_disp<-rbind(Santini_se, Sutherland_se)
empirical_disp<-empirical_disp
empirical_disp_with_estimate<-merge(empirical_disp, mammals_trait_g_full, by.x="sp", by.y="Scientific", all.x=T, all.y=F)


empirical_disp_with_estimate<-empirical_disp_with_estimate[!is.na(empirical_disp_with_estimate$estimated_disp)]
p<-ggplot(empirical_disp_with_estimate)+geom_point(aes(x=max_dis, y=estimated_disp))
ggsave(p, filename = "../../Figures/Estimate_Disp/estimate_disp_distance_for_others.png", width=10, height=4)

