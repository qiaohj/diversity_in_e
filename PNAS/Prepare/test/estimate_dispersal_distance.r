library(data.table)
library(ggplot2)
library(randomForest)
library(ggpubr)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
#Mammals
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
#Mammals Trait
mammals_trait<-read.table("../../Data/Dispersal distance/EltonTraits/MamFuncDat.txt", sep="\t", head=T, stringsAsFactors = F)
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

unique(mammals_trait$ForStrat.Value)
#M - marine: 126
#G - ground level, including aquatic foraging (see ForStrat-Comment): 3098
#S - scansorial: 370
#Ar- arboreal: 1069
#A - aerial: 737


plot(lmm)

table(mammals_trait$ForStrat.Value)

table(mammals_trait$diet_type)

#long
#Dcarn=40.7*M^0.81
#Dherb_omn=3.31*M^0.65

#median
#Dcarn=3.45*M^0.89
#Dherb_omn=1.45*M^0.54


mammals_trait$BodyMass.Value<-as.numeric(mammals_trait$BodyMass.Value)
mammals_trait$BodyMass.Value<-mammals_trait$BodyMass.Value/1000
mammals_trait$estimated_disp<-ifelse(mammals_trait$diet_type=="carnivore", 
                                     3.45*mammals_trait$BodyMass.Value^0.89,
                                     1.45*mammals_trait$BodyMass.Value^0.54)
mammals_trait<-mammals_trait[!is.na(mammals_trait$estimated_disp),]
mammals_trait[which(mammals_trait$estimated_disp==max(mammals_trait$estimated_disp)),]

hist(mammals_trait$estimated_disp[which(mammals_trait$estimated_disp<5e3)])
mammals_trait<-data.table(mammals_trait)

disp_df<-merge(empirical_disp, mammals_trait, by.x="sp", by.y="Scientific", all.x=T, all.y=F)
disp_df<-disp_df[!is.na(disp_df$estimated_disp),]
plot(disp_df$estimated_disp, disp_df$max_dis)
ggplot(disp_df)+geom_point(aes(x=max_dis, y=estimated_disp, color=factor(diet_type)))+
  facet_wrap(~source, scale="free", nrow=2)+theme_bw()
cor(disp_df$estimated_disp, disp_df$max_dis)

plot(disp_df$estimated_disp, disp_df$max_dis)
disp_df[which(disp_df$BodyMass.Value==max(disp_df$BodyMass.Value)),]

disp_df[which(disp_df$sp=="Alces alces"),]

glm(max_dis~sqrt(BodyMass.Value), data=disp_df[which(disp_df$source=="Santini, L., et. al., 2013"),])
gg<-glm(max_dis~sqrt(BodyMass.Value), data=disp_df[which(disp_df$source=="Sutherland, G. D., et. al., 2000"),])
pred<-predict(gg)

plot(disp_df[source=="Sutherland, G. D., et. al., 2000"]$max_dis, pred)
cor(disp_df[source=="Sutherland, G. D., et. al., 2000"]$max_dis, pred)
