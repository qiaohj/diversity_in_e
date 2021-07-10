
library(ggplot2)
library(ggpubr)
predicted_birds<-readRDS("../../Objects/estimate_disp_dist/models/predicted_birds.rda")
predicted_mammals<-readRDS("../../Objects/estimate_disp_dist/models/predicted_mammals.rda")

evaluated_metrics_birds<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")
evaluated_metrics_mammals<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_mammals.rda")

labels_birds<-c("max_dis~Diet"="Max disp ~ Diet",
                "max_dis~HWI"="Max disp ~ HWI",
                "max_dis~log_body_mass"="Max disp ~ Log(body mass)",
                "max_dis~HWI+log_body_mass"="Max disp ~ HWI + Log(body mass)",
                "max_dis~HWI+Diet"="Max disp ~ HWI + Diet",
                "max_dis~log_body_mass+Diet"="Max disp ~ Log(body mass) + Diet",
                "max_dis~HWI+log_body_mass+Diet"="Max disp ~ HWI + Log(body mass) + Diet")
color7<-rainbow(7)
colors_birds<-c("max_dis~Diet"= color7[1],
                  "max_dis~HWI"= color7[2],
                  "max_dis~log_body_mass"= color7[3],
                  "max_dis~HWI+log_body_mass"= color7[4],
                  "max_dis~HWI+Diet"= color7[5],
                  "max_dis~log_body_mass+Diet"= color7[6],
                  "max_dis~HWI+log_body_mass+Diet"= color7[7])

linetype_birds<-c("max_dis~Diet"= 1,
                "max_dis~HWI"= 2,
                "max_dis~log_body_mass"= 3,
                "max_dis~HWI+log_body_mass"= 4,
                "max_dis~HWI+Diet"= 5,
                "max_dis~log_body_mass+Diet"= 6,
                "max_dis~HWI+log_body_mass+Diet"= 7)

p1<-ggplot(predicted_birds)+
  geom_smooth(aes(x=max_dis, y=pred, 
                  color=factor(model), linetype=factor(formulas)),
              alpha=0.2)+
  geom_point(aes(x=max_dis, y=pred, 
                 color=factor(model), shape=factor(formulas)))+
  xlab("Empirical max natal dispersal distance")+
  ylab("Estimated max natal dispersal distance")+
  labs(linetype="Model formulas", shape="Model formulas", color="Model algorithms")+
  scale_shape_manual(values=linetype_birds,
    labels=labels_birds)+
  scale_linetype_manual(values=linetype_birds, labels=labels_birds)+
  ggtitle("Birds")+
  theme_bw()
p1

labels_mammals<-c("max_dis~Diet"= "Max disp ~ Diet",
                "max_dis~ForStrat"= "Max disp ~ Foraging strategy",
                "max_dis~log_body_mass"= "Max disp ~ Log(body mass)",
                "max_dis~ForStrat+log_body_mass"= "Max disp ~ Foraging strategy + Log(body mass)",
                "max_dis~ForStrat+Diet"= "Max disp ~ Foraging strategy + Diet",
                "max_dis~log_body_mass+Diet"= "Max disp ~ Log(body mass) + Diet",
                "max_dis~ForStrat+log_body_mass+Diet"= "Max disp ~ Foraging strategy + Log(body mass) + Diet")
color7<-rainbow(7)
colors_mammals<-c("max_dis~Diet"= color7[1],
                  "max_dis~ForStrat"= color7[2],
                  "max_dis~log_body_mass"= color7[3],
                  "max_dis~ForStrat+log_body_mass"= color7[4],
                  "max_dis~ForStrat+Diet"= color7[5],
                  "max_dis~log_body_mass+Diet"= color7[6],
                  "max_dis~ForStrat+log_body_mass+Diet"= color7[7])

linetype_mammals<-c("max_dis~Diet"= 1,
                  "max_dis~ForStrat"= 2,
                  "max_dis~log_body_mass"= 3,
                  "max_dis~ForStrat+log_body_mass"=4,
                  "max_dis~ForStrat+Diet"= 5,
                  "max_dis~log_body_mass+Diet"= 6,
                  "max_dis~ForStrat+log_body_mass+Diet"= 7)


p2<-ggplot(predicted_mammals)+
  geom_smooth(aes(x=max_dis, y=pred, 
                  color=factor(model), linetype=factor(formulas)),
              alpha=0.2)+
  geom_point(aes(x=max_dis, y=pred, 
                 color=factor(model), shape=factor(formulas)))+
  xlab("Empirical max natal dispersal distance")+
  ylab("Estimated max natal dispersal distance")+
  labs(linetype="Model formulas", shape="Model formulas", color="Model algorithms")+
  scale_shape_manual(values=linetype_mammals,
                     labels=labels_birds)+
  scale_linetype_manual(values=linetype_mammals,
                        labels=labels_birds)+
  #scale_color_discrete(guide=FALSE)+
  ggtitle("Mammals")+
  theme_bw()
p2

p<-ggarrange(p1, p2, ncol=1)
p
ggsave(p, filename="../../Figures/Estimate_Disp/Model_Results.png", width=10, height=12)

evaluated_metrics_birds[which(evaluated_metrics_birds$RMSE_Full==min(evaluated_metrics_birds$RMSE_Full)),]

p1<-ggplot(evaluated_metrics_birds)+
  geom_point(aes(x=factor(model), y=RMSE_Full, color=factor(formulas), size=cor))+
  xlab("Model algorithms")+
  ylab("RMSE")+
  labs(color="Model formulas", size="Correlation")+
  scale_color_manual(values=colors_birds,
                        labels=labels_birds)+
  ggtitle("Birds")+
  theme_bw()

p1

evaluated_metrics_mammals[which(evaluated_metrics_mammals$RMSE_Full==min(evaluated_metrics_mammals$RMSE_Full)),]
p2<-ggplot(evaluated_metrics_mammals)+
  geom_point(aes(x=model, y=RMSE_Full, color=formulas, size=cor))+
  xlab("Model algorithms")+
  ylab("RMSE")+
  labs(color="Model formulas", size="Correlation")+
  scale_color_manual(values=colors_mammals,
                     labels=labels_mammals)+
  ggtitle("Mammals")+
  theme_bw()

p2

p<-ggarrange(p1, p2, ncol=1)
p

ggsave(p, filename="../../Figures/Estimate_Disp/RMSE.png", width=8, height=10)


best_model_info<-evaluated_metrics_birds[which(evaluated_metrics_birds$RMSE_Full==min(evaluated_metrics_birds$RMSE_Full)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", tolower(best_model_info$model)))
best_model<-best_model[[best_model_info$formulas]]
birds_trait<-read.table("../../Data/Dispersal_distance/EltonTraits/Dataset_HWI.csv", sep=",", head=T, stringsAsFactors = F)
birds_trait[which(birds_trait$Diet=="fruit"), "Diet"]<-"plants"
birds_trait[which(birds_trait$Diet=="seeds"), "Diet"]<-"plants"
birds_trait[which(birds_trait$Diet=="nectar"), "Diet"]<-"plants"
birds_trait[which(birds_trait$Diet=="scav"), "Diet"]<-"omnivore"
birds_trait[which(is.na(birds_trait$Diet)), "Diet"]<-"omnivore"

birds_trait[which(is.na(birds_trait$Migration_3)), "Migration_3"]<-"unknown"
birds_trait<-data.table(birds_trait)
cols<-c("HWI", "log_body_mass", "Diet", "Migration_3", "iucn_name")
new_df_birds<-birds_trait[,..cols]
new_df_birds[is.na(Migration_3)]
new_df_birds<-new_df_birds[(!is.na(HWI))&(!is.na(log_body_mass))&(!is.na(Diet))]
new_df_birds$estimated_disp<-predict(best_model, new_df_birds)

unique(new_df_birds$Diet)
table(new_df_birds$Migration_3)
p1<-ggplot(new_df_birds)+
  geom_boxplot(aes(x=Migration_3, y=estimated_disp))+
  xlab("Migration type")+
  ylab("Estimated dispersal distance")+
  ggtitle("Birds")+
  theme_bw()
p1
ggsave(p1, filename="../../Figures/Estimate_Disp/Predicted_birds_mig.png", width=8, height=6)

p1<-ggplot(new_df_birds)+
  geom_boxplot(aes(x=Diet, y=estimated_disp))+
  xlab("Diet type")+
  ylab("Estimated dispersal distance")+
  ggtitle("Birds")+
  theme_bw()
p1



best_model_info<-evaluated_metrics_mammals[which(evaluated_metrics_mammals$RMSE_Full==min(evaluated_metrics_mammals$RMSE_Full)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_mammals.rda", tolower(best_model_info$model)))
best_model<-best_model[[best_model_info$formulas]]

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
mammals_trait$log_body_mass<-log(mammals_trait$BodyMass.Value)



cols<-c("Scientific", "ForStrat.Value", "log_body_mass", "diet_type")

new_df_mammals<-data.table(mammals_trait)[,..cols]
colnames(new_df_mammals)<-c("Scientific", "ForStrat", "log_body_mass", "Diet")
new_df_mammals<-new_df_mammals[(!is.na(ForStrat))&(!is.na(log_body_mass))&(!is.na(Diet))]
new_df_mammals<-new_df_mammals[!(ForStrat %in% c("A", "M"))]
new_df_mammals$estimated_disp<-predict(best_model, new_df_mammals)
new_df_mammals$Diet
model_df_mammals$ForStrat
unique(new_df_mammals$Diet)
summary(new_df_mammals$estimated_disp)

model_df_mammals$estimated_disp<-predict(best_model, model_df_mammals)

p2<-ggplot(new_df_mammals)+
  geom_boxplot(aes(x=Diet, y=estimated_disp))+
  xlab("Diet type")+
  ylab("Estimated dispersal distance")+
  ggtitle("Mammals")+
  theme_bw()
p2
p<-ggarrange(p1, p2, ncol=1)
p

ggsave(p, filename="../../Figures/Estimate_Disp/Predicted_diet.png", width=8, height=8)
