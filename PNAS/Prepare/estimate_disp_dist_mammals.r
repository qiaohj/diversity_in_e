library(randomForest)
library(caret)
library(rgdal)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")


#Mammals
#Sutherland
Sutherland<-read.table("../../Data/Dispersal_distance/Bird_Mammals_Dispersal_Distance/Mammals.csv", head=T, sep=",")
colnames(Sutherland)<-c("sp", "body_mass", "obs_type", "dispersal_distance_median", 
                        "dispersal_distance_max", "x1", "x2", "n", "source")

Sutherland$max_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_max))))
Sutherland$max_float<-as.numeric(Sutherland$max_str)

Sutherland$median_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_median))))
Sutherland$median_float<-as.numeric(Sutherland$median_str)

Sutherland<-data.table(Sutherland)
Sutherland_se<-Sutherland[, .(max_dis=max(c(max_float,median_float), na.rm=T),
                              median_dis=max(median_float, na.rm=T)), by=list(sp)]
Sutherland_se$source<-"Sutherland, G. D., et. al., 2000"
length(unique(Sutherland_se$sp))
dim(unique((Sutherland_se[!is.infinite(Sutherland_se$median_dis), "sp"])))
dim(unique((Sutherland_se[!is.infinite(Sutherland_se$max_dis), "sp"])))

#Santini
Santini<-read.table("../../Data/Dispersal_distance/Mammals_Distance/data.csv", head=T, sep=",", stringsAsFactors = F)
Santini_max_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Maximum)), c(" "))
Santini_mean_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Mean)), c(" "))
Santini_median_str<-strsplit(gsub("f", "", gsub("m", "", Santini$Median)), c(" "))
Santini$max_dis<-NA
Santini$median_dis<-NA

i=1
for (i in c(1:nrow(Santini))){
  print(paste(i, nrow(Santini)))
  dist<-c(Santini_max_str[[i]], Santini_mean_str[[i]], Santini_median_str[[i]])
  Santini[i, "max_dis"]<-max(as.numeric(dist), na.rm=T)
  dist<-c(Santini_median_str[[i]])
  if (length(dist)>0){
    Santini[i, "median_dis"]<-max(as.numeric(dist), na.rm=T)
  }
}
Santini<-data.table(Santini)
Santini_se<-Santini[, .(max_dis=max(max_dis, na.rm=T),
                        median_dis=max(median_dis, na.rm = T)), by=list(Species)]
Santini_se<-Santini_se[Species!=""]
colnames(Santini_se)[1]<-"sp"
Santini_se$source<-"Santini, L., et. al., 2013"
length(unique(Santini_se$sp))
dim(unique((Santini_se[!is.infinite(Santini_se$median_dis), "sp"])))
dim(unique((Santini_se[!is.infinite(Santini_se$max_dis), "sp"])))



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

unique(mammals_trait$ForStrat.Value)
#M - marine: 126
#G - ground level, including aquatic foraging (see ForStrat-Comment): 3098
#S - scansorial: 370
#Ar- arboreal: 1069
#A - aerial: 737


table(mammals_trait$ForStrat.Value)

table(mammals_trait$diet_type)

mammals_trait$BodyMass.Value<-as.numeric(mammals_trait$BodyMass.Value)
mammals_trait$log_body_mass<-log(mammals_trait$BodyMass.Value)
hist(mammals_trait$log_body_mass)
empirical_disp<-rbind(Santini_se, Sutherland_se)

length(unique(empirical_disp$sp))
#59/152
dim(empirical_disp[!is.infinite(max_dis)])
empirical_disp[which(empirical_disp$sp=="Aotus azarai"), "sp"]<-"Aotus azarae"
empirical_disp[which(empirical_disp$sp=="Cervus canadensis"), "sp"]<-"Cervus canadensis"
empirical_disp[which(empirical_disp$sp=="Clethrionomys glareolus"), "sp"]<-"Myodes glareolus"
empirical_disp[which(empirical_disp$sp=="Didelphis virginianus"), "sp"]<-"Bubo virginianus"
empirical_disp[which(empirical_disp$sp=="Felis concolor"), "sp"]<-"Puma concolor"
empirical_disp[which(empirical_disp$sp=="Giraffa camelopardalis thornicrofti"), "sp"]<-"Giraffa camelopardalis"
empirical_disp[which(empirical_disp$sp=="Lepus europea"), "sp"]<-"Lepus europaeus"
empirical_disp[which(empirical_disp$sp=="Lutra canadensis"), "sp"]<-"Lontra canadensis"
empirical_disp[which(empirical_disp$sp=="Lynx lynx **"), "sp"]<-"Lynx lynx"
empirical_disp[which(empirical_disp$sp=="Microtus townsendi"), "sp"]<-"Microtus townsendii"
empirical_disp[which(empirical_disp$sp=="Mustela putorius furo"), "sp"]<-"Mustela putorius"
empirical_disp[which(empirical_disp$sp=="Mustela vison"), "sp"]<-"Neovison vison"
empirical_disp[which(empirical_disp$sp=="Odocoileus hemionus columbianus"), "sp"]<-"Odocoileus hemionus"
empirical_disp[which(empirical_disp$sp=="Odocoileus hemionus hemionus"), "sp"]<-"Odocoileus hemionus"
empirical_disp[which(empirical_disp$sp=="Panthera leo persica"), "sp"]<-"Panthera leo"
empirical_disp[which(empirical_disp$sp=="Perognathus formosa"), "sp"]<-"Chaetodipus formosus"
#empirical_disp[which(empirical_disp$sp=="Peromyscus longicaudus"), "sp"]<-"xxxxx"
empirical_disp[which(empirical_disp$sp=="Phascogale tapotafa"), "sp"]<-"Phascogale tapoatafa"
empirical_disp[which(empirical_disp$sp=="Spermophilus tridecemliniatus"), "sp"]<-"Ictidomys tridecemlineatus"
#empirical_disp[which(empirical_disp$sp=="Spermophilus leucopus"), "sp"]<-"xxxxx"
empirical_disp[which(empirical_disp$sp=="Sus scrofa **"), "sp"]<-"Sus scrofa"
empirical_disp[which(empirical_disp$sp=="Sylvilagus bachmani ubericolor"), "sp"]<-"Sylvilagus bachmani"
empirical_disp[which(empirical_disp$sp=="Vulpes mactrotis"), "sp"]<-"Vulpes macrotis"
empirical_disp<-empirical_disp[, .(max_dis=max(max_dis)), by=list(sp)]



model_df_mammals<-merge(empirical_disp, mammals_trait, by.x="sp", by.y="Scientific", all.x=T, all.y=F)
model_df_mammals<-model_df_mammals[!is.na(MSW3_ID)]
model_df_mammals<-model_df_mammals[!is.infinite(max_dis)]
model_df_mammals[sp=="Puma concolor"]
table(model_df_mammals$sp)
table(model_df_mammals$ForStrat.Value)
table(model_df_mammals$diet_type)

cols<-c("sp", "max_dis", "ForStrat.Value", "diet_type", "log_body_mass")
model_df_mammals<-model_df_mammals[, ..cols]
colnames(model_df_mammals)<-c("sp", "max_dis", "ForStrat", "Diet", "log_body_mass")
model_df_mammals<-unique(model_df_mammals)
formulas<-c("max_dis~Diet",
            "max_dis~ForStrat",
            "max_dis~log_body_mass",
            "max_dis~ForStrat+log_body_mass",
            "max_dis~ForStrat+Diet",
            "max_dis~log_body_mass+Diet",
            "max_dis~ForStrat+log_body_mass+Diet")


if (T){
  ########For random forest##########################
  # define training control for 10-fold Cross Validation
  evaluated_metrics_all<-NULL
  predicted_all<-NULL
  seeds_index<-1000
  setSeeds <- function(method = "cv", numbers = 1, repeats = 1, tunes = NULL, seed = 1) {
    #B is the number of resamples and integer vector of M (numbers + tune length if any)
    B <- if (method == "cv") numbers
    else if(method == "repeatedcv") numbers * repeats
    else NULL
    
    if(is.null(length)) {
      seeds <- NULL
    } else {
      set.seed(seed = seed)
      seeds <- vector(mode = "list", length = B)
      seeds <- lapply(seeds, function(x) sample.int(n = 1000000, size = numbers + ifelse(is.null(tunes), 0, tunes)))
      seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
    }
    # return seeds
    seeds
  }
  
  
  i=6
  rf_no_rank_list<-list()
  glm_no_rank_list<-list()
  svm_no_rank_list<-list()
  brnn_no_rank_list<-list()
  rpart_no_rank_list<-list()
  evtree_no_rank_list<-list()
  
  for (i in c(1:length(formulas))){
    f<-formulas[i]
    print(f)
    # train the model
    #df_with_family$Migration_3<-as.factor(df_with_family$Migration_3)
    tunegrid <- expand.grid(.mtry=c(1:5))
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    rf_no_rank <- train(as.formula(f), data=model_df_mammals, 
                        trControl=train_control, method="rf", 
                        tuneGrid=tunegrid, ntree=1000)
    
    rf_no_rank_list[[f]]<-rf_no_rank
    # summarize results
    print(rf_no_rank)
    predicted<-model_df_mammals
    predicted$model<-"RF"
    predicted$formulas<-f
    predicted$clade<-"Mammal"
    predicted$pred<-predict(rf_no_rank, model_df_mammals)
    cor<-cor(model_df_mammals$max_dis, predicted$pred)
    
    evaluated_metrics<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"RF"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Mammal"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    if (is.null(evaluated_metrics_all)){
      evaluated_metrics_all<-evaluated_metrics  
      predicted_all<-predicted
    }else{
      evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
      predicted_all<-rbind(predicted_all, predicted)
    }
    
    #######################end of random forest ################################
    
    
    ########For GLM##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    glm_no_rank <- train(as.formula(f), 
                         data=model_df_mammals, trControl=train_control, 
                         method="glm")
    glm_no_rank_list[[f]]<-glm_no_rank
    # summarize results
    print(glm_no_rank)
    predicted<-model_df_mammals
    predicted$model<-"GLM"
    predicted$formulas<-f
    predicted$clade<-"Mammal"
    predicted$pred<-predict(glm_no_rank, model_df_mammals)
    cor<-cor(model_df_mammals$max_dis, predicted$pred)
    
    evaluated_metrics<-glm_no_rank$results[which(glm_no_rank$results$RMSE==min(glm_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"GLM"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Mammal"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of GLM ################################
    
    ########For Support Vector Machines with Linear Kernel ##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    svm_no_rank <- train(as.formula(f), data=model_df_mammals, trControl=train_control, 
                         method="svmLinear")
    svm_no_rank_list[[f]]<-svm_no_rank
    # summarize results
    print(svm_no_rank)
    predicted<-model_df_mammals
    predicted$model<-"SVM"
    predicted$formulas<-f
    predicted$clade<-"Mammal"
    predicted$pred<-predict(svm_no_rank, model_df_mammals)
    cor<-cor(model_df_mammals$max_dis, predicted$pred)
    
    evaluated_metrics<-svm_no_rank$results[which(svm_no_rank$results$RMSE==min(svm_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"SVM"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Mammal"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of SVM ################################
    
    
    
    ########For Bayesian Regularized Neural Networks##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    brnn_no_rank <- train(as.formula(f), data=model_df_mammals, trControl=train_control, 
                          method="brnn")
    brnn_no_rank_list[[f]]<-brnn_no_rank
    # summarize results
    print(brnn_no_rank)
    predicted<-model_df_mammals
    predicted$model<-"BRNN"
    predicted$formulas<-f
    predicted$clade<-"Mammal"
    predicted$pred<-predict(brnn_no_rank, model_df_mammals)
    cor<-cor(model_df_mammals$max_dis, predicted$pred)
    
    evaluated_metrics<-brnn_no_rank$results[which(brnn_no_rank$results$RMSE==min(brnn_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"BRNN"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Mammal"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of BRNN ################################
    
    ########For CART##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    rpart_no_rank <- train(as.formula(f), data=model_df_mammals, trControl=train_control, 
                           method="rpart")
    rpart_no_rank_list[[f]]<-rpart_no_rank
    # summarize results
    print(rpart_no_rank)
    predicted<-model_df_mammals
    predicted$model<-"CART"
    predicted$formulas<-f
    predicted$clade<-"Mammal"
    predicted$pred<-predict(rpart_no_rank, model_df_mammals)
    cor<-cor(model_df_mammals$max_dis, predicted$pred)
    
    evaluated_metrics<-rpart_no_rank$results[which(rpart_no_rank$results$RMSE==min(rpart_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"CART"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Mammal"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_mammals$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of CART ################################
    
  }
  saveRDS(rf_no_rank_list, "../../Objects/estimate_disp_dist/models/rf_no_rank_mammals.rda")
  
  saveRDS(glm_no_rank_list, "../../Objects/estimate_disp_dist/models/glm_no_rank_mammals.rda")
  
  saveRDS(svm_no_rank_list, "../../Objects/estimate_disp_dist/models/svm_no_rank_mammals.rda")
  
  saveRDS(brnn_no_rank_list, "../../Objects/estimate_disp_dist/models/brnn_no_rank_mammals.rda")
  
  saveRDS(rpart_no_rank_list, "../../Objects/estimate_disp_dist/models/rpart_no_rank_mammals.rda")
  
  saveRDS(evtree_no_rank_list, "../../Objects/estimate_disp_dist/models/evtree_no_rank_mammals.rda")
  
  saveRDS(evaluated_metrics_all, "../../Objects/estimate_disp_dist/models/evaluated_metrics_mammals.rda")
  
  saveRDS(predicted_all, "../../Objects/estimate_disp_dist/models/predicted_mammals.rda")
  
}

library(ggplot2)
predicted_all<-readRDS("../../Objects/estimate_disp_dist/models/predicted_mammals.rda")

evaluated_metrics_all<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_mammals.rda")
evaluated_metrics_all[evaluated_metrics_all$RMSE==min(evaluated_metrics_all$RMSE), ]
p<-ggplot(predicted_all)+
  geom_smooth(aes(x=max_dis, y=pred, 
                  color=factor(model), linetype=factor(formulas)),
              alpha=0.2)+
  geom_point(aes(x=max_dis, y=pred, 
                 color=factor(model), shape=factor(formulas)))+
  theme_bw()
p
ggsave(p, filename="../../Figures/Estimate_Disp/Model_Results_mammals.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=RMSE, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/RMSE_mammals.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=Rsquared, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/Rsquared_mammals.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=MAE, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/MAE_mammals.png", width=8, height=6)


evaluated_metrics_all[which(evaluated_metrics_all$cor==max(evaluated_metrics_all$cor, na.rm=T)),]
evaluated_metrics_all[which(evaluated_metrics_all$Rsquared==max(evaluated_metrics_all$Rsquared, na.rm=T)),]

threshold_rmse<-quantile(evaluated_metrics_all$RMSE_Full, 0.1, na.rm=T)
best_model_info<-
  evaluated_metrics_all[
    which(evaluated_metrics_all$RMSE_Full<=threshold_rmse),]

#best_model_info<-best_model_info[2,]
best_model_info<-evaluated_metrics_all[which(evaluated_metrics_all$RMSE_Full==min(evaluated_metrics_all$RMSE_Full)),]
#best_model_info<-evaluated_metrics_all[which(evaluated_metrics_all$cor==max(evaluated_metrics_all$cor, na.rm=T)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_mammals.rda", tolower(best_model_info$model)))
#best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))

best_model<-best_model[[best_model_info$formulas]]
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

p1<-ggplot(model_df_mammals)+geom_point(aes(x=max_dis, y=estimated_disp, color=factor(Diet)))+theme_bw()
p1
p2<-ggplot(model_df_mammals)+geom_point(aes(x=log_body_mass, y=max_dis, color=factor(Diet)))+theme_bw()
p2
p<-ggpubr::ggarrange(p2, p1, nrow=2)
ggsave(p, filename="../../Figures/Estimate_Disp/emperical_predicted_mammals.png", width=8, height=8)


p<-ggplot(new_df_mammals)+geom_boxplot(aes(y=estimated_disp, color=factor(Diet)))+theme_bw()
p
ggsave(p, filename="../../Figures/Estimate_Disp/Predicted_mammals.png", width=8, height=6)

p1<-ggplot(new_df_mammals)+geom_point(aes(x=log_body_mass, y=estimated_disp, color=factor(Diet)))+theme_bw()
p1
p2<-ggplot(new_df_birds)+geom_point(aes(x=HWI, y=estimated_disp, color=factor(Diet)))+theme_bw()
p2
p<-ggpubr::ggarrange(p1, p2, nrow=2)
p
ggsave(p, filename="../../Figures/Estimate_Disp/Final_results.png", width=8, height=6)

saveRDS(new_df_mammals, "../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")

ggplot(model_df_mammals) + geom_boxplot(aes(x=Diet, y=max_dis)) +geom_point(aes(x=Diet, y=max_dis), color="red")


emperical_disp<-rbind(data.frame(clade="Bird", max_dis=model_df_birds$max_dis), 
                      data.frame(clade="Mammal", max_dis=model_df_mammals$max_dis))

p<-ggplot(emperical_disp)+geom_boxplot(aes(x=clade, y=max_dis))
p

